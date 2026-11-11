"""Análise do imageamento mono+biestático em solos estratificados (Capítulo 4 da tese).

Reutiliza a hélice e os parâmetros do motor analise_1610 (params/common.json) e calcula:
  (1) traçado de raios exato em N camadas (invariante de Snell, bisseção em p);
  (2) PSF exata (retroprojeção sem ruído, soma coerente em pulsos e frequências) de um
      alvo pontual, para o modelo de focalização correto e para um modelo homogêneo
      equivalente, nos canais monoestático, biestático (90 graus) e na fusão;
  (3) larguras a -3 dB frente à forma fechada do Capítulo 3 com n -> n_L;
  (4) transmitância TM/TE de pilhas com perdas (matriz característica) e sua variação de
      amplitude e de fase ao longo da abertura.

Uso:  python -I simulations/analise_multicamada/analise_multicamada.py
Saídas: simulations/io/analise_multicamada/*.csv e doc/figures/multicamada/*.pdf
"""
import json
import os
import sys

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, "..", ".."))
sys.path.insert(0, os.path.join(RAIZ, "simulations", "analise_1610"))
from bisar.geometry import tx_trajectory, rx_from_tx  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT_IO = os.path.join(RAIZ, "simulations", "io", "analise_multicamada")
OUT_FIG = os.path.join(RAIZ, "doc", "figures", "multicamada")
os.makedirs(OUT_IO, exist_ok=True)
os.makedirs(OUT_FIG, exist_ok=True)

with open(os.path.join(RAIZ, "simulations", "analise_1610", "params", "common.json"), encoding="utf-8") as fh:
    P = json.load(fh)
C = P["c"]
FC = P["fc"]
B = P["B"]
NF = P["Nf"]
FREQS = FC + B * (np.arange(NF) / (NF - 1) - 0.5)
LAMBDA0 = C / FC
XJ0 = 1.1264
CW_RECT = 0.886
ALVO = np.array([0.0, 0.0, -2.0])

# ── Cenários: lista de camadas (eps' , eps'', espessura); a última é semi-infinita ──
CENARIOS = {
    "H0": {"rot": "homogêneo, $\\varepsilon=3,99$", "camadas": [(3.99, 0.40, None)]},
    "M1": {"rot": "seca 1 m ($\\varepsilon=3$) sobre $\\varepsilon=3,99$",
           "camadas": [(3.0, 0.30, 1.0), (3.99, 0.40, None)]},
    "M2": {"rot": "úmida 0,5 m ($\\varepsilon=6,5$) sobre $\\varepsilon=3,99$",
           "camadas": [(6.5, 0.90, 0.5), (3.99, 0.40, None)]},
    "M3": {"rot": "$\\varepsilon=3,99$ 1 m sobre úmida ($\\varepsilon=6,5$)",
           "camadas": [(3.99, 0.40, 1.0), (6.5, 0.90, None)]},
}


def estrutura(camadas):
    """Índices reais, topo e base (profundidades positivas) de cada camada."""
    n = np.array([np.sqrt(c[0]) for c in camadas])
    topo, base, z = [], [], 0.0
    for c in camadas:
        topo.append(z)
        z = np.inf if c[2] is None else z + c[2]
        base.append(z)
    return n, np.array(topo), np.array(base)


def espessuras(prof, topo, base):
    """Espessura atravessada em cada camada por um raio que desce até a profundidade prof."""
    prof = np.asarray(prof, float)[..., None]
    return np.clip(np.minimum(prof, base) - topo, 0.0, None)


def tracar(D, h, esp, n, iters=64):
    """Resolve D = h p/sqrt(1-p^2) + sum_i esp_i p/sqrt(n_i^2-p^2) para p (bisseção).
    D, h: (...); esp: (..., L). Devolve p, caminho óptico e comprimento geométrico."""
    usadas = esp > 0
    pmax = np.min(np.where(usadas, n, np.inf), axis=-1)
    pmax = np.minimum(pmax, 1.0) * (1 - 1e-12)
    lo = np.zeros_like(D)
    hi = pmax.copy()
    for _ in range(iters):
        p = 0.5 * (lo + hi)
        x = h * p / np.sqrt(1 - p**2) + np.sum(esp * p[..., None] / np.sqrt(n**2 - p[..., None] ** 2), axis=-1)
        maior = x > D
        hi = np.where(maior, p, hi)
        lo = np.where(maior, lo, p)
    p = 0.5 * (lo + hi)
    cos_ar = np.sqrt(1 - p**2)
    cos_i = np.sqrt(n**2 - p[..., None] ** 2) / n
    optico = h / cos_ar + np.sum(n * esp / cos_i, axis=-1)
    return p, optico


def retardo(plat, vox, n, topo, base):
    """Retardo óptico de ida (plataforma->voxel). plat (M,3), vox (V,3) -> (M,V)."""
    D = np.hypot(plat[:, None, 0] - vox[None, :, 0], plat[:, None, 1] - vox[None, :, 1])
    h = np.broadcast_to(plat[:, None, 2], D.shape)
    esp = np.broadcast_to(espessuras(-vox[:, 2], topo, base)[None], D.shape + (len(n),))
    _, opt = tracar(D, h, esp, n)
    return opt


def psf(tx, rx, verdade, modelo, vox):
    """Imagem retroprojetada sem ruído de um alvo em ALVO: verdade e modelo são (n,topo,base)."""
    tau_v = (retardo(tx, ALVO[None], *verdade) + retardo(rx, ALVO[None], *verdade)) / C  # (M,1)
    I = np.zeros(len(vox), complex)
    for a in range(0, len(vox), 100):
        v = vox[a:a + 100]
        tau_m = (retardo(tx, v, *modelo) + retardo(rx, v, *modelo)) / C  # (M,V)
        d = tau_v - tau_m
        I[a:a + 100] = np.exp(2j * np.pi * FREQS[None, None, :] * d[..., None]).sum(axis=(0, 2))
    return I / (tx.shape[0] * NF)


def largura(eixo, I):
    """Largura a -3 dB (amplitude) com interpolação linear em dB; pico e sua posição."""
    a = np.abs(I)
    k = int(np.argmax(a))
    db = 20 * np.log10(a / a[k] + 1e-30)
    def cruza(passo):
        j = k
        while 0 <= j + passo < len(a) and db[j + passo] > -3.0103:
            j += passo
        if not 0 <= j + passo < len(a):
            return np.nan
        x0, x1, y0, y1 = eixo[j], eixo[j + passo], db[j], db[j + passo]
        return x0 + (-3.0103 - y0) * (x1 - x0) / (y1 - y0)
    return cruza(1) - cruza(-1), a[k], eixo[k]


def fechada(tx, rx, verdade):
    """Forma fechada do Capítulo 3 com n -> n_L e incidências no ar dadas pelo traçado em camadas."""
    n, topo, base = verdade
    esp = espessuras(-ALVO[2], topo, base)
    nL = n[np.nonzero(esp)[0][-1]]
    out = []
    for plat in (tx, rx):
        D = np.hypot(plat[:, 0] - ALVO[0], plat[:, 1] - ALVO[1])
        p, _ = tracar(D, plat[:, 2], np.broadcast_to(esp, D.shape + esp.shape), n)
        az = np.arctan2(plat[:, 1] - ALVO[1], plat[:, 0] - ALVO[0])
        out.append((p, az))
    (sT, aT), (sR, aR) = out
    S = 0.5 * np.sqrt(sT**2 + sR**2 + 2 * sT * sR * np.cos(aR - aT))
    q = np.sqrt(nL**2 - sT**2) + np.sqrt(nL**2 - sR**2)
    dxy = XJ0 * LAMBDA0 / (2 * np.pi * S.mean())
    dz = CW_RECT * C / np.sqrt((B * q.mean()) ** 2 + 12 * FC**2 * q.var())
    return dxy, dz, nL, np.degrees(np.arcsin(sT)).min(), np.degrees(np.arcsin(sT)).max()


def n_equivalente(n, topo, base, prof):
    """Índice homogêneo que reproduz o retardo vertical até a profundidade do alvo."""
    esp = espessuras(prof, topo, base)
    return float(np.sum(n * esp) / prof)


# ── Matriz característica (Born & Wolf) para TE e TM, com perdas ───────────────────
def coef_pilha(theta0, camadas, f, pol):
    """Coeficientes de reflexão r e de transmissão t (campo) do ar para a última camada,
    e transmitância de potência T. Convenção e^{-i w t}: n = n' + i n'' (módulo idêntico
    ao da convenção e^{+j w t} da tese)."""
    k0 = 2 * np.pi * f / C
    s0 = np.sin(theta0)
    ns = [np.sqrt(complex(e1, e2)) for e1, e2, _ in camadas]
    def adm(nc):
        cz = np.sqrt(nc**2 - s0**2 + 0j)
        if np.imag(cz) < 0:
            cz = -cz
        return (cz if pol == "TE" else cz / nc**2), cz  # TE: n cos; TM: cos/n (= n cos/n^2)
    Y0, _ = adm(1.0 + 0j)
    m11, m12, m21, m22 = 1 + 0j, 0j, 0j, 1 + 0j
    for (e1, e2, d), nc in zip(camadas[:-1], ns[:-1]):
        Y, cz = adm(nc)
        beta = k0 * cz * d
        a11, a12, a21, a22 = np.cos(beta), -1j * np.sin(beta) / Y, -1j * Y * np.sin(beta), np.cos(beta)
        m11, m12, m21, m22 = m11 * a11 + m12 * a21, m11 * a12 + m12 * a22, m21 * a11 + m22 * a21, m21 * a12 + m22 * a22
    YL, _ = adm(ns[-1])
    den = (m11 + m12 * YL) * Y0 + (m21 + m22 * YL)
    r = ((m11 + m12 * YL) * Y0 - (m21 + m22 * YL)) / den
    t = 2 * Y0 / den
    T = np.real(YL) / np.real(Y0) * np.abs(t) ** 2
    return r, t, T


def verifica_coef():
    """Verificações: Fresnel (interface única) e conservação R+T=1 com camadas sem perdas."""
    th = np.radians(np.linspace(0, 85, 18))
    erro_cons, erro_fres = 0.0, 0.0
    for pol in ("TE", "TM"):
        for t0 in th:
            r, _, T = coef_pilha(t0, [(2.0, 0.0, 0.3), (6.0, 0.0, 0.17), (3.99, 0.4, None)], FC, pol)
            erro_cons = max(erro_cons, abs(abs(r) ** 2 + T - 1))
            r1, _, _ = coef_pilha(t0, [(3.99, 0.0, None)], FC, pol)
            eps = 3.99
            cz = np.sqrt(eps - np.sin(t0) ** 2)
            rf = ((np.cos(t0) - cz) / (np.cos(t0) + cz)) if pol == "TE" else ((eps * np.cos(t0) - cz) / (eps * np.cos(t0) + cz))
            erro_fres = max(erro_fres, abs(abs(r1) - abs(rf)))
    return erro_cons, erro_fres


def principal():
    t, tx = tx_trajectory(P)
    rx = rx_from_tx(tx, 90.0)
    ex = np.arange(-0.6, 0.6 + 1e-9, 0.002)
    ez = np.arange(-2.5, 2.5 + 1e-9, 0.005)
    linhas = []
    cortes = {}
    for chave, cen in CENARIOS.items():
        verdade = estrutura(cen["camadas"])
        modelos = {"correto": verdade}
        if chave != "H0":
            neq = n_equivalente(*verdade, -ALVO[2])
            modelos["homogêneo eq."] = (np.array([neq]), np.array([0.0]), np.array([np.inf]))
            modelos["camada superior"] = (verdade[0][:1], np.array([0.0]), np.array([np.inf]))
        for nome_mod, modelo in modelos.items():
            res = {}
            for canal, (a, b) in {"mono": (tx, tx), "bi90": (tx, rx)}.items():
                vz = np.stack([np.zeros_like(ez), np.zeros_like(ez), ALVO[2] + ez], 1)
                Iz = psf(a, b, verdade, modelo, vz)
                res[canal] = {"Iz": Iz}
            Iz_f = 0.5 * (res["mono"]["Iz"] + res["bi90"]["Iz"])
            res["fusao"] = {"Iz": Iz_f}
            for canal in ("mono", "bi90", "fusao"):
                dz, pico_z, zpk = largura(ez, res[canal]["Iz"])
                res[canal].update(dz=dz, pico=pico_z, zpk=zpk)
            # corte em x na profundidade do pico de cada canal
            for canal, (a, b) in {"mono": (tx, tx), "bi90": (tx, rx)}.items():
                zc = ALVO[2] + res[canal]["zpk"]
                vx = np.stack([ex, np.zeros_like(ex), np.full_like(ex, zc)], 1)
                res[canal]["Ix"] = psf(a, b, verdade, modelo, vx)
            zc = ALVO[2] + res["fusao"]["zpk"]
            vx = np.stack([ex, np.zeros_like(ex), np.full_like(ex, zc)], 1)
            res["fusao"]["Ix"] = 0.5 * (psf(tx, tx, verdade, modelo, vx) + psf(tx, rx, verdade, modelo, vx))
            for canal in ("mono", "bi90", "fusao"):
                dx, pico_x, _ = largura(ex, res[canal]["Ix"])
                res[canal]["dx"] = dx
                res[canal]["pico_x"] = pico_x
            fm = fechada(tx, tx, verdade)
            fb = fechada(tx, rx, verdade)
            for canal in ("mono", "bi90", "fusao"):
                f = fm if canal == "mono" else fb
                linhas.append({
                    "cenario": chave, "modelo": nome_mod, "canal": canal,
                    "n_L": round(fm[2], 4),
                    "psi_min": round(fm[3], 2), "psi_max": round(fm[4], 2),
                    "dx_m": res[canal]["dx"], "dz_m": res[canal]["dz"],
                    "dxy_fech": f[0] if canal != "fusao" else np.nan,
                    "dz_fech": f[1] if canal != "fusao" else np.nan,
                    "pico_rel_dB": 20 * np.log10(max(res[canal]["pico"], res[canal]["pico_x"])),
                    "desloc_z_m": res[canal]["zpk"],
                })
            cortes[(chave, nome_mod)] = res
            print(chave, nome_mod, {c: (round(res[c]["dx"], 4), round(res[c]["dz"], 4), round(res[c]["zpk"], 3)) for c in res})

    # Normaliza o pico de cada caso pelo pico do modelo correto do mesmo cenário/canal
    for ln in linhas:
        ref = [r for r in linhas if r["cenario"] == ln["cenario"] and r["canal"] == ln["canal"] and r["modelo"] == "correto"][0]
        ln["perda_pico_dB"] = ln["pico_rel_dB"] - ref["pico_rel_dB"]

    import csv
    campos = ["cenario", "modelo", "canal", "n_L", "psi_min", "psi_max", "dx_m", "dxy_fech", "dz_m", "dz_fech",
              "desloc_z_m", "perda_pico_dB"]
    with open(os.path.join(OUT_IO, "psf_multicamada.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for ln in linhas:
            w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in ln.items()})

    # ── Transmitância ──────────────────────────────────────────────────────────────
    th = np.radians(np.linspace(0, 85, 341))
    tab_T = []
    fig, eixos = plt.subplots(1, 2, figsize=(10.5, 3.9), sharey=True)
    psi_tx = None
    for chave, cen in CENARIOS.items():
        verdade = estrutura(cen["camadas"])
        esp = espessuras(-ALVO[2], verdade[1], verdade[2])
        D = np.hypot(tx[:, 0], tx[:, 1])
        p, _ = tracar(D, tx[:, 2], np.broadcast_to(esp, D.shape + esp.shape), verdade[0])
        psi = np.arcsin(p)
        if chave == "H0":
            psi_tx = np.degrees(psi)
        # Camadas finitas sem perdas: a absorção nelas é contada à parte, como atenuação do
        # percurso (como no modelo de sinal da tese); o substrato conserva as perdas.
        camadas_alvo = [(e1, 0.0, d) for e1, _, d in cen["camadas"][:-1]] + [cen["camadas"][-1]]
        k0 = 2 * np.pi * FC / C
        for ax, pol in zip(eixos, ("TM", "TE")):
            T = np.array([coef_pilha(a, camadas_alvo, FC, pol)[2] for a in th])
            ax.plot(np.degrees(th), 10 * np.log10(T), label=cen["rot"])
            Th = np.array([coef_pilha(a, camadas_alvo, FC, pol)[2] for a in psi])
            # fase excedente do coeficiente de transmissão: fase de t menos a fase de propagação
            # vertical k0*sum(n_i cos(theta_i) d_i), já contida no retardo do raio
            def fase_exc(a):
                t_ = coef_pilha(a, camadas_alvo, FC, pol)[1]
                prop = sum(k0 * np.sqrt(e1 - np.sin(a) ** 2) * d for e1, _, d in camadas_alvo[:-1])
                return np.angle(t_ * np.exp(-1j * prop))
            ph = np.array([fase_exc(a) for a in psi])
            ph2 = np.degrees(np.unwrap(2 * ph))
            T2 = 2 * 10 * np.log10(Th)
            tab_T.append({"cenario": chave, "pol": pol,
                          "T_ida_volta_med_dB": float(np.mean(T2)),
                          "var_T_dB": float(T2.max() - T2.min()),
                          "var_fase_graus": float(ph2.max() - ph2.min())})
    for ax, pol in zip(eixos, ("TM", "TE")):
        ax.axvspan(psi_tx.min(), psi_tx.max(), color="0.85", zorder=0)
        ax.set_xlabel("Incidência no ar [°]")
        ax.set_title(f"Polarização {pol}")
        ax.set_ylim(-6, 0.3)
        ax.grid(alpha=0.3)
    eixos[0].set_ylabel("Transmitância de ida [dB]")
    eixos[1].legend(fontsize=7.5, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIG, "mc_transmitancia.pdf"))
    plt.close(fig)
    with open(os.path.join(OUT_IO, "transmitancia_multicamada.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tab_T[0]))
        w.writeheader()
        for ln in tab_T:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in ln.items()})

    # ── Figura dos cortes da PSF ───────────────────────────────────────────────────
    fig, eixos = plt.subplots(2, 3, figsize=(11, 6.2))
    for j, chave in enumerate(("M1", "M2", "M3")):
        for nome_mod, est in (("correto", "-"), ("homogêneo eq.", "--"), ("camada superior", ":")):
            if (chave, nome_mod) not in cortes:
                continue
            res = cortes[(chave, nome_mod)]
            ref = cortes[(chave, "correto")]
            for canal, cor in (("mono", "C0"), ("fusao", "C2")):
                a = np.abs(res[canal]["Iz"]) / np.abs(ref[canal]["Iz"]).max()
                eixos[0, j].plot(ez, 20 * np.log10(a + 1e-12), est, color=cor,
                                 label=f"{'mono' if canal == 'mono' else 'fusão'}, {nome_mod}")
                a = np.abs(res[canal]["Ix"]) / np.abs(ref[canal]["Ix"]).max()
                eixos[1, j].plot(ex, 20 * np.log10(a + 1e-12), est, color=cor)
        eixos[0, j].set_title(CENARIOS[chave]["rot"], fontsize=9)
        eixos[0, j].set_xlabel("$z-z_P$ [m]")
        eixos[1, j].set_xlabel("$x$ [m]")
        for i in range(2):
            eixos[i, j].set_ylim(-30, 1)
            eixos[i, j].grid(alpha=0.3)
    eixos[0, 0].set_ylabel("Corte vertical [dB]")
    eixos[1, 0].set_ylabel("Corte horizontal [dB]")
    eixos[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_FIG, "mc_psf_cortes.pdf"))
    plt.close(fig)

    ec, ef = verifica_coef()
    print(f"Verificação: |R+T-1| max = {ec:.2e}; |r| vs Fresnel max = {ef:.2e}")
    with open(os.path.join(OUT_IO, "verificacao.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"max|R+T-1| (pilha com camadas sem perdas e substrato com perdas) = {ec:.3e}\n")
        fh.write(f"max||r|-|r_Fresnel|| (interface única) = {ef:.3e}\n")
    for ln in tab_T:
        print(ln)


if __name__ == "__main__":
    principal()
