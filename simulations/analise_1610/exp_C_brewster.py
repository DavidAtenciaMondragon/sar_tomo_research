"""Experimento C - Incidencia con y sin pseudo-Brewster: supresion del retorno de superficie.

Parte 1 (enlace, analitico): reflectancia TM/TE de los suelos S0..S3, errores angulares y de
          atitud (roll), mapa del sigma0 biestatico SPM (VV y HH) en funcion de (theta, dphi).
Parte 2 (imagen): SCR en el voxel del blanco = |I_blanco(r0)|^2 / E|I_superficie|^2 (vecindad
          del blanco, superficie congelada compartida por los canales mono y bi).
Parte 3 (calidad de imagen vs incidencia): PSF 3D mono/bi/fusion para las dos familias de
          trayectoria (F1 helice girada, F2 cono de incidencia constante = restriccion Brewster).
Uso: python exp_C_brewster.py
"""
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from bisar import config, em, geometry, imaging, metrics, plotting as P, runner, scene

NAME = "expC"
p = config.load(NAME)
config.save_resolved(p, NAME)
out = config.outdir(NAME)
tgt = np.array(p["blanco"])
t, Ptx0 = geometry.tx_trajectory(p)
M = len(Ptx0)
th0 = geometry.incidence_deg(Ptx0)
az0 = np.arctan2(Ptx0[:, 1], Ptx0[:, 0])
r0 = np.linalg.norm(Ptx0, axis=1)
SOILS = p["suelos"]
eps = {k: em.eps_complex(v["eps_r"]) for k, v in SOILS.items()}
thB = {k: em.pseudo_brewster_deg(e)[0] for k, e in eps.items()}
print("pseudo-Brewster:", {k: round(v, 3) for k, v in thB.items()})


def family(name, theta_c):
    return geometry.familia_C(Ptx0, name, theta_c)


# =========================================================== Parte 1: enlace
th = np.linspace(0, 89.5, 1791)
rows = []
for k, e in eps.items():
    for d in [0] + p["errores_deg"]:
        a = thB[k] + d
        rows.append({"suelo": k, "eps_r": str(e), "theta_pB": thB[k], "err_deg": d, "theta": a,
                     "Gamma_TM_dB": 10 * np.log10(em.reflectance(np.radians(a), e, "TM")),
                     "Gamma_TE_dB": 10 * np.log10(em.reflectance(np.radians(a), e, "TE")),
                     **{f"Gamma_TM_roll{r}_dB": 10 * np.log10(em.reflectance_mixed(np.radians(a), e, r))
                        for r in p["roll_list_deg"] if r > 0},
                     **{f"chi_coh_sh{int(sh*100)}cm_dB": 20 * np.log10(em.coherent_factor(
                         2 * np.pi * p["fc"] / p["c"], sh, np.radians(a))) for sh in p["rugosidad_lista"]}})
dl = pd.DataFrame(rows).drop_duplicates(["suelo", "err_deg"])
dl.to_csv(f"{out}/expC_enlace_reflectancia.csv", index=False)

# --- Fig. 7 reflectancia TM y TE
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
cols = [P.SERIES[i] for i in (0, 1, 2, 3)]
for (k, e), col in zip(eps.items(), cols):
    axs[0].plot(th, 10 * np.log10(em.reflectance(np.radians(th), e, "TM")), color=col,
                label=f"{SOILS[k]['nombre']}  $\\varepsilon_r$={e.real:g}-j{-e.imag:g}  ($\\theta_{{pB}}$={thB[k]:.2f}°)")
    axs[0].plot(th, 10 * np.log10(em.reflectance(np.radians(th), e, "TE")), color=col, ls=":")
axs[0].set_ylim(-45, 0); axs[0].set_xlim(0, 90)
axs[0].set_xlabel("incidencia $\\theta_i$ [°]"); axs[0].set_ylabel("reflectancia [dB]")
axs[0].set_title("TM (continuo) y TE (punteado)")
axs[0].legend(fontsize=6.5, loc="lower left")
axs[0].axvspan(th0.min(), th0.max(), color=P.C["grid"], alpha=0.5, lw=0)
axs[0].text(th0.mean(), -43, "helice\nbase", ha="center", fontsize=7, color=P.C["ink2"])
e1 = eps["S1"]
d = np.linspace(-6, 6, 241)
for r, col in zip(p["roll_list_deg"], [P.SERIES[i] for i in (0, 1, 2, 3)]):
    axs[1].plot(d, 10 * np.log10(em.reflectance_mixed(np.radians(thB["S1"] + d), e1, r)), color=col,
                label=f"roll {r}°")
axs[1].set_xlabel("error angular respecto a $\\theta_{pB}$ [°]")
axs[1].set_ylabel("reflectancia TM efectiva [dB]")
axs[1].set_title("S1: robustez del minimo (error de incidencia y de atitud)")
for e_ in [0.5, 1, 2, 5]:
    axs[1].axvline(e_, color=P.C["grid"], lw=0.6); axs[1].axvline(-e_, color=P.C["grid"], lw=0.6)
axs[1].legend(fontsize=7)
P.save(fig, "C_reflectancia_basico", out)

# --- mapa sigma0 SPM biestatico (Rx a la misma elevacion que el Tx)
k0 = 2 * np.pi * p["fc"] / p["c"]
ru = p["rugosidad"]
TH, DP = np.meshgrid(np.linspace(20, 80, 241), np.linspace(0, 180, 361), indexing="ij")
phis = np.radians(180 - DP)               # phi_s = 0 especular (dphi = 180), pi retro (dphi = 0)
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
for a, pol in zip(axs, ["TM", "TE"]):
    s0 = em.spm_sigma0(np.radians(TH), np.radians(TH), phis, e1, pol, k0, ru["sigma_h"], ru["corr_len"])
    im = a.pcolormesh(DP, TH, 10 * np.log10(s0), cmap=P.SEQ, vmin=-60, vmax=-10, shading="auto", rasterized=True)
    a.set_xlabel("$\\Delta\\phi$ Tx-Rx [°]"); a.set_ylabel("incidencia $\\theta$ [°]")
    a.set_title(f"$\\sigma^0_{{{'VV' if pol == 'TM' else 'HH'}}}$ SPM, S1, $s$={ru['sigma_h']*100:g} cm, $\\ell$={ru['corr_len']*100:g} cm")
    a.grid(False)
    a.axhspan(th0.min(), th0.max(), color="w", alpha=0.15, lw=0)
    a.axhline(thB["S1"], color=P.C["s8"], lw=0.8, ls="--")
fig.colorbar(im, ax=axs, label="$\\sigma^0$ [dB]")
# locus del nulo VV (eps real aprox.): cos(phi_s) = eps sin^2 / (eps - sin^2)
tt = np.radians(np.linspace(20, thB["S1"], 200))
er = e1.real
cphi = er * np.sin(tt) ** 2 / (er - np.sin(tt) ** 2)
axs[0].plot(180 - np.degrees(np.arccos(np.clip(cphi, -1, 1))), np.degrees(tt), color=P.C["s8"], lw=1.2,
            label="nulo VV ('Brewster biestatico')")
axs[0].legend(fontsize=7, loc="lower left")
P.save(fig, "C_sigma0_mapa", out)
dphi_null = lambda theta: 180 - np.degrees(np.arccos(np.clip(er * np.sin(np.radians(theta)) ** 2 /
                                                             (er - np.sin(np.radians(theta)) ** 2), -1, 1)))
print("dphi del nulo VV para theta = 53.39: %.1f°" % dphi_null(53.39))

# =========================================================== Parte 2: SCR en imagen
vn = p["vecindad_scr"]
gN = imaging.Grid.around(tgt, vn["h"], vn["dx"], -vn["hz"], vn["hz"], vn["dz"])
i0 = np.argmin(np.linalg.norm(gN.vox - tgt, axis=1))
rough = dict(ru)
fac = scene.surface_facets(rough["patch_half"], rough["facet"], rough["taper"])
R = p["realizaciones"]
sh_ref = rough["sigma_h"]


def eval_channel(sysm, Ptx, Prx, seeds):
    s_t, _ = scene.target_echo(sysm, Ptx, Prx, tgt)
    It = runner.bp(sysm, Ptx, Prx, s_t, gN.vox, gate=True)[i0]
    out_ = {"Pt_W": np.mean(np.abs(s_t) ** 2), "It": It}
    for sh in p["rugosidad_lista"]:
        s_s, _ = scene.specular_echo(sysm, Ptx, Prx, sh)
        out_[f"Ps_W_{sh}"] = np.mean(np.abs(s_s) ** 2)
        out_[f"Is_{sh}"] = runner.bp(sysm, Ptx, Prx, s_s, gN.vox, gate=True)
    Id, Pd = [], []
    for sd in seeds:
        s_d = scene.diffuse_echo(sysm, Ptx, Prx, rough, np.random.default_rng(sd), fac)
        Pd.append(np.mean(np.abs(s_d) ** 2))
        Id.append(runner.bp(sysm, Ptx, Prx, s_d, gN.vox, gate=True))
    out_["Pd_W_ref"] = np.mean(Pd)
    out_["Id_ref"] = np.array(Id)              # (R, V) para sigma_h de referencia
    return out_


def scr_rows(ch_res, meta):
    """Combina blanco + superficie para cada rugosidad y canal (incluida la fusion)."""
    rows_ = []
    chans = dict(ch_res)
    for dphi in [k for k in ch_res if k != "mono"]:
        chans[f"fus{dphi}"] = None
    for sh in p["rugosidad_lista"]:
        sc = sh / sh_ref                                   # sigma0 ~ sigma_h^2 -> amplitud ~ sigma_h
        for ch in chans:
            if isinstance(ch, str) and ch.startswith("fus"):
                b = ch_res[int(ch[3:])]
                m = ch_res["mono"]
                It = 0.5 * (m["It"] + b["It"])
                Ic = 0.5 * ((m[f"Is_{sh}"] + b[f"Is_{sh}"])[None] + sc * (m["Id_ref"] + b["Id_ref"]))
                Pt = Ps = Pd = np.nan
                dphi = int(ch[3:])
            else:
                r = ch_res[ch]
                It = r["It"]
                Ic = r[f"Is_{sh}"][None] + sc * r["Id_ref"]
                Pt, Ps, Pd = r["Pt_W"], r[f"Ps_W_{sh}"], r["Pd_W_ref"] * sc ** 2
                dphi = 0 if ch == "mono" else ch
            Pc_img = np.mean(np.abs(Ic) ** 2)
            rows_.append({**meta, "sigma_h": sh, "canal": "mono" if ch == "mono" else ("fus" if str(ch).startswith("fus") else "bi"),
                          "dphi": dphi, "Pt_dBm": 10 * np.log10(Pt) + 30, "Psup_esp_dBm": 10 * np.log10(Ps + 1e-30) + 30,
                          "Psup_dif_dBm": 10 * np.log10(Pd + 1e-30) + 30,
                          "SCR_enlace_dB": 10 * np.log10(Pt / (Ps + Pd)) if np.isfinite(Pt) else np.nan,
                          "It_dB": 20 * np.log10(abs(It)), "Ic_img_dB": 10 * np.log10(Pc_img),
                          "SCR_img_dB": 10 * np.log10(abs(It) ** 2 / Pc_img)})
    return rows_


cases = []
for fam in ["F1", "F2"]:
    for thc in p["theta_fuera"]:
        cases.append(("S1", fam, thc, np.nan, "fuera"))
    for d in p["errores_deg"]:
        cases.append(("S1", fam, thB["S1"] + d, d, "pB"))
for sk in p["suelos_extra"]:
    for d in [-2, 0, 2]:
        cases.append((sk, "F2", thB[sk] + d, d, "pB"))
    cases.append((sk, "F2", 53.39, np.nan, "fuera"))
    cases.append((sk, "F1", thB[sk], 0, "pB"))
    cases.append((sk, "F1", 53.39, np.nan, "fuera"))

srows = []
T0 = time.time()
for ic, (sk, fam, thc, err, kind) in enumerate(cases):
    Ptx = family(fam, thc)
    seeds = [p["semilla"] + 1000 * ic + r for r in range(R)]
    for pol in ["TM", "TE"]:
        sysm = scene.System(p, soil=SOILS[sk]["eps_r"], pol=pol)
        res = {"mono": eval_channel(sysm, Ptx, Ptx, seeds)}
        for dphi in p["dphi_canales"]:
            res[dphi] = eval_channel(sysm, Ptx, geometry.rx_from_tx(Ptx, dphi), seeds)
        meta = {"suelo": sk, "familia": fam, "theta_c": thc, "err_deg": err, "tipo": kind, "pol": pol,
                "theta_min": geometry.incidence_deg(Ptx).min(), "theta_max": geometry.incidence_deg(Ptx).max()}
        srows += scr_rows(res, meta)
    print(f"[{ic+1}/{len(cases)}] {sk} {fam} theta={thc:.2f}  ({time.time()-T0:.0f}s)", flush=True)
# atitud (roll) en el cono de Brewster, S1, TM
for roll in p["roll_list_deg"]:
    Ptx = family("F2", thB["S1"])
    seeds = [p["semilla"] + 77 + r for r in range(R)]
    sysm = scene.System(p, soil=SOILS["S1"]["eps_r"], pol="TM", roll_deg=roll)
    res = {"mono": eval_channel(sysm, Ptx, Ptx, seeds)}
    for dphi in p["dphi_canales"]:
        res[dphi] = eval_channel(sysm, Ptx, geometry.rx_from_tx(Ptx, dphi), seeds)
    srows += scr_rows(res, {"suelo": "S1", "familia": "F2", "theta_c": thB["S1"], "err_deg": 0.0,
                            "tipo": f"roll{roll}", "pol": "TM", "theta_min": np.nan, "theta_max": np.nan})
ds = pd.DataFrame(srows)
ds.to_csv(f"{out}/expC_scr_imagen.csv", index=False)

# --- Fig. 8 SCR vs incidencia y error angular (S1, rugosidad moderada y superficie lisa)
for surf, sh in [("rugosa", sh_ref), ("lisa", 0.0)]:
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.4), sharey="row")
    for j, fam in enumerate(["F1", "F2"]):
        for i, pol in enumerate(["TM", "TE"]):
            a = axs[i, j]
            d = ds[(ds.suelo == "S1") & (ds.familia == fam) & (ds.pol == pol) & (ds.sigma_h == sh)
                   & ds.tipo.isin(["fuera", "pB"])].sort_values("theta_c")
            for (ch, dphi), col, lab in [(("mono", 0), P.C["mono"], "Mono"), (("bi", 90), P.SERIES[6], "Bi 90°"),
                                         (("bi", 135), P.SERIES[3], "Bi 135°"), (("bi", 180), P.C["bi"], "Bi 180° (Ref1)"),
                                         (("fus", 135), P.C["fus"], "Mono+Bi 135°")]:
                dd = d[(d.canal == ch) & (d.dphi == dphi)]
                a.plot(dd.theta_c, dd.SCR_img_dB, color=col, marker="o", ms=3.5, label=lab)
            a.axvline(thB["S1"], color=P.C["s8"], lw=0.8, ls="--")
            a.set_title(f"{'F1 helice girada' if fam == 'F1' else 'F2 cono de incidencia constante'} | {pol}")
            a.set_xlabel("incidencia central $\\theta_c$ [°]")
            if j == 0:
                a.set_ylabel("SCR en imagen [dB]")
    axs[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle(f"Exp. C: SCR en el voxel del blanco (S1, superficie {surf}"
                 f"{'' if sh == 0 else f', s = {sh*100:g} cm'}); linea roja: pseudo-Brewster", y=1.0)
    fig.tight_layout()
    P.save(fig, f"C_scr_incidencia_{surf}", out)

# --- SCR vs error angular (zoom) para F2, TM y TE, rugosa
fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.4), sharey=True)
for a, pol in zip(axs, ["TM", "TE"]):
    d = ds[(ds.suelo == "S1") & (ds.familia == "F2") & (ds.pol == pol) & (ds.tipo == "pB") & (ds.sigma_h == sh_ref)]
    for (ch, dphi), col, lab in [(("mono", 0), P.C["mono"], "Mono"), (("bi", 135), P.SERIES[3], "Bi 135°"),
                                 (("bi", 180), P.C["bi"], "Bi 180° (Ref1)"), (("fus", 135), P.C["fus"], "Mono+Bi 135°")]:
        dd = d[(d.canal == ch) & (d.dphi == dphi)].sort_values("err_deg")
        a.plot(dd.err_deg, dd.SCR_img_dB, color=col, marker="o", ms=3.5, label=lab)
    a.set_xlabel("error de incidencia respecto a $\\theta_{pB}$ [°]")
    a.set_title(f"F2 (cono de Brewster), S1, s = 2 cm, {pol}")
axs[0].set_ylabel("SCR en imagen [dB]")
axs[0].legend(fontsize=7)
P.save(fig, "C_scr_error_angular", out)

# =========================================================== Parte 3: costo en imagen
rg = p["rejilla_psf"]
g = imaging.Grid.around(tgt, rg["hx"], rg["dx"], rg["z_lo"], rg["z_hi"], rg["dz"])
sysm = scene.System(p)
prow = []
for fam in ["F1", "F2"]:
    for thc in sorted(set(p["theta_fuera"] + [thB["S1"]])):
        Ptx = family(fam, thc)
        s_m, _ = scene.target_echo(sysm, Ptx, Ptx, tgt)
        I_m = runner.bp(sysm, Ptx, Ptx, s_m, g.vox, gate=True)
        for dphi in [90, 135]:
            Prx = geometry.rx_from_tx(Ptx, dphi)
            s_b, _ = scene.target_echo(sysm, Ptx, Prx, tgt)
            imgs = runner.four_images(sysm, Ptx, Prx, s_m, s_b, g.vox, I_m, equal_budget=False)
            prow += runner.psf_rows(sysm, g, imgs, tgt, M, {"familia": fam, "theta_c": thc, "dphi_bi": dphi})
        print(f"PSF {fam} theta={thc:.2f}", flush=True)
dp = pd.DataFrame(prow)
dp.to_csv(f"{out}/expC_psf_vs_incidencia.csv", index=False)

fig, axs = plt.subplots(1, 4, figsize=(13, 3.2))
for fam, ls in [("F1", "-"), ("F2", "--")]:
    for ch, col in [("mono", P.C["mono"]), ("fus", P.C["fus"])]:
        d = dp[(dp.familia == fam) & (dp.canal == ch) & (dp.dphi_bi == 135)].sort_values("theta_c")
        lab = f"{'Mono' if ch == 'mono' else 'Mono+Bi 135°'} {fam}"
        axs[0].plot(d.theta_c, d.fwhm_x, color=col, ls=ls, marker="o", ms=3, label=lab)
        axs[1].plot(d.theta_c, d.fwhm_z, color=col, ls=ls, marker="o", ms=3)
        axs[2].plot(d.theta_c, d.vol3dB_m3 * 1e3, color=col, ls=ls, marker="o", ms=3)
        axs[3].plot(d.theta_c, d.PSLR_dB, color=col, ls=ls, marker="o", ms=3)
for a, yl in zip(axs, ["$\\delta_x$ (-3 dB) [m]", "$\\delta_z$ (-3 dB) [m]", "volumen -3 dB [l]", "PSLR 3D [dB]"]):
    a.axvline(thB["S1"], color=P.C["s8"], lw=0.8, ls="--")
    a.set_xlabel("$\\theta_c$ [°]"); a.set_ylabel(yl)
axs[0].legend(fontsize=6.5)
fig.suptitle("Exp. C: costo en calidad de imagen de la restriccion de incidencia (F2 = cono de pseudo-Brewster)", y=1.03)
fig.tight_layout()
P.save(fig, "C_costo_imagen", out)
print("Exp. C terminado ->", out)
