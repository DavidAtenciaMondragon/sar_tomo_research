function [Gamma_TM, Tau_TM, R_TM, T_TM] = multislabTMcoef(f, theta_in, n_in, n_slabs, d_slabs, n_out)
% MULTISLABTMCOEF  Coeficientes de reflexion/transmision TM para un sistema
% multicapa, generalizado para medios de entrada y salida distintos.
%
% Implementa el formalismo de matriz ABCD descrito en models/multiple_slab_tm.md,
% pero sin la simplificacion Z_in = Z_out de su Seccion 6 (valida solo para
% entrada/salida en vacio). Aqui el medio de salida (n_out) puede ser distinto
% del medio de entrada (n_in), como ocurre cuando el target queda embebido en
% una capa de suelo intermedia (Seccion 11 del documento).
%
% SYNOPSIS:
%   [Gamma_TM, Tau_TM, R_TM, T_TM] = multislabTMcoef(f, theta_in, n_in, n_slabs, d_slabs, n_out)
%
% INPUT:
%   f        - Frecuencia de la onda [Hz] (escalar)
%   theta_in - Angulo de incidencia en el medio n_in [rad] [1 x N]
%   n_in     - Indice de refraccion del medio de entrada (escalar)
%   n_slabs  - Indices de refraccion de las capas intermedias [1 x M] (puede ser [])
%   d_slabs  - Espesores de las capas intermedias [m] [1 x M] (puede ser [])
%   n_out    - Indice de refraccion del medio de salida, semi-infinito (escalar)
%
% OUTPUT:
%   Gamma_TM - Coeficiente de reflexion complejo [1 x N]
%   Tau_TM   - Coeficiente de transmision complejo [1 x N]
%   R_TM     - Reflectancia de potencia |Gamma_TM|^2 [1 x N]
%   T_TM     - Transmitancia de potencia [1 x N], R_TM + T_TM = 1

c    = physconst('lightspeed');
eps0 = 8.8541878128e-12;
mu0  = 4*pi*1e-7;
eta0 = sqrt(mu0/eps0);

k0 = 2*pi*f/c;

N_slabs = numel(n_slabs);
n_seq   = [n_in, n_slabs(:)', n_out];

N = numel(theta_in);
Gamma_TM = zeros(1,N);
Tau_TM   = zeros(1,N);
R_TM     = zeros(1,N);
T_TM     = zeros(1,N);

for idx = 1:N

    % --- Ley de Snell generalizada (Seccion 2) ---
    thetas = zeros(1, N_slabs+2);
    thetas(1) = theta_in(idx);
    TIR = false;
    for i = 1:(N_slabs+1)
        arg = n_seq(i)/n_seq(i+1) * sin(thetas(i));
        if abs(arg) > 1
            TIR = true;
            break;
        end
        thetas(i+1) = asin(arg);
    end

    if TIR
        Gamma_TM(idx) = 1;
        Tau_TM(idx)   = 0;
        R_TM(idx)     = 1;
        T_TM(idx)     = 0;
        continue;
    end

    % --- Matriz ABCD de la cascada de slabs (Secciones 3-5) ---
    M_tot = eye(2);
    for i = 1:N_slabs
        n_i      = n_slabs(i);
        theta_i1 = thetas(i+1);
        eta_i    = eta0/n_i;

        k_z  = k0 * n_i * cos(theta_i1);
        Z_TM = eta_i * cos(theta_i1);

        M_i = [cos(k_z*d_slabs(i)),      1j*Z_TM*sin(k_z*d_slabs(i));
               1j/Z_TM*sin(k_z*d_slabs(i)), cos(k_z*d_slabs(i))];

        M_tot = M_tot * M_i;
    end

    A = M_tot(1,1); B = M_tot(1,2);
    C = M_tot(2,1); D = M_tot(2,2);

    % --- Impedancias TM de entrada y salida (generalizacion Seccion 6/11) ---
    Z_in  = (eta0/n_in)  * cos(thetas(1));
    Z_out = (eta0/n_out) * cos(thetas(end));

    % --- Coeficientes de reflexion y transmision (Secciones 7-8) ---
    num = A*Z_out + B - Z_in*(C*Z_out + D);
    den = A*Z_out + B + Z_in*(C*Z_out + D);
    Gamma_TM(idx) = num/den;

    Tau_TM(idx) = 2*Z_out / (A*Z_in + B + C*Z_in*Z_out + D*Z_out);

    R_TM(idx) = abs(Gamma_TM(idx))^2;
    T_TM(idx) = abs(Tau_TM(idx))^2 * (Z_in/Z_out);
end

end
