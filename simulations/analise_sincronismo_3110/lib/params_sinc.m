function P = params_sinc()
% PARAMS_SINC  Parametros por defecto de la simulacion de sincronismo BiSAR.
%
% Convencion de tiempos
%   * El reloj del Radar 1 (Tx/Rx monostatico, dron) se toma como referencia
%     ("tiempo verdadero" t). Su error absoluto frente al tiempo fisico solo
%     escala fc y la PRF en ~ppm y no afecta al problema relativo.
%   * El Radar 2 (receptor pasivo) tiene reloj local t2 = t + delta(t),
%     delta(t) = delta0 + eps*t + gamma*t^2/2 + w(t) + jitter, y su LO deriva
%     del mismo oscilador -> fase de demodulacion 2*pi*fc*delta(t) + phi0.
%
% Banda, PRF y alturas heredadas de analise_1610/params/common.json
% (fc = 425 MHz, B = 50 MHz, PRF = 400 Hz, Tp = 1.3 us).

P.c   = 299792458;
P.fc  = 425e6;
P.B   = 50e6;
P.Tp  = 1.3e-6;
P.fs  = 100e6;               % muestreo complejo en banda base (ambos radares)
P.PRF = 400;
P.PRT0 = 1/P.PRF;

% --- Geometria: vuelo circular de ambos radares ------------------------
% Los dos drones vuelan el MISMO circulo (centro en el origen, radio R, altura
% h) a la misma velocidad, con el Radar 2 adelantado dphi en azimut. La
% distancia Tx-Rx2 es constante: 2R sin(dphi/2) = 70.7 m.
P.radio = 50;                % m
P.alt   = 100;               % m
P.v     = 10;                % m/s (ambos)
P.th0   = -pi/2;             % azimut del Tx en t = 0
P.dphi_rx = pi/2;            % Rx2 adelantado 90 grados (sentido antihorario)
P.arco_deg = 360;            % arco de la apertura escuchada (una vuelta)
% Instantes de emision t_n = n*PRT0 + s_n con s_n ~ U(-S,S) conocido por el
% Tx (log). Es la "firma" que hace observable el indice de pulso.
P.stagger_max = 5e-6;

% --- Escena: reflectores puntuales en z = 0 -----------------------------
% Dentro del circulo de vuelo. Con apertura circular completa la IRF es de
% tipo J0 (lobulo principal ~0.4 m, anillos que decaen lentamente), de modo
% que los blancos se separan >= 32 m para que sus anillos no dominen.
P.targets = [  0   0 0;
              25  20 0;
             -25 -20 0;
              25 -20 0;
             -25  20 0];

% --- SNR por pulso tras compresion en rango -----------------------------
P.SNR_d_dB = 30;             % camino directo en el Rx2
P.SNR_e_dB = 20;             % eco bistatico de cada blanco
P.SNR_m_dB = 20;             % eco monostatico de cada blanco
P.det_umbral_dB = 13;        % umbral sobre ruido para detectar el directo

% --- Escucha del Rx2 ----------------------------------------------------
P.rx.t_start = 0.40;         % instante (reloj Tx) en que el Rx2 empieza
P.rx.M_slots = round(deg2rad(P.arco_deg)*P.radio/P.v/P.PRT0);  % 12566 (31.4 s)
P.rx.obs     = [];           % mascara logica M_slots x 1 ([] = continua)

% --- Oscilador del Rx2 y offset inicial ---------------------------------
P.osc = osc_clase('TCXO');
P.delta0_rango = [5 20];     % s; el Rx2 arranca con un reloj arbitrario

P.N_tx = P.rx.M_slots + round((P.rx.t_start + 1)*P.PRF);   % pulsos en el log del Tx

% --- Error de navegacion (proceso AR(1) por eje) ------------------------
P.nav.sigma    = 0;          % m rms por eje, Tx
P.nav.sigma_rx = 0;          % m rms por eje, Rx2
P.nav.tcorr = 1.0;           % s

% --- Procesado e imagen -------------------------------------------------
P.U        = 8;              % sobremuestreo del dato comprimido para BP
P.seg_pre  = 10;             % muestras antes del directo en el segmento BP
P.seg_post = 100;            % muestras despues del directo (cubre la escena)
P.parche_semi = 3;           % m, semiancho del parche de metricas por blanco
P.parche_d    = 0.075;       % m (lobulo principal ~0.4 m)
P.ml_radio    = 0.5;         % m, radio maximo del lobulo principal (1er nulo J0: 0.43 m)
P.escena_semi = 32;          % m, imagen completa (solo figuras)
P.escena_d    = 0.1;
P.mono = false;              % generar tambien el dato monostatico

% --- Criterios de exito (imagen sincronizada vs. reloj ideal) -----------
P.crit.perdida_pico_dB = 1.0;   % |perdida de pico| maxima
P.crit.err_pos_m       = 0.05;  % desplazamiento maximo del pico (~0.1 resolucion)
P.crit.dPSLR_dB        = 3.0;
P.crit.dISLR_dB        = 3.0;
P.crit.rho_min         = 0.95;  % correlacion compleja con la imagen ideal
end
