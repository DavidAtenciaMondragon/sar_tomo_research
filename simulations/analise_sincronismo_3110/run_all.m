function run_all()
% RUN_ALL  Ejecuta todos los experimentos en orden (~30 min en un PC de escritorio).
%   exp1  tolerancias (precision NECESARIA)
%   exp2  demostracion de extremo a extremo (una realizacion, figuras)
%   exp3  precision ALCANZADA vs SNR del directo y oscilador (Monte Carlo)
%   exp4  pulsos perdidos / patrones de escucha y puenteo de fase
%   exp5  observabilidad de la asociacion de pulsos (stagger, J)
%   exp6  acoplamiento con el error de navegacion del Tx
%   fig_informe  figuras adicionales de informe_sincronismo.tex
here = fileparts(mfilename('fullpath'));
for e = {'exp1_tolerancias', 'exp2_cadena_completa', 'exp3_snr_oscilador', ...
         'exp4_pulsos_perdidos', 'exp5_asociacion', 'exp6_navegacion', 'fig_informe'}
    fprintf('\n================ %s ================\n', e{1});
    t0 = tic;
    ejecuta(fullfile(here, [e{1} '.m']));
    fprintf('%s: %.0f s\n', e{1}, toc(t0));
end
end

function ejecuta(script)
% Espacio de trabajo propio: el "clear" de cada experimento no afecta al bucle
run(script);
end
