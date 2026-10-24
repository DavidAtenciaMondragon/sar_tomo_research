function valida_bp_matlab()
% VALIDA_BP_MATLAB  Ejecuta el nucleo de back-projection de proc/PROC_snell_test_script.m
% (mismo bucle, misma interpolacion y misma compensacion de fase) sobre cortes 1D finos
% que pasan por el blanco, usando el raw data ya generado por GS_snell_test_script.m
% (io/snell/raw/bistatic_sar_snell_raw_espiral_mono0.mat).  Sirve para comparar con el
% motor Python de analise_1610 (tests/t_bp_matlab.py) sin el efecto de la rejilla gruesa
% (dxy = 0.04 m, dz = 0.08 m) de gridSnell.
here = fileparts(mfilename('fullpath'));
sim  = fileparts(fileparts(here));
addpath(fullfile(sim,'proc','src')); addpath(genpath(fullfile(sim,'tools')));
suf = '_espiral_mono0';
S  = json2struct(fullfile(sim,'parametros',['system' suf '.json'])); S = S.system;
Tx = json2struct(fullfile(sim,'parametros',['radarTx' suf '.json'])); Tx = Tx.radar;
raw = load(fullfile(sim,'io','snell','raw',['bistatic_sar_snell_raw' suf '.mat']));
c = S.VelocidadeLuz; n2 = S.IndiceRefracaoSolo; lamb = c/Tx.FreqPortadora; fs = Tx.fs;
rc = raw.rootData.';
x0 = raw.x0; y0 = raw.y0; z0 = raw.z0;
cuts = struct('name', {'x','z'}, ...
    'vox', {[x0+(-0.03:0.0005:0.03)', y0+0*(-0.03:0.0005:0.03)', z0+0*(-0.03:0.0005:0.03)'], ...
            [x0+0*(-1:0.005:1)', y0+0*(-1:0.005:1)', z0+(-1:0.005:1)']});
for k = 1:numel(cuts)
    V = cuts(k).vox; out = zeros(size(V,1),1);
    for m = 1:size(V,1)
        [R1t,R2t] = slantRange(raw.PxT,raw.PyT,raw.PzT,V(m,1),V(m,2),V(m,3),n2,1e-11);
        [R1r,R2r] = slantRange(raw.PxR,raw.PyR,raw.PzR,V(m,1),V(m,2),V(m,3),n2,1e-11);
        t = (R1t + n2*R2t + n2*R2r + R1r)/c; rb = 1 + t*fs;
        phi = (2*pi/lamb)*(R1t + n2*R2t + n2*R2r + R1r);
        il = floor(rb); q = rb - il;
        idx1 = sub2ind(size(rc), 1:numel(rb), il); idx2 = sub2ind(size(rc), 1:numel(rb), il+1);
        d = (1-q).*rc(idx1) + q.*rc(idx2);
        out(m) = sum(d .* exp(1i*phi));
    end
    T = table(V(:,1),V(:,2),V(:,3),abs(out),'VariableNames',{'x','y','z','amp'});
    writetable(T, fullfile(here, sprintf('matlab_bp_corte_%s.csv', cuts(k).name)));
    a = abs(out)/max(abs(out)); ax = V(:, 1 + 2*(k==2));
    [~,ip] = max(a); w = compute3dBWidth(ax, 20*log10(a), ip, -3);
    fprintf('MATLAB corte %s: FWHM = %.5f m\n', cuts(k).name, w);
end
end
