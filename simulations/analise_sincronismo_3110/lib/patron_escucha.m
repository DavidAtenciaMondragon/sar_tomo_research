function obs = patron_escucha(M, tipo, par, rs)
% PATRON_ESCUCHA  Mascara de ranuras PRT en las que el Rx2 graba (1a = adquisicion).
%   'continuo'            : todas
%   'periodico', G        : una de cada G
%   'aleatorio', p        : cada ranura con probabilidad p
%   'rafagas', [Non Noff] : Non ranuras grabadas, Noff ignoradas
obs = false(M, 1);
switch tipo
    case 'continuo',  obs(:) = true;
    case 'periodico', obs(1:par:M) = true;
    case 'aleatorio', obs = rand(rs, M, 1) < par;
    case 'rafagas',   obs(mod((0:M-1).', sum(par)) < par(1)) = true;
    otherwise, error('patron desconocido');
end
obs(1) = true;
end
