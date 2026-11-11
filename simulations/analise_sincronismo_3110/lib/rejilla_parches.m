function R = rejilla_parches(P)
% REJILLA_PARCHES  Parches cuadrados de metricas centrados en cada blanco (z = 0).
u = -P.parche_semi:P.parche_d:P.parche_semi;
[X, Y] = meshgrid(u, u);
K = size(P.targets, 1); n = numel(X);
G = zeros(K*n, 3);
for k = 1:K
    G((k-1)*n + (1:n), :) = [X(:) + P.targets(k,1), Y(:) + P.targets(k,2), zeros(n,1)];
end
R.G = G; R.u = u; R.n = numel(u); R.K = K; R.npix = n;
end
