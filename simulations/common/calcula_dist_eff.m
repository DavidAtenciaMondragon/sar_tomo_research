function dist = calcula_dist_eff(A,B)

dist = sqrt(sum(A.^2, 2) - 2 * (mtimes(A,B')) + sum(B.^2, 2)');

end