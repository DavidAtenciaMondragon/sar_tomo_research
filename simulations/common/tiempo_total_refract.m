function T = tiempo_total_refract(A, C, Qxy, interfaz_obj, v1, v2)
    z_interp = interfaz_obj(Qxy(2), Qxy(1));
    if isnan(z_interp), T = inf; return; end
    Q_surf = [Qxy(1), Qxy(2), z_interp];
    T = norm(A - Q_surf) / v1 + norm(Q_surf - C) / v2;
end