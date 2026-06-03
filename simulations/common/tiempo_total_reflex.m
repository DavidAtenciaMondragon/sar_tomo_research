function T = tiempo_total_reflex(A, B, Pxy, interfaz_obj, v)
    z_interp = interfaz_obj(Pxy(2), Pxy(1));
    if isnan(z_interp), T = inf; return; end
    P_surf = [Pxy(1), Pxy(2), z_interp];
    T = norm(A - P_surf) / v + norm(P_surf - B) / v;
end