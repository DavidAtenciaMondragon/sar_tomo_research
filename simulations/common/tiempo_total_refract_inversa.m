function T = tiempo_total_refract_inversa(C, B, Rxy, interfaz_obj, v2, v1)
    z_interp = interfaz_obj(Rxy(2), Rxy(1));
    if isnan(z_interp), T = inf; return; end
    R_surf = [Rxy(1), Rxy(2), z_interp];
    T = norm(C - R_surf) / v2 + norm(R_surf - B) / v1;
end