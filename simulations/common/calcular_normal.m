function N_unit = calcular_normal(Pxy, interfaz_obj)
    x = Pxy(1); y = Pxy(2); h = 1e-5;
    dzdx = (interfaz_obj(y, x + h) - interfaz_obj(y, x - h)) / (2*h);
    dzdy = (interfaz_obj(y + h, x) - interfaz_obj(y - h, x)) / (2*h);
    N = [-dzdx, -dzdy, 1];
    N_unit = N / norm(N);
end