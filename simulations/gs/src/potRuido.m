function powerNoise = potRuido(T,B)
    k = 1.38e-23; % Constante de Boltzmann
    powerNoise = k*T*B;
end