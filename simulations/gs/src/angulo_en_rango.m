function dentro = angulo_en_rango(angulo, minimo, maximo)

% Normalizar todos los ángulos al rango [0, 360)
angulo = mod(angulo, 360);
minimo = mod(minimo, 360);
maximo = mod(maximo, 360);

if minimo <= maximo
    dentro = (angulo >= minimo) && (angulo <= maximo);
else
    % Intervalo cruzando el 0°, como por ejemplo de 300° a 60°
    dentro = (angulo >= minimo) || (angulo <= maximo);
end

end