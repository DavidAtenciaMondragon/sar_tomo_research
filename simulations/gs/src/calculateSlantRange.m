function [R1, R2, theta_inc_deg] = calculateSlantRange(Px,Py,Pz,x0,y0,z0,c,n,threshold)

u = sqrt((Px-x0).^2 + (Py-y0).^2);
u1  = u./(1-z0./(n*Pz));
u2 = u - u1;

error = 1;
while error > threshold
    R1 = sqrt(u1.^2+Pz.^2);
    R2 = sqrt(u2.^2+z0.^2);
    dR = abs(R2-R1.*n.*u2./u1);
    error = 2*n*dR/c;
    u1  = u./(1+R2./(n*R1));
    u2 = u - u1;
end

% --- NUEVO CÁLCULO DEL ÁNGULO DE INCIDENCIA ---
% Una vez que el bucle converge, R1 y Pz tienen los valores correctos.
% El ángulo de incidencia (theta) se calcula con respecto a la normal
% (la línea vertical). En el triángulo rectángulo formado por R1 (hipotenusa)
% y Pz (cateto adyacente), se cumple que: cos(theta) = Pz / R1.
% Usamos acosd para obtener el ángulo directamente en grados.
theta_inc_deg = acos(Pz ./ R1);

end