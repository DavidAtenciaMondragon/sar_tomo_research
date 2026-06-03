function [Px, Py, Pz, t] = funcao_espiral(N, R_top, R_bottom, H_high,H_low, Vt, PRF, NorthOffset)
% Essa função retorna os pontos em x, y e z de uma trajetória em espiral
% dado as variáveis de entrada:
% N = Número de voltas desejado
% R_top = Raio superior da espiral [m]
% R_bottom = Raio inferior da espiral [m]
% H_high = Altura máxima da espiral [m]
% H_low = Altura mínima da espiral [m]
% Vt = Velocidade tangencial [m/s]
% PRF = Frequência de repetição de pulso [Hz]

% Variáveis de saída ------------------------------------------------------
% Px = Coordenada x [m]
% Py = Coordenada y [m]
% Pz = Coordenada z [m]

%--------------------------------------------------------------------------

Vr = Vt*(log(R_top)-log(R_bottom))/(N*2*pi);	% Taxa de redução do raio [m/s]
T = (R_top-R_bottom)/Vr;                        % Tempo total [s]
Vd = (H_high-H_low)/T;                          % Velocidade de descida [m/s]

t = 0:1/PRF:T;                                  % Tempo [s]
theta = (Vt/Vr)*(log(R_top)-log(R_top-Vr*t)) + deg2rad(NorthOffset);	% Ângulo [rad]
R = R_top-Vr*t;                                 % Raio [m]


Px = R.*cos(theta);                            
Py = R.*sin(theta);                             
Pz = H_high-Vd*t;                               
end

