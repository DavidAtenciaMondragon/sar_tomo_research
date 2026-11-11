function o = osc_clase(nombre)
% OSC_CLASE  Modelo de oscilador del Rx2 (valores tipicos de catalogo).
%   eps_max   : offset de frecuencia fraccional maximo (|eps| <= eps_max)
%   gamma_max : deriva de frecuencia (envejecimiento/termica), 1/s
%   sigma_y1  : desviacion de Allan a 1 s (modelo FM blanca -> el error de
%               tiempo es un paseo aleatorio con incremento sigma_y1*sqrt(dt))
%   jit_rms   : jitter blanco de tiempo por ventana (comun a toda la ventana)
%   phi_lo_rms: ruido de fase blanco extra del LO por ventana (rad)
switch upper(nombre)
    case 'OCXO'
        o = struct('eps_max',5e-8,'gamma_max',1e-12,'sigma_y1',1e-11, ...
                   'jit_rms',1e-12,'phi_lo_rms',0.005);
    case 'TCXO'
        o = struct('eps_max',1e-6,'gamma_max',1e-10,'sigma_y1',1e-10, ...
                   'jit_rms',2e-12,'phi_lo_rms',0.01);
    case 'XO'
        o = struct('eps_max',2e-5,'gamma_max',1e-9,'sigma_y1',1e-9, ...
                   'jit_rms',5e-12,'phi_lo_rms',0.02);
    case 'XO_BARATO'
        o = struct('eps_max',5e-5,'gamma_max',1e-8,'sigma_y1',1e-8, ...
                   'jit_rms',1e-11,'phi_lo_rms',0.05);
    otherwise
        error('osc_clase: clase desconocida %s', nombre);
end
o.nombre = upper(nombre);
end
