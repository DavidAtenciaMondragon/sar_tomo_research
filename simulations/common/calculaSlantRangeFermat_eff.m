function [strReflexao, strRefraccoes] = calculaSlantRangeFermat_eff(terrain, Tx, Rx, p, n1, n2, factor_interp)

c = physconst('LightSpeed');

x = terrain.X_vec;
y = terrain.Y_vec;
z = terrain.Z_DEM;

sizeGrid = length(x);

x_interp = linspace(min(x(:)),max(x(:)),sizeGrid*factor_interp);
y_interp = linspace(min(y(:)),max(y(:)),sizeGrid*factor_interp);

[x_grid_interp,y_grid_interp] = meshgrid(x_interp,y_interp);

z_interp = interp2(x,y,z,x_grid_interp,y_grid_interp);

pos_interface = [x_grid_interp(:), y_grid_interp(:), z_interp(:)];

% dist 
Tx_above = pdist2(pos_interface,Tx);
Tx_below = pdist2(pos_interface,p);
Rx_above = pdist2(pos_interface,Rx);
Rx_below = pdist2(pos_interface,p);

% Pre alocate
coder.varsize('tempo_Tx_P', [20000, 20000], [true, true])
coder.varsize('tempo_Rx_P', [20000, 20000], [true, true])
coder.varsize('tempo_Tx_Rx', [20000, 20000], [true, true])

% Calcular tiempo de viaje
tempo_Tx_P  = (Tx_above) / (c/n1) + (Tx_below) / (c/n2);
tempo_Rx_P  = (Rx_above) / (c/n1) + (Rx_below) / (c/n2);
tempo_Tx_Rx = (Tx_above + Rx_above) / (c/n1);

% --- Calcular el mínimo para cada trayectoria (para cada columna) ---
[~, idx_min_Tx_P]  = min(tempo_Tx_P, [], 1);
[~, idx_min_Rx_P]  = min(tempo_Rx_P, [], 1);
[~, idx_min_Tx_Rx] = min(tempo_Tx_Rx, [], 1);

x_flat = x_grid_interp(:);
y_flat = y_grid_interp(:);
z_flat = z_interp(:);

Preflex  = [x_flat(idx_min_Tx_Rx), y_flat(idx_min_Tx_Rx), z_flat(idx_min_Tx_Rx)];
Qrefract = [x_flat(idx_min_Tx_P),  y_flat(idx_min_Tx_P),  z_flat(idx_min_Tx_P)];
R_PRx    = [x_flat(idx_min_Rx_P),  y_flat(idx_min_Rx_P),  z_flat(idx_min_Rx_P)];

strReflexao.P_reflex = Preflex.';
strReflexao.Ang_incidencia_deg = 0;
strReflexao.Ang_reflexion_deg = 0;

strRefraccoes.P_refrac_ida = Qrefract.';
strRefraccoes.P_refrac_volta = R_PRx.';
strRefraccoes.Ang_incidencia_ida_deg = 0;
strRefraccoes.Ang_refrac_ida_deg = 0;
strRefraccoes.Ang_incidencia_volta_deg = 0;
strRefraccoes.Ang_refrac_volta_deg = 0;

end