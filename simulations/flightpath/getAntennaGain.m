function G = getAntennaGain(tg, Tx_pos, az_Tx_ref, el_Tx_ref, radPattern)

    PxT = Tx_pos(1,:);
    PyT = Tx_pos(2,:);
    PzT = Tx_pos(3,:);

    delta_Tx_tg =  tg - [PxT; PyT; PzT];
    
    [az_Tx_tg, el_Tx_tg, ~] = cart2sph(delta_Tx_tg(1,:), delta_Tx_tg(2,:), delta_Tx_tg(3,:));
    
    % Preproc az
    az_Tx_tg(az_Tx_tg < 0) = az_Tx_tg(az_Tx_tg < 0) + 2*pi;
    
    az_delta_pattern = az_Tx_ref - az_Tx_tg;
    az_delta_pattern(az_delta_pattern > pi) = az_delta_pattern(az_delta_pattern > pi) - 2*pi;
    
    idx_az = (az_delta_pattern + pi/2) * 180/pi;
    
    % Preproc el 
    el_delta_pattern = el_Tx_ref - el_Tx_tg;
    
    idx_el = (el_delta_pattern + pi/2) * 180/pi;
    
    % Limit indices
    idx_az(idx_az < 1) = 1;
    idx_el(idx_el < 1) = 1;
    idx_az(idx_az > size(radPattern,2)) = size(radPattern,2);
    idx_el(idx_el > size(radPattern,1)) = size(radPattern,1);

    % Obtener valores de Gt com sub2ind
    idx = sub2ind(size(radPattern), round(idx_el), round(idx_az));

    % Get value 
    G = radPattern(idx);

end