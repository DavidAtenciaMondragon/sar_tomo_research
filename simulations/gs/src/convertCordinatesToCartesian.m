function cartesianPosition = convertCordinatesToCartesian(position,tipoTrajetoria)

cartesianPosition = zeros(3,1); % Coordenadas cartesianas x,y,z

switch tipoTrajetoria
    case 'lineal'
        cartesianPosition(1) = position(1);
        cartesianPosition(2) = position(2);
        cartesianPosition(3) = position(3);
    case 'circular'
        [cartesianPosition(1),cartesianPosition(2),cartesianPosition(3)] = sph2cart(position(2),position(3),position(1));
    otherwise
        
end
        


end