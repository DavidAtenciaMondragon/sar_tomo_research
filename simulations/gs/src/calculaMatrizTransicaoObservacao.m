function [F, H] = calculaMatrizTransicaoObservacao(strRadar)

switch strRadar.tipoTrajetoria 
    case 'lineal'

        F = [1 0 0 strRadar.PRT 0 0; ...
             0 1 0 0 strRadar.PRT 0; ...
             0 0 1 0 0 strRadar.PRT; ...
             0 0 0 1 0 0; ...
             0 0 0 0 1 0; ...
             0 0 0 0 0 1];
        
        H = [1 0 0 0 0 0; ...
             0 1 0 0 0 0; ...
             0 0 1 0 0 0];

    case 'circular'
        
        F = [1, 0, 0, 0;  
             0, 1, 0, strRadar.PRT; 
             0, 0, 1, 0;
             0, 0, 0, 1];
         
         H = [1, 0, 0, 0;
              0, 1, 0, 0;
              0, 0, 1, 0];
         
         
    otherwise
        
end

end