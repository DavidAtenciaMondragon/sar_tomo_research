function strRadar = converterPosGrados2Rad(strRadar)

switch strRadar.tipoTrajetoria
    case 'lineal'
        
    case 'circular'
        strRadar.posInicial(2) = strRadar.posInicial(2)*pi/180;
        strRadar.posInicial(3) = strRadar.posInicial(3)*pi/180;
    otherwise
        
end

end