function powerOutGain = EqRadar(strRadarTx, strRadarRx, strTarget, R1 , R2)

Gt   = strRadarTx.GananciaAnt;
Pt   = strRadarTx.PotenciaTx;
lamb = strRadarTx.lamb;
RCS  = strTarget.rcs;

Gr   = strRadarRx.GananciaAnt;

F    = 1;

powerOutGain = Pt * Gt *Gr * RCS * (lamb ^ 2) * (F ^ 4) / ...
                (((4 * pi) ^ 3) * (R1 ^ 2) * (R2 ^ 2));

end
