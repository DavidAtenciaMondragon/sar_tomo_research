function [chirp] = createCHRIP(strRadar)

    f0 = strRadar.FreqMenor;
    f1 = strRadar.FreqMayor;
    T  = strRadar.DuracaoPulso;
    fs = strRadar.fs;

    t = 0:1/fs:T;
    k = (f1 - f0)/T;
    chirp = exp(1j*2*pi*(f0*t + k/2*t.^2));
end
