function [out, eSigma] = mmseNoiseReduce(in, prev, eSigmaPrev, prior)
    amp = abs(in);
    %phase = angle(in);
    amp2 = amp*amp;
    % step 1. compute xi(k, i)
    y = amp/eSigmaPrev;
    xi = max(y - 1, 0);
    % step 2. compute E(N^2|y;xi(k, i))
    E = ((1/(1+xi)^2) + (xi/((1+xi)*prior)))*amp2;
    % step 3. get xiDD using decision-directed approach https://ieeexplore.ieee.org/abstract/document/1164453
    eA = (E/(1+E))*abs(prev);
    x = 0.96; % according to slides this should be between 0.96 - 0.99
    xidd = x*(eA/eSigmaPrev) + (1-x)*max(0, y - 1);
    % step 4. B(xidd)
    g = gammainc(2, 1/(1+Edd));
    B = 1/((1+xidd)*g + exp(-(1/(1+xidd))));
    % step 5. sigma estimate
    eSigma = E * B;
    % step 6. smooth it out
    betha = 0.8;
    eSigma = betha*eSigmaPrev + (1 - betha)*eSigma;
    % posterior = amp2/eSigma;
    out = in - eSigma; %?
end

