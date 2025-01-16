function [thetas, phis] = createthetaphigrid(nTheta, nPhi, gridType)
% [thetas, phis] = createthetaphigrid(nTheta, nPhi)
% 
% Output thetas (nTheta, 1) and phis (1, nPhi) in the interval from 0 to pi
% Note: avoids theta = 0 and theta = pi. Avoids phi = pi.

if strcmpi(gridType, "zech")
    costhetas = linspace(1 - 1/nTheta, -1 + 1/nTheta, nTheta);
    thetas = acos(costhetas);
else
    thetas = linspace(0, pi, nTheta + 2);
    thetas = thetas(2:end - 1);
end
thetas = thetas';

phis = linspace(0, 180, nPhi + 1)*pi/180;
phis = phis(1:end - 1);

end

