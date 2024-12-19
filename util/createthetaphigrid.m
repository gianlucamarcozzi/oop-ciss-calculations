function [thetas, phis] = createthetaphigrid(nTheta, nPhi)
% [thetas, phis] = createthetaphigrid(nTheta, nPhi)
% 
% Output thetas (nTheta, 1) and phis (nPhi, 1) in the interval from 0 to pi
% Note: avoids theta = 0 and theta = pi. Avoids phi = pi.
    
thetas = linspace(0, 180, nTheta + 2)*pi/180;
thetas = thetas(2:end - 1);
thetas = thetas';
phis = linspace(0, 180, nPhi + 1)*pi/180;
phis = phis(1:end - 1);
end

