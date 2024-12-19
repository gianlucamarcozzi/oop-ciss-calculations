clearvars
addpath(genpath('../util'))

%% Definitions

Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

B0 = 0.35;
JJ = mt2mhz(0.001);
dip = -mt2mhz(0.170);

nTheta = 91;
nPhi = 91;
[thetas, phis] = createthetaphigrid(nTheta, nPhi);
nSolidAngle = nTheta*nPhi;

rVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
size(rVers)
% Effective g-values
g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
g1 = sqrt( sum( (g1Tensor*rVers).^2, 1));
g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
g2 = sqrt( sum( (g2Tensor*rVers).^2, 1));
% Dipolar interaction
zD = erot(Sys.eeFrame)*[0, 0, 1]';
dd = dipinteraction(dip, rVers, zD);

for ii = 1:nSolidAngle
    dipolar = dd(ii);
    % angTheta = thetas(ii);
    % angPhi = phis(ii);
    gvalue1 = g1(ii);
    gvalue2 = g2(ii);
    xi = atan((dipolar + 2*JJ)/...
        (bmagn*B0/hbar*(gvalue1 - gvalue2)));

    integral2(@(x,y) x.^2 + y.^2, 0, 1, 0, 1);
end
