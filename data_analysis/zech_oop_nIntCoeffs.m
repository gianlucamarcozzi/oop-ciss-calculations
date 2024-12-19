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
% Effective g-values
g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
g1 = sqrt( sum( (g1Tensor*rVers).^2, 1));
g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
g2 = sqrt( sum( (g2Tensor*rVers).^2, 1));
% Dipolar interaction
zD = erot(Sys.eeFrame)*[0, 0, 1]';
dd = dipinteraction(dip, rVers, zD);
deltaw = bmagn*B0/hbar*(g1 - g2);
xi = mixinganglexi(deltaw, JJ*1e6, dd*1e6);

imagesc(reshape(xi, [nTheta, nPhi]))
colorbar

cBeta0 = @(tau) sin(thetas)'.*sin(2*(JJ - dd).*tau).*sin(2*xi).^2;
cTwoBeta0 = @(tau) sin(thetas)'.*sin(2*(JJ - dd).*tau).*cos(xi).^4;

nTau = 1000;
taus = linspace(0, 10, nTau);
cBeta = zeros(1, nTau);
cTwoBeta = zeros(1, nTau);
cBetaN = zeros(1, nTau);
cTwoBetaN = zeros(1, nTau);
for ii = 1:nTau
    cBetaN(ii) = averageoversolidangle(cBeta0(taus(ii)), 91, 91, 2);
    cTwoBetaN(ii) = averageoversolidangle(cTwoBeta0(taus(ii)), 91, 91, 2);
    cTwoBetaN2(ii) = sum(sin(thetas).*..cTwoBeta0(taus(ii))/...
        sum(sin(thetas)*ones(1, nPhi), 'all'), 'all');
    cBeta(ii) = sum(cBeta0(taus(ii)), 'all');
    cTwoBeta(ii) = sum(cTwoBeta0(taus(ii)), 'all');
end

figure(1)
clf
plot(taus, cBeta)
hold on
plot(taus, cTwoBeta)
yyaxis right
% plot(taus, cTwoBetaN2)
% hold on
% plot(taus, cTwoBeta/90)



figure(2)
clf
plot(taus, cBetaN + cTwoBetaN)