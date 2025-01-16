%
clearvars

% Import
cd('/net/storage/gianlum33/projects/oop_ciss_calculations/')
expName = 'data_extracted/zech_p46_oopEseem_expData.csv';
fitName = 'data_extracted/zech_p46_oopEseem_fit.csv';

importedData = readtable(expName);
x = importedData{:, 1}*pi/180;
y = importedData{:, 2};
importedData = readtable(fitName);
xx = importedData{:, 1}*pi/180;
yy = importedData{:, 2};

crystalEseem = @(p, beta) -( ...
    0.5*sin(2*beta) * sin(2*p)^4 + ...
    sin(beta) * cos(2*p)^2 * sin(2*p)^2);  % p = alpha

p0Z = 31.5*pi/180;
initFit = rescaledata(crystalEseem(p0Z, xx), 'maxabs');

crystalEseem = @(p, beta) -(...
    1/2 * sin(2*beta) * cos(p)^4 + ...
    1/4 * sin(beta) * sin(2*p)^2);  % p = 2alpha

p0B = pi/2 - 2*p0Z;  % xi_Bittl = pi/2 - 2*alpha_Zech
initFit2 = rescaledata(crystalEseem(p0B, xx), 'maxabs');

clf
plot(x*180/pi, rescaledata(y, 'maxabs'), 'o-')
hold on
plot(xx*180/pi, rescaledata(yy, 'maxabs'))
plot(xx*180/pi, initFit)
plot(xx*180/pi, initFit2)
% esfit()


%%

% xi = atan((d+2*J)/(OmegaA - OmegaB));
Rz = @(x) [cos(x) sin(x)    0; 
            -sin(x) cos(x)  0;
            0       0       1];
Ry = @(x) [cos(x)   0   -sin(x); 
           0        1   0;
           sin(x)   0   cos(x)];
rotMatrixZyz = @(euAngles) Rz(euAngles(3))*Ry(euAngles(2))*Rz(euAngles(1));
       
eulerAnglesZech = [81, 126, 182]*pi/180;  % Zech, 'A structural model ...'
eulerAnglesZech = [-182, 126, 81]*pi/180;  % Zech, 'A structural model ...'

Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-
eulerAnglesEasyspin = Sys.gFrame(1, :);
% Sys.J = -unitconvert(1e-3,'mT->MHz'); % MHz
% Sys.dip = unitconvert(+0.177,'mT->MHz'); % MHz
% Sys.eeFrame = [0 90 0]*pi/180;

B0 = 0.35;  % T, static magnetic field
JJ = unitconvert(-1e-3,'mT->MHz');
dip = unitconvert(0.177,'mT->MHz');

% Dipolar interaction, p = [theta, phi]
dFunc = @(p) dip*((cos(p(2))*sin(p(1)))^2 - 1/3);

% Inhomogenous line broadening due to hfi (isotropic)
% Prisner et al., Time-resolved W-band..., 1995
lw1 = 15e6;  % MHz
lw2 = 15e6;  % MHz
% The g-factor linewidth
glw1 = planck/bmagn/B0*lw1;
glw2 = planck/bmagn/B0*lw2;
% Spacing between g-values
dgAxis = 0.0001;

thetas = 0:0.5:180;
nTheta = numel(thetas);
phis = 0:0.5:360;
nPhi = numel(phis);
[g1, g2, dd, valueOfCosThetaD] = deal(zeros(nTheta, nPhi));
[nVers, nVersInFrame1, xi] = deal({});

wbarTh = waitbar(0, '1', 'Name', 'xi', ...
    'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
% wbarg1 = waitbar(0, '1', 'Name', 'g1', ...
%     'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');

for ith = 1:nTheta        
    % Update waitbar and message
    waitbar(ith/nTheta, wbarTh, sprintf('%d out of %d', ith, nTheta))
    
    theta_ = thetas(ith)*pi/180;
    % theta_ = theta_/multTheta;
    % theta_ = theta_ * 10;
    for iph = 1:nPhi
        phi_ = phis(iph)*pi/180;
        
        % Direction of B0 in the frame of spin 2
        nVers{ith, iph} = ...
            [cos(phi_)*sin(theta_), sin(phi_)*sin(theta_), cos(theta_)];

        % Angle between zD and B0, the sign is not important
        valueOfCosThetaD(ith, iph) = cos(phi_)*sin(theta_);
        dd(ith, iph) = dFunc([phi_, theta_]);
        g2(ith, iph) = sqrt( sum( (Sys.g(2, :).*nVers{ith, iph}).^2));
        
        % g1InFrame2 = Rz(eulerAngles(3))*Ry(eulerAngles(2))*Rz(eulerAngles(1))*Sys.g(1, :)';
%         g1InFrame2 = ...
%             Rz(eulerAngles2(3))*Ry(eulerAngles2(2))*Rz(eulerAngles2(1))*g1In1*( ...
%             Rz(eulerAngles2(1))^-1 *Ry(eulerAngles2(2))^-1 *Rz(eulerAngles2(3))^-1);
%         g1(ith, iph) = sqrt( sum( (diag(g1InFrame2).*nVers').^2));
        nVersInFrame1{ith, iph} = ...
            rotMatrixZyz(eulerAnglesZech)*nVers{ith, iph}';
        g1(ith, iph) = sqrt( sum( (Sys.g(1, :).*nVersInFrame1{ith, iph}').^2));

        xi(ith, iph) = atan( (dd(ith, iph) + 2*JJ)*1e6 * ...
            planck/(2*pi*(bmagn*0.35)) / (g1(ith, iph) - g2(ith, iph)) );
    end
end

delete(wbarTh)

% TODO maybe set same colorbar extremes for both
figure(1)
clf
tL = tiledlayout(1, 2);
nexttile
imagesc(phis, thetas, g1)
colorbar
title('g1')
nexttile
imagesc(phis, thetas, g2)
colorbar
labelaxesfig(tL, 'Phi', 'Theta')
title('g2')

figure(2)
clf
% imagesc(dip*(valueOfCosThetaD.^2 - 1/3))
% imagesc(thetas, phis, dd)
imagesc(phis, thetas, dd)
labelaxesfig(gca, 'Phi', 'Theta')
% imagesc(thetas, phis, xi)
cbar = colorbar;
% set(cbar, 'YDir', 'reverse');
title('Dipolar interaction dd')


figure(3)
clf
imagesc(phis, thetas, g1 - g2)
labelaxesfig(gca, 'Phi', 'Theta')
cbar = colorbar;
title('g1 minus g2')

figure(4)
clf
% imagesc(dip*(valueOfCosThetaD.^2 - 1/3))
imagesc(phis, thetas, xi)
cbar = colorbar;
labelaxesfig(gca, 'Phi', 'Theta')
% set(cbar, 'YDir', 'reverse');
title('xi')
%
% 
% clf
% plot(thetas, xi(:, 360), '-o')

%

% [nTheta, nPhi] = size(xi);
signalPowder = zeros(numel(xx), 1);

wbarTh = waitbar(0, '1', 'Name', 'signalPowder');
% wbarPh = waitbar(0, '1', 'Name', 'In the for cycle (Phi)');
% 
for ith = 1:nTheta
    % Update waitbar and message
    waitbar(ith/nTheta, wbarTh, sprintf('%d out of %d', ith, nTheta))
    
    for iph = 1:nPhi
%         waitbar(iph/nPhi, wbarPh, sprintf('%d out of %d', iph, nPhi))
        signal_{ith, iph} = crystalEseem(xi(ith, iph), xx);
        % Solid angle normalization
        theta_ = thetas(ith)*pi/180;
        signal_{ith, iph} = signal_{ith, iph}*sin(theta_);
        signalPowder = signalPowder + signal_{ith, iph};
    end
end

delete(wbarTh)
% delete(wbarPh)


figure(5)
clf
plot(x*180/pi, rescaledata(y, 'maxabs'), 'o-')
hold on
plot(xx*180/pi, rescaledata(yy, 'maxabs'))
plot(xx*180/pi, initFit)
plot(xx*180/pi, initFit2)
plot(xx*180/pi, rescaledata(signalPowder, 'maxabs'), 'ko-')
legend('Exp. data', 'Fit Zech', 'Fit Gianluca using alpha', ...
    'Fit Gianluca using xi', 'Sum of ESEEM signals', ...
    'Location', 'northwest')
labelaxesfig(gca, 'Flip angle beta', 'ESEEM intensity')
% plot(xx*180/pi, -rescaledata(sin(2*xx), y))

%% Mean xi

% for ith = 1:nTheta
%     theta_ = thetas(ith)*pi/180;
%     xitest(ith, :) = repmat(cos(theta_)^4, 1, nPhi);
%     xitest(ith, :) = repmat(1, 1, nPhi);
% end

for ith = 1:nTheta
    theta_ = thetas(ith)*pi/180;
    xiMean = xiMean + sum(xi(ith, :))*sin(theta_);
end
normTh = pi/2/nTheta;
normPh = 1/nPhi;
normTot = normPh*normTh;
xiMean = xiMean*normTot*180/pi % /nTheta*pi/2
% xiMean = xiMean/normTh

%%

g1_ = g1(1, 1);
g2_ = g2(1, 1);
lw1 = 15e6;  % Hz
lw2 = 15e6;  % Hz
% The g-factor linewidth for each one is (delta_nu_i = lw_i):
% delta_g_i mu_B B0 = h delta_nu_i
glw1 = planck/bmagn/0.35*lw1;
glw2 = planck/bmagn/0.35*lw2;

xx = creategaxis(g1_, g2_, glw1, glw2, 0.0001);
% xx = 1.995:0.0001:2.015;

gd1 = gaussian(xx, g1_, glw1);
gd2 = gaussian(xx, g2_, glw2);

figure()
plot(xx, gd1, xx, gd2)


%%

function xx = creategaxis(g1, g2, glw1, glw2, dg)
    xxFactor = 1.5;
    xmin = min([g1 - xxFactor*glw1, g2 - xxFactor*glw2]);
    xmax = max([g1 + xxFactor*glw1, g2 + xxFactor*glw2]);
    xx = xmin:dg:xmax;
end















