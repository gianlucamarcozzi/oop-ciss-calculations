%
clearvars

%% Easyspin simulation

% Taken from Easyspin example

Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.J = (-1/2)*unitconvert(1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(0.170,'mT->MHz')/1.5; % MHz
Sys.initState = 'singlet';

freqAxis = linspace(9.63, 9.73, 301);  % GHz
B0 = 345;  % mT

Opt.GridSize = [91 0];  
Opt.GridSymmetry = '';

% Frequency sweep
Exp.mwCenterSweep = [mean(freqAxis), max(freqAxis) - min(freqAxis)];
Exp.Field = B0;
Exp.Harmonic = 0;
Exp.nPoints = numel(freqAxis);

Sys.lwpp = 10;  % MHz

ySim = pepper(Sys, Exp, Opt);
ySim = ySim/max(abs(ySim));

clf
plot(freqAxis, ySim)

%% My powder average (only for equivalent nuclei)

mySys.g = Sys.g;
mySys.gFrame = Sys.gFrame;
mySys.dip = mt2mhz(-0.170);  % MHz
mySys.J = mt2mhz(1e-3);  % MHz
mySys.eeFrame = Sys.eeFrame;
mySys.lwpp = Sys.lwpp;
mySys.Field = B0;
mySys.x = freqAxis;
mySys.nTheta = 91;
mySys.nPhi = 91;

myySim = mytrepr(mySys);

hold on
plot(freqAxis, myySim)

%%

function ySim = mytrepr(Sys)

    freqAxis = Sys.x;
    dip = Sys.dip;
    JJ = Sys.J;
    B0 = Sys.Field;
    nTheta = Sys.nTheta;
    nPhi = Sys.nPhi;
    
    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi);

    % Direction of B0
    nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
    % Effective g-values
    g1Tensor = rotatematrix(diag(Sys.g(1, :)), Sys.gFrame(1, :));
    g1 = sqrt( sum( (g1Tensor*nVers).^2, 1));
    g2Tensor = rotatematrix(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
    g2 = sqrt( sum( (g2Tensor*nVers).^2, 1));
    % Dipolar interaction
    zD = erot(Sys.eeFrame)*[0, 0, 1]';
    dd = dipFunc(dip, nVers, zD);
    
    w0 = g2wFunc(B0, (g1 + g2)/2);  % MHz
    deltaw = g2wFunc(B0, (g1 - g2)/2);  % MHz
    Omega = hypot(deltaw, JJ + dd/2);  % MHz

    wReson = ...
        [1; 1; 1; 1]*w0 + [-1; 1; -1; 1]*(JJ - dd) + [-1; -1; 1; 1]*Omega;

    intensityReson = [1; -1; 1; -1]*1/8*(deltaw.^2)./(Omega.^2);
    intensityReson = reshape(intensityReson, [1, size(intensityReson)]);

    trSignalTrans = gaussiantransitions(...
        freqAxis', wReson*1e-3, Sys.lwpp*1e-3, "lwpp");

    trSignalTrans = trSignalTrans.*intensityReson;
    trSignal = squeeze(sum(trSignalTrans, 2));  % Sum over transitions

    solidAngleWeight = repmat(sin(thetas), [nPhi, 1])/...
        sum(repmat(sin(thetas), [nPhi, 1]), 'all');
    ySim = sum(squeeze(sum(trSignal.*solidAngleWeight', 3)), 2);
    ySim = ySim/max(abs(ySim));

end

function freq = g2wFunc(B0, g)
    % Expected input: B0 in mT
    % Output: frequency nu in MHz
    freq = bmagn*B0/planck*g*1e-9;
end

function g = w2gFunc(B0, freq)
    % Expected input: B0 in mT
    % Output: frequency nu in MHz
    g = planck/bmagn/B0*freq*1e9;
end

function dd = dipFunc(dip, nVers, zD)
    % Dipolar interaction
    % Makes use of the fact that cos(thetaD) = dotProduct(B0, zD)
    % Expected input:
    %   dip:    1 x 1
    %   nVers:  3 x nTheta x nPhi
    %   zD:     3 x 1
    dd = dip*((sum(nVers.*zD)).^2 - 1/3);
end

function [thetas, phis] = createthetaphigrid(nTheta, nPhi)
    
    thetas = linspace(0, 180, nTheta + 2)*pi/180;
    thetas = thetas(2:end - 1);
    thetas = thetas';
    phis = linspace(0, 180, nPhi + 1)*pi/180;
    phis = phis(1:end - 1);
    
end

function newMat = rotatematrix(mat, euAngles)
    % newMat = R*mat*transpose(R)
    % euAngles are the ones given in the Easyspin convention
    
    euMatrix = erot(euAngles)';
    newMat = euMatrix*mat*transpose(euMatrix);
    
end

function y1 = gaussiantransitions(xx, x0, c, mode)
    % Output:
    % y1:   nAx x 4 x nSolidAngle
    arguments
        xx  (:, 1) double       % nA x 1 equally spaced values
        x0  (4, :) double       % 4 x nSolidAngle
        c   (1, 1) double
        mode string = "lwpp"
    end
    
    if strcmp("mode", "lw") || strcmp("mode", "fwhm")
        c = c/(sqrt(2*log(2)));
    end
    nAx = numel(xx);
    xx = repmat(xx, [1, size(x0)]);
    x0 = repmat( ...
            reshape(x0, [1, size(x0)]), ...
            [nAx, 1, 1]);

    normTerm = sqrt(2/pi)/c;
    y1 = normTerm*exp(-2 * (xx - x0).^2 ./ c^2);

end





