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

fieldAxis = linspace(343.6, 347.6, 501);  % GHz
mwFreq = 9.7;  % GHz

Opt.GridSize = [91 0];
Opt.GridSymmetry = '';

% Field sweep
Exp.CenterSweep = [mean(fieldAxis), max(fieldAxis) - min(fieldAxis)];
Exp.mwFreq = mwFreq;
Exp.Harmonic = 0;
Exp.nPoints = numel(fieldAxis);

Sys.lwpp = 0.035;  % mT

ySim = pepper(Sys, Exp, Opt);
ySim = ySim/max(abs(ySim));

clf
plot(fieldAxis, ySim)

%% My powder average (only for equivalent nuclei)

mySys.g = Sys.g;
mySys.gFrame = Sys.gFrame;
mySys.dip = mt2mhz(-0.170);  % MHz
mySys.J = mt2mhz(1e-3);  % MHz
mySys.eeFrame = Sys.eeFrame;
mySys.lwpp = Sys.lwpp;
mySys.mwFreq = mwFreq;
mySys.x = fieldAxis;
mySys.nTheta = 91;
mySys.nPhi = 91;

[myySim, saw1, aa] = mytrepr(mySys);

hold on
plot(fieldAxis, myySim)

%%

function ySim = mytrepr(Sys)

    fieldAxis = Sys.x;
    dip = Sys.dip;
    JJ = Sys.J;
    mwFreq = Sys.mwFreq;
    nTheta = Sys.nTheta;
    nPhi = Sys.nPhi;
    
    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi);
   
    nSolidAngle = nTheta*nPhi;

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
    
    w0 = g2wFunc(fieldAxis', (g1 + g2)/2);  % MHz
    deltaw = g2wFunc(fieldAxis', (g1 - g2)/2);  % MHz
    Omega = hypot(deltaw, JJ + dd/2);  % MHz

    ddRep = repmat(dd, [numel(fieldAxis), 1]);
    reshapeSize = [numel(fieldAxis), 1, nSolidAngle];
    wReson = ...
        pagemtimes(reshape(w0, reshapeSize), [1, 1, 1, 1]) + ...
        pagemtimes(JJ - reshape(ddRep, reshapeSize), [-1, 1, -1, 1]) + ...
        pagemtimes(reshape(Omega, reshapeSize), [-1, -1, 1, 1]);

    intensityReson = 1/8*(deltaw.^2)./(Omega.^2);
    intensityReson = ...
        pagemtimes(reshape(intensityReson, reshapeSize), [1, -1, 1, -1]);

    trSignalTrans = gaussianbs(...
        wReson*1e-3, mwFreq, mt2mhz(Sys.lwpp)*1e-3, "lwpp");

    trSignalTrans = trSignalTrans.*intensityReson;
    trSignal = squeeze(sum(trSignalTrans, 2));  % Sum over transitions

    % Solid angle average
    solidAngleWeight = repmat(sin(thetas), [nPhi, 1])/...
        sum(repmat(sin(thetas), [nPhi, 1]), 'all');
    ySim = squeeze(sum(trSignal.*solidAngleWeight', 2));
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

function y1 = gaussianbs(wReson, mwFreq, c, mode)
    arguments
        wReson
        mwFreq  
        c    (1, 1) double     
        mode string = "lwpp"
    end    

    normTerm = 1; %sqrt(2/pi)/c;
    y1 = normTerm*exp(-2 * (wReson - mwFreq).^2 ./ c^2);

end





