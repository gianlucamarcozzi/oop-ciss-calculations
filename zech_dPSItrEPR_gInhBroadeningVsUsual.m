%
clearvars

%% Easyspin simulation

% Taken from Easyspin example
Sys.S = [1/2 1/2];
Sys.g = [2.00308 2.00264 2.00226; ... % P700+
         2.00622 2.00507 2.00218];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.J = (-1/2)*unitconvert(1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(0.170,'mT->MHz')/1.5; % MHz
Sys.initState = 'singlet';

fieldAxis = linspace(343.635, 347.635, 301);  % GHz
mwFreq = 9.7;  % GHz

Opt.GridSize = [91 0];
Opt.GridSymmetry = '';

% Field sweep
Exp.CenterSweep = [mean(fieldAxis), max(fieldAxis) - min(fieldAxis)];
Exp.mwFreq = mwFreq;
Exp.Harmonic = 0;
Exp.nPoints = numel(fieldAxis);

Sys.lwpp = 0.35;  % mT

ySim = pepper(Sys, Exp, Opt);
ySim = ySim/max(abs(ySim));

clf
plot(fieldAxis, ySim)

%%

% pathFolder = "/net/storage/gianlum33/";
pathFolder = "S:\";
pathFolder = pathFolder + "projects/oop_ciss_calculations/data/digitized/";
pathX = pathFolder + "expData_zech_aStructuralModelFor_xBandDeut.mat";
pathQ = pathFolder + "expData_zech_aStructuralModelFor_qBandDeut.mat";
pathW = pathFolder + "expData_zech_aStructuralModelFor_wBandDeut.mat";

importedData = load(pathX);
xdataXdeut = importedData.x;
ydataXdeut = importedData.y;
ydataXdeut = ydataXdeut - (ydataXdeut(1) + ydataXdeut(end))/2;
importedData = load(pathQ);
xdataQdeut = importedData.x;
ydataQdeut = importedData.y;
ydataQdeut = ydataQdeut - (ydataQdeut(1) + ydataQdeut(end))/2;
importedData = load(pathW);
xdataWdeut = importedData.x;
ydataWdeut = importedData.y;

importedData = load('S:\projects\zech_psi\data\processed\ZePSI-E-007015.mat');
% clf
% h = ScrollableAxes();
% plot(h, x{2}, x{1}, y');
xdataX = importedData.x{2}'/10;  % mT
ydataX = importedData.y(815, :);
ydataX = ydataX/max(abs(ydataX));
pathQ2 = pathFolder + "expData_zech_aStructuralModelFor_qBand.mat";
importedData = load(pathQ2);
xdataQ = importedData.x;
ydataQ = importedData.y;

%% My powder average (only for equivalent nuclei)

mySys.g = Sys.g;
mySys.gFrame = Sys.gFrame;
mySys.dip = mt2mhz(-0.170);  % MHz
mySys.J = mt2mhz(1e-3);  % MHz
mySys.eeFrame = Sys.eeFrame;
mySys.nNuc = 0;
mySys.A = [0, 0, 0, 9, 9.4, 12.8];  % MHz
mySys.AFrame = [0 0 0 60 -90 0]*pi/180;
mySys.lw1 = 0.3;  % mT
% mySys.lw1 = 0.7499;  % mT
% mySys.lw2 = 0.25;  % mT
mySys.lw2 = 0.3;  % mT

% FitOpt.maxTime = 0.1;
FitOpt.TolFun = 1e-4;
FitOpt.TolEdgeLength = 1e-4;
FitOpt.CalculateUncertainties = 0;
FitOpt.Method = 'levmar';
FitOpt.Gradient = 1e-3;
FitOpt.AutoScale = 'lsq';

Vary.lw1 = mySys.lw1 - 0.005;
Vary.lw2 = mySys.lw2 - 0.005;
Vary.xshift = 0.1;

ydata = {ydataXdeut, ydataQdeut, ydataWdeut};
mySys.mwFreq_MF = {9.7, 34, 94.11};
 % xdataX + 5.5, xdataQ, xdataW + 3.5
mySys.xshift = -0.0245;
mySys.x_MF = {xdataXdeut - 1, xdataQdeut, xdataWdeut + 0.};
mySys.nTheta_MF = {40, 20, 50}; % {40, 40, 20}; %, 50};
mySys.nPhi_MF = {20, 10, 30}; % {20, 20, 20}; %, 20};
mySys.nSamplegInh1_MF = {40, 20, 20};
mySys.nSamplegInh2_MF = {40, 20, 20};
mySys.coeffResLw_MF = {0.05, 0.1, 0.1};

nSim = 15;
xshifts = linspace(-0.04, 0.04, nSim);
wb = waitbar(0, 'Simulting');
ydata = {ydataXdeut};
for ii = 1
    waitbar(ii/nSim, wb, ...
        "Simulation " + string(ii) + " out of " + string(nSim) + "...")
    mySys.xshift = -0.0245;
    tic
    ySim = mytreprmf(mySys);
    Fit = esfit(ydata, @mytreprmf, {mySys}, {Vary}, FitOpt);
    toc

    figure()
    clf
    tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')
    for ifig = 1
        x_ = mySys.x_MF{ifig};
        nexttile
        plot(x_, ydata{ifig})
        hold on
        plot(x_ + mySys.xshift, Fit.fit)
        plot(x_ + mySys.xshift, rescaledata(ySim{ifig}, ydata{ifig}, 'lsq'))
%      xlim(setaxlim(mySys.x_MF{ii}))
%     xlim(setaxlim(1:numel(ydata{ii})))
%     ylim(setaxlim(ydata{ii}, 1.05))
        title(string(xshifts(ii)) + " : " + ...
            string(Fit.pfit(1)) + "   " + string(Fit.pfit(2)))

    end
end
close(wb); clear('wbar')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(1))*1e3, 'r')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(3))*1e3, 'r')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(5))*1e3, 'r')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(2))*1e3, 'b')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(4))*1e3, 'b')
% xline(mhz2mt(mySys.mwFreq_MF{3}, mySys.g(6))*1e3, 'b')
% tic
% ydata = {yX, yQ, yW};
% Fit = esfit(ydata, @mytreprmf, {mySys}, {Vary}, FitOpt);

%%
ySim0 = mytreprmf(mySys);
%%

function data = mytreprmf(Sys)
    
    nSpc = numel(Sys.mwFreq_MF);
    fields = fieldnames(Sys);
    for ispc = 1
        for ifield = 1:numel(fields)
            f = fields{ifield};
            if endsWith(f, '_MF')
                f_ = f(1:end - 3);
                Sys.(f_) = Sys.(f){ispc};
            end
        end
        data{ispc} = mytrepr(Sys);
    end
end

% function [ySim, gax, gInh1, gInh2, gax1, gInh11, gax2, gInh22] = mytrepr(Sys)
function ySim = mytrepr(Sys)
    fieldAxis = Sys.x + Sys.xshift;
    Bmean = mean(fieldAxis);
    dip = Sys.dip;
    JJ = Sys.J;
    mwFreq = Sys.mwFreq;
    nTheta = Sys.nTheta;
    nPhi = Sys.nPhi;
    if isfield(Sys, 'nNuc')
        nNuc = Sys.nNuc;
    else
        nNuc = 0;
    end
    nHfiLine = nNuc + 1;
    
    % Parameters for g broadening
    lw1 = mt2mhz(Sys.lw1);
    lw2 = mt2mhz(Sys.lw2);
    % The g-factor linewidth
    glw1 = freq2gvalue(Bmean, lw1);
    glw2 = freq2gvalue(Bmean, lw2);
    % Spacing between g-values
    dgax1 = glw1 / Sys.nSamplegInh1;
    dgax2 = glw2 / Sys.nSamplegInh2;
    resLw = min(lw1, lw2)*Sys.coeffResLw*1e-3;  % GHz
    
    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi);
    nSolidAngle = nTheta*nPhi;

    % Direction of B0
    nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
    % Effective g-values
    g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
    g1 = sqrt( sum( (g1Tensor*nVers).^2, 1));
    g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
    g2 = sqrt( sum( (g2Tensor*nVers).^2, 1));
    % Dipolar interaction
    zD = erot(Sys.eeFrame)*[0, 0, 1]';
    dd = dipinteraction(dip, nVers, zD);
    
    if nNuc > 0
        [APlus, AMinus] = calculateaplusaminusnohfi(...
            Sys.A, Sys.AFrame, nVers);
    else
        APlus = zeros(nSolidAngle, 1);
        AMinus = zeros(nSolidAngle, 1);
    end
    
    ySim = zeros(numel(fieldAxis), 1);
    parfor ii = 1:nSolidAngle
        nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
        
        [gax1, gInh1] = creategaxisinhbroadening(g1(ii), glw1, dgax1);
        [gax2, gInh2] = creategaxisinhbroadening(g2(ii), glw2, dgax2);
        ngax1 = numel(gax1);
        ngax2 = numel(gax2);
        
        gaussianFactor = gInh1'*gInh2;
        gaussianFactor = gaussianFactor/sum(gaussianFactor, 'all');
        gaussianFactor = gaussianFactor(:)';

        gPlus = 1/2*(repmat(gax1', [1, ngax2]) + repmat(gax2, [ngax1, 1]));
        gPlus = gPlus(:)';
        gMinus = 1/2*(repmat(gax1', [1, ngax2]) - repmat(gax2, [ngax1, 1]));
        gMinus = gMinus(:)';
        %}
        
        w0_ = gvalue2freq(fieldAxis', gPlus);
        deltaw_ = gvalue2freq(fieldAxis', gMinus);
        sintheta_ = sqrt( nVers(1, ii)^2 + nVers(2, ii)^2 );

        quantNumNuc = -nNuc/2:1:nNuc/2;
        signJmd = [-1, 1, -1, 1];
        signOmega = [-1, -1, 1, 1];
        signIntens = [1, -1, 1, -1];
        if nNuc > 0
            pascalMatrix = pascal(nHfiLine);
            % Antidiag
            pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
        else
            pascalFactor = 1;
        end

        for itrans = 1:4
            for ihfi = 1:nHfiLine
                w0__ = w0_ + quantNumNuc(ihfi)*APlus(ii);
                deltaw__ = deltaw_ + quantNumNuc(ihfi)*AMinus(ii);
                Omega__ = hypot(deltaw__, JJ+dd(ii)/2);
                wReson = w0__ + signJmd(itrans)*(JJ - dd(ii)) + ...
                    signOmega(itrans)*Omega__;

                intensityReson = signIntens(itrans)*1/8*...
                    (deltaw__.^2)./(Omega__.^2);
                trSignal__ = gaussianresonancebsweep(...
                    wReson*1e-3, mwFreq, resLw, "lwpp");
                trSignal__ = intensityReson.*trSignal__*...
                    pascalFactor(ihfi).*gaussianFactor;
                ySim = ySim + sum(trSignal__, 2)*sintheta_;
            end
        end
        %}
    end
    %}
%     ySim = ySim/max(abs(ySim));

end
