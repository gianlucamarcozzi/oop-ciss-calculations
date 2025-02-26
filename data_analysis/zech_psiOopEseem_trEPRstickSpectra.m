%
clearvars
addpath(genpath('S:\soft\matlab\'))
addpath(genpath('/net/storage/gianlum33/soft/matlab/'))

%% Easyspin simulation

tic
% Taken from Easyspin example
clear('Sys', 'Exp', 'Vary', 'VaryExp', "Opt")

Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.J = (-1/2)*unitconvert(1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(0.170,'mT->MHz')/1.5; % MHz
Sys.lwpp = 0.35;  % mT
Sys.initState = 'singlet';

nNuc = 0;
Sys.Nucs = '1H';
Sys.Nucs = repmat(append(Sys.Nucs, ', '), [1, nNuc]);
Sys.Nucs = Sys.Nucs(1:end - 2);
Sys.A = [zeros(1, 3), 9, 9.4, 12.8];
% Sys.A = 10*ones(1, 6);
% Sys.A = 0*ones(1, 6);
Sys.A = repmat(Sys.A, [nNuc, 1]);
Sys.AFrame = [0 0 0 60 -90 0]*pi/180;
Sys.AFrame = repmat(Sys.AFrame, [nNuc, 1]);

% gPlus = (Sys.g(1) + Sys.g(2))/2;
mwFreq = 9.6;
% B0 = planck/bmagn/gfree*mwFreq*1e12;  % mT
% xxSim = linspace(mwFreq - 100e-3, mwFreq + 100e-3, 1024);  % GHz
xxSim1 = linspace(339.5, 344.5, 301);  % mT

% Frequency sweep
% Exp.mwCenterSweep = [mean(xxSim), max(xxSim) - min(xxSim)];
% Exp.Field = B0;

% Magnetic field sweep
Exp.CenterSweep = [mean(xxSim1), max(xxSim1) - min(xxSim1)]; % mT
Exp.mwFreq = mwFreq; % GHz

Exp.nPoints = numel(xxSim1);
Exp.Harmonic = 0;

Opt.GridSize = [10 0];  
Opt.GridSymmetry = '';

% Opt.separate = 'transitions';
FitOpt.Verbosity = 0;
FitOpt.maxTime = 0.5;

% Exp.CrystalSymmetry = 1;
% Exp.MolFrame = [0 -0 -0]*pi/180;      % molecular frame aligned with crystal/sample frame
% Exp.SampleFrame = [-0 -0 -0]*pi/180; 
% % 
% nRot_L = [0;0;1];          % sample rotation axis in lab frame (lab x axis here)
% rho = deg2rad(0:1:180);              % 10 degree increments
% Exp.SampleRotation = {'y',rho};
% Opt.separate = 'orientations';
% 

% 
% % ySimNoHfi = pepper(Sys, Exp);

% ySim = pepper(Sys, Exp, Opt);
% yySim = ySim;
% ySim = ySim'/max(abs(ySim));

% Ham = ham(Sys, [0, 0, B0]);
% [Ham0, mux, muy, muz] = ham(Sys);
toc
clf
plot(xxSim1, ySim) %, xxSim, yySim)
% clf

%% Functions
% Expected input: B0 in mT
% Output: frequency nu in MHz
g2wFunc = @(B0, g) bmagn*B0/planck*g*1e-9;

% Expected input: deltaw in MHz, J in MHz, d in MHz
% Output: Omega in MHz
OmegaFunc = @(deltaw, J, d) sqrt(deltaw.^2 + (J + d/2).^2);
% OmegaFuncEasy = @(deltaw, J, d) sqrt(deltaw.^2 + (J/2 + d/2).^2);

% Dipolar interaction
% Makes use of the fact that cos(thetaD) = dotProduct(B0, zD)
% Expected input:
%   dip:    1 x 1
%   nVers:  3 x nTheta x nPhi
%   zD:     3 x 1
dipFunc = @(dip, nVers, zD) dip*((sum(nVers.*zD)).^2 - 1/3);
% dipFuncEasy = @(dip, nVers, zD) -3/2*dip*((sum(nVers.*zD)).^2 - 1/3);

% Mixing angle, p = [dd, JJ, deltaw] in [MHz, MHz, MHz]
xiFunc = @(dd, JJ, deltaw) atan( (dd/2 + JJ) ./ deltaw);

%% My powder average (only for equivalent nuclei)

tic
dip = -mt2mhz(0.170);
JJ = mt2mhz(0.001);
nHfiLine = nNuc + 1;
trlwpp = 0.35;  % mT

% Theta, phi grid
nTheta = 20;
nPhi = 20;
thetas = linspace(0, 180, nTheta)*pi/180;
thetas = thetas';
phis = linspace(0, 180, nPhi)*pi/180;

% Direction of B0
nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
% Effective g-values
euMatrixg1 = erot(Sys.gFrame(1, :))';
g1Tensor = euMatrixg1*diag(Sys.g(1, :))*euMatrixg1';
g1 = sqrt( sum( (g1Tensor*nVers).^2, 1));
euMatrixg2 = erot(Sys.gFrame(2, :))';
g2Tensor = euMatrixg2*diag(Sys.g(2, :))*euMatrixg2';
g2 = sqrt( sum( (g2Tensor*nVers).^2, 1));
% Dipolar interaction
zD = erot(Sys.eeFrame)*[0, 0, 1]';
dd = dipFunc(dip, nVers, zD);

% Hyperfine interaction with a number nNuc of equal nuclei
if nNuc ~= 0
    euMatrixHfi1 = erot(Sys.AFrame(1, 1:3));
    Ahfi1Tensor = euMatrixHfi1'*diag(Sys.A(1, 1:3))*euMatrixHfi1;
    Ahfi1 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi1Tensor', nVers)).^2)));
    % For spin 2
    euMatrixHfi2 = erot(Sys.AFrame(1, 4:6));
    Ahfi2Tensor = euMatrixHfi2'*diag(Sys.A(1, 4:6))*euMatrixHfi2;
    Ahfi2 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi2Tensor', nVers)).^2)));

    % Shift resonances due to hyperfine interaction
    Aplus = (Ahfi1 + Ahfi2)/2;  % MHz
    Aminus = (Ahfi1 - Ahfi2)/2;  % MHz

    % Pascal triangle for hfi interaction
    pascalMatrix = pascal(nHfiLine);
    pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);  % Antidiag
else
    Aplus = 0;
    Aminus = 0;
    pascalFactor = 1;
end

w0 = g2wFunc(xxSim1', (g1 + g2)/2);  % MHz
deltaw = g2wFunc(xxSim1', (g1 - g2)/2);  % MHz
% Signs for different transitions
signJmd = [-1, 1, -1, 1];
signOmega = [-1, -1, 1, 1];
signIntens = [1, -1, 1, -1];
% Possible values of sum of nuclear spin quantum numbers
quantNumNuc = -nNuc/2:1:nNuc/2;
% Initialization
trSignal = zeros(numel(xxSim1), nTheta*nPhi);
xi_ = zeros(4 + nHfiLine, nTheta*nPhi);
% Spectra
for itrans = 1:4
    for ihfi = 1:nHfiLine
        w0_ = w0 + quantNumNuc(ihfi)*Aplus;
        deltaw_ = deltaw + quantNumNuc(ihfi)*Aminus;
        Omega = OmegaFunc(deltaw_, JJ, dd);
        wReson = w0_ + signJmd(itrans)*(JJ - dd) + signOmega(itrans)*Omega;
        % intensityReson correct only for case of equal lws
        intensityReson = signIntens(itrans)*1/8*(deltaw_.^2)./(Omega.^2);
        trSignal_ = ...
            gaussianbs(wReson*1e-3, mwFreq, mt2mhz(trlwpp)*1e-3, "lwpp");
        trSignal = trSignal + intensityReson.*trSignal_*pascalFactor(ihfi);
        % xiWeight = abs(trSignal_)/sum(abs(trSignal_));
        % xi_(itrans + ihfi, :) = sum(xiFunc(dd, JJ, deltaw_).*xiWeight);
    end 
end
% xi = mean(xi_, 1);
trSignal = reshape(trSignal, [numel(xxSim1), nTheta, nPhi]);
% Sum over solid angle distribution
solidAngleWeight = sin(thetas')/sum(sin(thetas)*ones(1, nPhi), 'all');
trSignalPowder = sum(squeeze(sum(trSignal.*solidAngleWeight, 3)), 2);
% Normalized signal to max = 1
ntrSignal = trSignalPowder/abs(max(trSignalPowder));

zfile = load('/net/storage/gianlum33/projects/oop_ciss_calculations/from_bittl/src/test2_spec.asc');
% zfile = load('S:\projects\oop_ciss_calculations\from_bittl\src\test2_spec.asc');
xxz = zfile(:, 1);
yyz = zfile(:, 2);
% figure()
clf
% plot(xxz/10, yyz/max(yyz), 'DisplayName', 'Zech')
hold on
plot(xxSim1, ntrSignal, 'DisplayName', 'Gianluca')
yyaxis right
plot(xxSim1, xixi(1, :))
plot(xxSim1, xixi(2, :))
plot(xxSim1, xixi(3, :))
plot(xxSim1, xixi(4, :))
hold on
% plot(xxSim1, ySim/max(ySim), 'DisplayName', 'Easyspin')
legend()
xlim([min(xxSim1) max(xxSim1)])
toc

%%
%{
% if numel(trlw) == 1
%     trlw = repmat(repelem(trlw, 4), [1, nHfiLine*nTheta*nPhi]);
% elseif numel(trlw) == 2
%     trlw = repmat(repelem(trlw, 2), [1, nHfiLine*nTheta*nPhi]);
% else
%     error('trlw should be of size (1, 1) or (1, 2)')
% end



% Intensities
intensityReson = 1/8*(deltaw(1, :).^2)./(Omega(1, :).^2);
intensityReson = ...
    repmat([1, -1, 1, -1], [1, max(size(dd))]).* ...
    repelem(intensityReson, 4);  
intensityNorm = intensityReson ./ sum(abs(trSignal), 1);
pascalFactor = repmat(repelem(pascalFactor, 4), [1, nTheta*nPhi]);

trSignal = trSignal.*intensityNorm.*pascalFactor;
% trSignalSumHfi = squeeze(sum(trSignal, 3));  % Sum over hfi lines
% trSignalSum = squeeze(sum(trSignalSumHfi, 2));  % Sum over transitions

trSignal = reshape(trSignal, [numel(xxSim1), 4, nHfiLine, nTheta, nPhi]);
if nNuc ~= 0
    trSignalSumHfi = squeeze(sum(trSignal, 2));  % Sum over hfi lines
    trSignalSum = squeeze(sum(trSignalSumHfi, 2));  % Sum over transitions
else
    trSignalSum = squeeze(sum(trSignal, 2));  % Sum over transitions
end

solidAngleWeight = sin(thetas')/sum(sin(thetas)*ones(1, nPhi), 'all');
trSignalPowder = sum(squeeze(sum(trSignalSum.*solidAngleWeight, 3)), 2);

% [~, ith] = min(abs(Exp.SampleFrame(2) + thetas));
% ntrSignal = trSignalSum(:, ith, 1)/abs(max(trSignalSum(:, ith, 1)));
ntrSignal = trSignalPowder/abs(max(trSignalPowder));

diff = ySim - ntrSignal;
% toc
% 
% clf
% hold on
% plot(xxSim, ySim, xxSim, ntrSignal)
% yyaxis right
% plot(xxSim, diff)

% if abs(round(sum(diff), 5)) == abs(round(sum(ySim), 5))
%     disp("Same areas. Area of my norm. signal is: " + string(sum(ntrSignal)))
% else
%     disp("Not the same areas!!!!!!!")
% end

toc



%%
clf
plot(Js, sum(ySim, 2), 'o')
hold on
% plot(Js, sum(diff, 2), 'o')
% plot(Js, sum(ntrSignalPowder, 2), 'o')

%%

% figure()
ij = 1;
clf
plot(xxSim, ySim(ij, :), xxSim, trSignalPowder/max(trSignalPowder))
yyaxis right
plot(xxSim, ySim(ij, :)' - trSignalPowder/max(trSignalPowder))

%%
clf
for ij = 1

    plot(xxSim, diff(ij, :))
    hold on
end
legend()
%%

[~, ith] = min(abs(Exp.SampleFrame(2) - thetas));
% figure(1)
clf
plot(xxSim, ySim/max(ySim))
hold on
plot(xxSim, trSignalSum(:, ith, 1)/max(trSignalSum(:, ith, 1)))
% plot(xxSim, trSignal(:, 1, 1, ith, 1)/max(trSignal(:, 1, 1, ith, 1)))
% plot(xxSim, trSignal(:, 2, 1, ith, 1)/max(trSignal(:, 1, 1, ith, 1)))
% plot(xxSim, trSignal(:, 3, 1, ith, 1)/max(trSignal(:, 1, 1, ith, 1)))
% plot(xxSim, trSignal(:, 4, 1, ith, 1)/max(trSignal(:, 1, 1, ith, 1)))
yyaxis right
plot(xxSim, ySim/max(ySim) - trSignalSum(:, ith, 1)/max(trSignalSum(:, ith, 1)))

%%
[Vary, VaryExp] = deal(struct());
% Vary.gFrame = [1.5, 1.5, 1.5; 0, 0, 0];
% Vary.lw = 5;
% Vary.dip = 1;
% Vary.J = 0.1;
% Vary.A = [zeros(1, 3), 2, 2, 2];
% Vary = {};
% VaryExp.MolFrame = [0 5e-1, 0];
VaryExp.mwCenterSweep = [5e-3, 0];

% ydata = trSignalSum(:, 1, 1);
ydata = trSignalPowder;
ydata = ydata/max(ydata);
% Fit = esfit(ydata, @pepper, {Sys, Exp}, {Vary, VaryExp}, FitOpt);
esfit(ydata, @pepper, {Sys, Exp}, {Vary, VaryExp}, FitOpt);
% clf
% plot(xxSim, ySim/max(ySim), xxSim, ydata) %, xxSim, Fit.fit)

%%

HamTest = [-Sys.dip/2 + Sys.J/4, 0, 0, 0; ...
            0, Sys.dip/2 - Sys.J/4, +Sys.dip/2 + Sys.J/2, 0;
            0, +Sys.dip/2 + Sys.J/2, Sys.dip/2 - Sys.J/4, 0;
            0, 0, 0,-Sys.dip/2 + Sys.J/4 ];
% - g2wFunc(B0, (g1(3) + g2(3))/2) + HamTest(1)
+ g2wFunc(B0, (Sys.g(1, 3) - Sys.g(2, 3))/2) + HamTest(6)

%%

rotation = 1/sqrt(2)*[1,  1;  -1,  1];
hamm = muz(2:3, 2:3);
newHam = rotation'*hamm*rotation;

Om = g2wFunc(B0, (Sys.g(1, 3) - Sys.g(2, 3))/2)^2 + (Sys.J/2 + Sys.dip/2)^2;
Om = sqrt(Om);
eigen1 = Om - (Sys.J/4 - Sys.dip/2)
eigen2 = - Om - (Sys.J/4 - Sys.dip/2)

load("vik.mat")
clf
imagesc(squeeze(reshapeThPh(ddn)) - dd)
colormap(vik)
colorbar

%%
% Theta, phi grid
nTheta = 91*2 - 1;
thetas = linspace(0, 180, nTheta)*pi/180;
% thetas = 10*pi/180;
nTheta = numel(thetas);
nPhi = 91*2*2 - 1;
phis = linspace(0, 360, nPhi)*pi/180;

% Direction of B0
clear('nVers')
nVers(1, :, :) = sin(thetas')*cos(phis);
nVers(2, :, :) = sin(thetas')*sin(phis);
nVers(3, :, :) = cos(thetas')*ones(1, nPhi);
% Effective g-values
g2 = squeeze(sqrt( sum( (Sys.g(2, :)'.*nVers).^2)));
% nVersInFrame1 = pagemtimes(erot(eAngles), nVers);
% g1 = squeeze(sqrt( sum( (Sys.g(1, :)'.*nVersInFrame1).^2)));
euMatrixg1 = erot(Sys.gFrame(1, :));
g1TensorInFrame2 = euMatrixg1'*diag(Sys.g(1, :))*euMatrixg1;
g1 = squeeze(sqrt( sum( (pagemtimes(g1TensorInFrame2, nVers)).^2)));
% Dipolar interaction
zD = erot(Sys.eeFrame)*[0, 0, 1]';
dd = squeeze(dipFunc(dip, nVers, zD));

%
% Calculate w0, deltaw, Omega
%
w0 = g2wFunc(B0, (g1 + g2)/2);  % MHz
deltaw = g2wFunc(B0, (g1 - g2)/2);  % MHz

%
% hfi 
%
if isfield(Sys, 'Nucs')
    nNuc = numel(strsplit(Sys.Nucs, ','));  % Number of equivalent nuclei
    nHfiLine = nNuc + 1;
else
    nNuc = 0;
    nHfiLine = 0;
end

%
% Hyperfine interaction with a number nNuc of equal nuclei
%

%
% hfi 
%
if isfield(Sys, 'Nucs')
    nNuc = numel(strsplit(Sys.Nucs, ','));  % Number of equivalent nuclei
    nHfiLine = nNuc + 1;
else
    nNuc = 0;
    nHfiLine = 0;
end

%
% Hyperfine interaction with a number nNuc of equal nuclei
%
if nNuc ~= 0
    % For spin 1 (for PSI, Sys.A(1, 1:3) = [0, 0, 0])
    euMatrixHfi1 = erot(Sys.AFrame(1, 1:3));
    Ahfi1TensorInFrame2 = euMatrixHfi1'*diag(Sys.A(1, 1:3))*euMatrixHfi1;
    Ahfi1 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi1TensorInFrame2', nVers)).^2)));
    % For spin 2
    euMatrixHfi2 = erot(Sys.AFrame(1, 4:6));
    Ahfi2TensorInFrame2 = euMatrixHfi2'*diag(Sys.A(1, 4:6))*euMatrixHfi2;
    Ahfi2 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi2TensorInFrame2', nVers)).^2)));

    % Shift resonances due to hyperfine interaction
    Aplus = (Ahfi1 + Ahfi2)/2;  % MHz
    % Make w0 have size [nHfiLine, nTheta, nPhi]qrt(deltaw.^2 + (J + d/2).^2);
    w0 = insertDimensionPos1(w0, nHfiLine);
    % Make Aplus have size [1, nTheta, nPhi]
    Aplus = reshape(Aplus, [1, size(Aplus)]);
    % Multiply by all possible values of sum of quantum number of the nuclei
    % Final size = [nHfiLine, nTheta, nPhi]
    Aplus = pagemtimes((-nNuc/2:1:nNuc/2)', Aplus);
    w0 = w0 + Aplus;
    % Same for deltaw and Aminus
    Aminus = (Ahfi1 - Ahfi2)/2;  % MHz
    deltaw = insertDimensionPos1(deltaw, nHfiLine); 
    Aminus = reshape(Aminus, [1, size(Aminus)]);
    Aminus = pagemtimes((-nNuc/2:1:nNuc/2)', Aminus);
    deltaw = deltaw + Aminus;

    % Pascal factors because some values of the sum of the quantum number
    % come up more than others
    pascalMatrix = pascal(nHfiLine);
    pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);  % Antidiag

else
    pascalFactor = 1;
end

Omega = OmegaFunc(deltaw, JJ, insertDimensionPos1(dd, nHfiLine));  % MHz

%
% Calculate the resonance positions in GHz
%
clear('wReson')
dipInteraction = (JJ - insertDimensionPos1(dd, nHfiLine));  % MHz

wReson(1, :, :, :) = w0 - dipInteraction - (Omega);  % w12
wReson(2, :, :, :) = w0 + dipInteraction - (Omega);  % w34
wReson(3, :, :, :) = w0 - dipInteraction + (Omega);  % w13
wReson(4, :, :, :) = w0 + dipInteraction + (Omega) ;  % w24

% BReson = ...
%     unitconvert(wReson*1e-6, 'MHz->mT', ...
%     permute(repmat((g1 + g2)/2, [1, 1, 4]), [3, 1, 2]));

% Lineshapes for each transition  - (xxSim(2) - xxSim(1))/8
trSignal = ...
    gaussiantransitions(xxSim', wReson*1e-3, trlw*1e-3, "var");
% trSignal = ...
%     gaussiantransitions(xxSim', wReson, trlw*1e-3, "fwhm");

% Intensities
intensityReson = 1/8*(deltaw.^2)./(Omega.^2);
% Sign of intensityReson alternates for the transitions w12, w34, w13, w24
intensityReson = [1; -1; 1; -1].*insertDimensionPos1(intensityReson, 4);  

% Normalization to account for possible different linewidths
intensityNorm = insertDimensionPos1(intensityReson, 1) ./ sum(trSignal);
trSignal = trSignal.*intensityNorm.*insertDimensionPos1(pascalFactor, 1);
if nNuc ~= 0
    trSignalSumHfi = squeeze(sum(trSignal, 3));  % Sum over hfi lines
    trSignalSum = squeeze(sum(trSignalSumHfi, 2));  % Sum over transitions
else
    trSignalSum = squeeze(sum(trSignal, 2));  % Sum over transitions
end

solidAngleWeight = sin(thetas)/sum(sin(thetas')*ones(1, nPhi), 'all');
trSignalPowder = sum(squeeze(sum(trSignalSum.*solidAngleWeight, 3)), 2);

% hold on
% clf
% % plot(xxSim1, ySim/max(ySim), xxSim1 - 0.21, flip(trSignalPowder/max(trSignalPowder)))
% plot(xxSim, ySim/max(ySim))
% hold on
% % plot(xxSim(1:end - 1), trSignalPowder(2:end)/max(trSignalPowder))
% plot(xxSim, trSignalPowder/max(trSignalPowder(:)))

% yyaxis right
% plot(xxSim, ySim'/max(ySim) - trSignalPowder/max(trSignalPowder))
% xlim([9.695, 9.708])
toc

%%
newsig = interp1(xxSim, trSignalPowder/max(trSignalPowder), xxSim + 0.00005);
% hold on
clf
% plot(xxSim1, ySim/max(ySim), xxSim1 - 0.21, flip(trSignalPowder/max(trSignalPowder)))
plot(xxSim, ySim/max(ySim))
hold on
% plot(xxSim(1:end - 1), trSignalPowder(2:end)/max(trSignalPowder))
% plot(xxSim, trSignalPowder/max(trSignalPowder(:)))
plot(xxSim, newsig)
yyaxis right
% plot(xxSim, ySim'/max(ySim) - trSignalPowder/max(trSignalPowder))
plot(xxSim, ySim/max(ySim) - newsig)
% xlim([9.695, 9.708])
toc

%% Plot

figure()
clf
plot(xxSim, ySim/max(ySim), xxSim, trSignalPowder/max(trSignalPowder))
ylim(setaxlim(trSignalPowder/max(trSignalPowder), 1.05))
yticks(0)
labelaxesfig(gca, 'Microwave frequency / GHz', 'trEPR signal')
legend('Easyspin', 'Gianluca')

% Using Sys.dip = unitconvert(+0.177,'mT->MHz')/1.5; % MHz
figPath = "/net/storage/gianlum33/projects/oop_ciss_calculations/images/";
figPath = figPath + "zech_psiOopEseem_trEPRstickSpectra_EasyspinComparison";
% exportgraphics(gcf, append(figPath, '.pdf'))

%% Reproduce experimental data

load('/net/storage/gianlum33/projects/zech_psi/data/processed/ZePSI-E-007015.mat');
% clf
% h = ScrollableAxes();
% plot(h, x{2}, x{1}, y');
xdata = x{2}/10;  % mT
ydata = y(500, :);
ydata = ydata/max(ydata);

%%

figure()
clf
plot(xdata, ydata, xxSim1-5.7, flip(trSignalPowder/max(trSignalPowder)))
xlim([min([xdata; xxSim1'-5.7]), max([xdata; xxSim1'-5.7])])
ylim(setaxlim(trSignalPowder/max(trSignalPowder), 1.05))
xticks(339:342)
yticks(0)
labelaxesfig(gca, 'Magnetic field / mT', 'trEPR signal')
legend('Exp. data', 'Gianluca')

% Using Sys.dip = unitconvert(+0.177,'mT->MHz')/1.5; % MHz
figPath = "/net/storage/gianlum33/projects/oop_ciss_calculations/images/";
figPath = figPath + "zech_psiOopEseem_trEPRstickSpectra_experimentComparison";
exportgraphics(gcf, append(figPath, '.pdf'))

%%

tan2alpha = squeeze(deltaw) ./ ((JJ + dd/2)*1e-3);

%%
%{
2.00300 0.00000 0.00000
0.00000 2.00260 0.00000
0.00000 0.00000 2.00230 G-TENSOR P700+

2.00620 0.00000 0.00000
0.00000 2.00510 0.00000
0.00000 0.00000 2.00220 G-TENSOR Qb-


tic
clear('Sys', 'Exp', 'Vary', 'VaryExp')

Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.J = unitconvert(-1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(+0.177,'mT->MHz'); % MHz
% Sys.lw = 9;  % MHz
Sys.lw = mhz2mt(10);  % MHz
Sys.initState = 'singlet';

nNuc = 3;
Sys.Nucs = '1H';
Sys.Nucs = repmat(append(Sys.Nucs, ', '), [1, nNuc]);
Sys.Nucs = Sys.Nucs(1:end - 2);
Sys.A = [zeros(1, 3), 9, 9, 12.8];
Sys.A = repmat(Sys.A, [nNuc, 1]);
Sys.AFrame = [0 0 0 -60 90 0]*pi/180;
Sys.AFrame = repmat(Sys.AFrame, [nNuc, 1]);

% Vary.dip = 1;
% Vary.J = 0.1;
Vary.lw = 5;
% VaryExp.mwCenterSweep = [0.1, 0];

B0 = 345.9;  % mT
mwFreq = 9.7;
xxSim = linspace(9.65, 9.75, 2048);  % GHz
xxSim1 = mhz2mt(xxSim*1e3);  % mT

% Magnetic field sweep
Exp.CenterSweep = [mean(xdata) + 1.95, max(xdata) - min(xdata)];
Exp.mwFreq = 9.6;
Exp.nPoints = numel(xdata);
Exp.Harmonic = 0;

VaryExp.CenterSweep = [0.5, 0];

FitOpt.Verbosity = 0;
FitOpt.maxTime = 5;

% 
% % Opt.GridSize = 10;
% FitOpt.Verbosity = 1;
% FitOpt.maxTime = 10;
% % FitOpt.delta = 1;

yFit0 = pepper(Sys, Exp);

plot(xdata, ydata, xdata, yFit0/max(yFit0))
toc

%%

esfit(ydata, @pepper, {Sys, Exp}, {Vary, VaryExp}, FitOpt)

%%
plot(xdata, ydata, xdata, Fit.fit)

%%
clf
plot(xxSim, ySim/max(ySim))
hold on
plot(xxSim, trSignalSum(:, 1, 1)/max(trSignalSum(:, 1, 1)))
% plot(xxSim, trSignal(:, 1, 1, 1)/max(trSignal(:, 1, 1, 1)))
% plot(xxSim, yLow(:, 1, 1, 1)/max(yLow(:, 1, 1, 1)))
% plot(xxSim, yHigh(:, 1, 1, 1)/max(yHigh(:, 1, 1, 1)))

%%
ith = 1;
iph = 1;
% clf
hold on
plot(xxSim, trSignalSum(:, ith, iph)/max(trSignalSum(:, ith, iph)))

%%
for itrans = 1:4
    plot(xxSim, trSignal(:, itrans, ith, iph))
    hold on
end



%%

yA = gaussiantransitions(xxSim', wReson, ...
    repmat(Sys.lw*1e-3, [4, nTheta, nPhi]), "fwhm");
yB = gaussiantransitions(xxSim', wReson, ...
    repmat(reshape(Atensor*1e-3, [1, nTheta, nPhi]), [4, 1, 1]), "fwhm");

yC = gaussiantransitionshfi(xxSim', wReson, ...
    repmat(Sys.lw*1e-3, [4, nTheta, nPhi]), ...
    repmat(reshape(Atensor*1e-3, [1, nTheta, nPhi]), [4, 1, 1]), "fwhm");

ith = 1;
iph = 1;
yplot1 = yA(:, itrans, ith, iph);
yplot1 = yplot1/max(yplot1);
yplot2 = yB(:, itrans, ith, iph);
yplot2 = yplot2/max(yplot2);
yplot3 = yC(:, itrans, ith, iph);
yplot3 = yplot3/max(yplot3);
clf
itrans = 1;
plot(xxSim, yplot1)
hold on
plot(xxSim, yplot2)

convSpectrum = conv(yplot1, yplot2, 'same');
plot(xxSim, convSpectrum/max(convSpectrum))
plot(xxSim - 0.0051, yplot3)
yline(0.5)

%% Reproduce experimental data

load('S:\projects\zech_psi\data\processed\ZePSI-E-007015.mat');
% clf
% h = ScrollableAxes();
% plot(h, x{2}, x{1}, y');
xdata = x{2}/10;  % mT
ydata = y(815, :);
ydata = ydata/max(ydata);

% Taken from Easyspin example
clear('Sys', 'Exp', 'Vary', 'VaryExp')
Sys.S = [1/2 1/2];
Sys.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.J = unitconvert(-1e-3,'mT->MHz'); % MHz
Sys.dip = unitconvert(0.177,'mT->MHz'); % MHz
Sys.lw = 0.3;
Sys.Nucs = 'H';
Sys.A = [0, 0, 0, 9, 9, 13];

Sys.initState = 'singlet';

% Vary.g = [0.001, 0.001, 0.001; 0.001, 0.001, 0.001];
Vary.lw = 5;
Vary.dip = 1;
Vary.J = 0.1;
Vary.A = [0, 0, 0, 9, 9, 13];


Exp.CenterSweep = [mean(xdata) + 2, max(xdata) - min(xdata)];
Exp.mwFreq = 9.6;
Exp.nPoints = numel(xdata);
Exp.Harmonic = 0;
VaryExp.CenterSweep = [3, 0];

% Opt.GridSize = 10;
FitOpt.Verbosity = 0;


clf
plot(xdata, ydata, xdata, pepper(Sys, Exp)/max(pepper(Sys, Exp)))

%%

esfit(ydata, @pepper, {Sys, Exp}, {Vary, VaryExp})

%%

figure(4)
s = surf(X, Y, Z, g2, 'EdgeColor', 'none');
axis equal
viridis = viridis();
colormap(viridis)
colorbar       
xlabel('X')
title('g2')

figure(6)
s = surf(X, Y, Z, gSys, 'EdgeColor', 'none');
axis equal
viridis = viridis();
colormap(viridis)
colorbar       
xlabel('X')
title('g1Sys')

figure(7)
s = surf(X, Y, Z, g11, 'EdgeColor', 'none');
axis equal
viridis = viridis();
colormap(viridis)
colorbar       
xlabel('X')
title('g11')

figure(8)
s = surf(X, Y, Z, g12, 'EdgeColor', 'none');
axis equal
viridis = viridis();
colormap(viridis)
colorbar       
xlabel('X')
title('g12')

nVersInFrame1 = pagemtimes(erot([81 126 50]*pi/180), nVers);
% nVersInFrame1 = nVers;
g1 = squeeze(sqrt( sum( (Sys.g(1, :)'.*nVersInFrame1).^2)));

figure(9)
s = surf(X, Y, Z, g1, 'EdgeColor', 'none');
axis equal
viridis = viridis();
colormap(viridis)
colorbar       
xlabel('X')
title('g1')

%}
%}
%%

function y1 = gaussianbs(wReson, mwFreq, c, mode)
    arguments
        wReson
        mwFreq  
        c    (1, 1) double     
        mode string = "lwpp"
    end

    if strcmp(mode, "fwhm")
        % If the input is fwhm, multiply the input parameters to have var
        c = c/(2*sqrt(2*log(2)));
    elseif strcmp(mode, "lwpp")
        c = c/2;
    elseif strcmp(mode, "var")
        % Do nothing
    else
        error("Indicate if c is 'lwpp', 'fwhm'(=lw) or" ...
            + "'var'(=variance) of the gaussian.")
    end
    
    % Adjust dimensions of xx if needed
    % [nAx, dim2] = size(xx);
    % if nAx < dim2
    %     xx = xx';
    %     nAx = dim2;
    % end
    % nTotLine = size(x0, 2);
    % 
    % xx = repmat(xx, [1, nTotLine]);
    % x0 = repmat(x0, [nAx, ones(1, nTotLine)]);
    % c = repmat(c, [nAx, ones(1, nTotLine)]);

    y1 = exp(-1/2 * (wReson - mwFreq).^2 ./ c^2);

end
function y1 = newlorentziantransitions(xx, x0, c, mode)
    % Output:
    % y1:
    arguments
        xx      % [nAx, 1] equally spaced values
        x0      
        c       
        mode string = "fwhm"
    end

    if strcmp(mode, "fwhm")
        % If the input is fwhm, multiply the input parameters to have var
        c = c/sqrt(3);
    end
    
    % Adjust dimensions of xx if needed
    [nAx, dim2] = size(xx);
    if nAx < dim2
        xx = xx';
        nAx = dim2;
    end
    nTotLine = size(x0, 2);

    xx = repmat(xx, [1, nTotLine]);
    x0 = repmat(x0, [nAx, ones(1, nTotLine)]);
    c = repmat(c, [nAx, ones(1, nTotLine)]);

    y1 = (1 + 4/3*(xx - x0).^2 ./ c.^2).^(-1);
end

function y1 = newgaussiantransitions(xx, x0, c, mode)
    % Output:
    % y1:
    arguments
        xx      % [nAx, 1] equally spaced values
        x0      
        c       
        mode string = "fwhm"
    end

    if strcmp(mode, "fwhm")
        % If the input is fwhm, multiply the input parameters to have var
        c = c/(2*sqrt(2*log(2)));
        disp('fwhm')
    end
    
    % Adjust dimensions of xx if needed
    [nAx, dim2] = size(xx);
    if nAx < dim2
        xx = xx';
        nAx = dim2;
    end
    nTotLine = size(x0, 2);

    xx = repmat(xx, [1, nTotLine]);
    x0 = repmat(x0, [nAx, ones(1, nTotLine)]);
    c = repmat(c, [nAx, ones(1, nTotLine)]);

    y1 = exp(-1/2 * (xx - x0).^2 ./ c.^2);
end

function y1 = gaussiantransitions(xx, x0, c, mode)
    % Output:
    % y1:
    arguments
        xx      % [nAx, 1] equally spaced values
        x0      % [4, nHfiLine, nTheta, nPhi]
        c       
        mode string = "fwhm"
    end
    
    if strcmp(mode, "fwhm")
        % If the input is fwhm, multiply the input parameters to have var
        c = c/(2*sqrt(2*log(2)));
    end
    
    % Adjust dimensions of xx if needed
    [nAx, dim2] = size(xx);
    if nAx < dim2
        xx = xx';
        nAx = dim2;
    end
    if min(size(xx)) ~= 1
        error('gaussiantransitions: xx should have size [n, 1] or [1, n].')
    end

    
    if numel(size(x0)) > 4  % nNuc ~= 0
        y1 = zeros([nAx, size(x0)]);  % Size [nAx, 4, nHfiLine, nTheta, nPhi]
        xx = repmat(xx, [1, size(x0, [1, 3, 4])]);  % Size [nAx, 4, nTheta, nPhi]
        for iNuc = 1:size(x0, 2)  % 1:(nHfiLine)
            x0_ = squeeze(x0(:, iNuc, :, :));
            x0_ = repmat( ...
                reshape(x0_, [1, size(x0_)]), ...
                [nAx, ones(1, numel(size(x0_)))]);
            xx_ = xx - x0_;
            y1(:, :, iNuc, :, :) = exp(-1/2 * (xx_).^2 ./ c.^2);
        end
    else
        xx = repmat(xx, [1, size(x0, [1, 2, 3])]);  % Size [nAx, 4, nTheta, nPhi]
        x0 = repmat( ...
                reshape(x0, [1, size(x0)]), ...
                [nAx, ones(1, numel(size(x0)))]);
        y1 = exp(-1/2 * (xx - x0).^2 ./ c.^2);
    end
    % y1 = exp(-1/2 * (xx - x0).^2 ./ c.^2);
    % Normalization over a range equal to 10*fwhm such that the area is
    % equal to 1
%     dxx = xx(2) - xx(1);
%     xxNorm = -5*max(fwhm(:)):dxx:5*max(fwhm(:));
%     xxNorm = repmat(xxNorm', size(x04D));
%     xxNorm = xxNorm + x04D;
%     y2 = exp(-4*log(2)*(xxNorm - x04D).^2 ./ fwhm4D.^2);
%     norm = sum(y2, 1)*dxx;
%     y1 = y1 ./ norm;

end


function vec = insertDimensionPos1(vec, n)
% Input vec of dimensions a1 x a2 x a3 x ... and return vec of dimensions 
% n x a1 x a2 x a3 x ... with vec(i, :, :, ...) = vec(j, :, :, ...) for
% every i and j
     if n ~= 0
         vec = repmat( ...
             reshape(vec, [1, size(vec)]), ...
             [n, ones(1, numel(size(vec)))]);
     end
end





