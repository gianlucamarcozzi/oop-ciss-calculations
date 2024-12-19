%
clearvars
addpath(genpath('/net/storage/gianlum33/projects/oop_ciss_calculations/util/'))
addpath(genpath('S:\projects/oop_ciss_calculations/util/'))

%%

% Import
% cd('/net/storage/gianlum33/projects/oop_ciss_calculations/')
expName = 'data/digitized/zech_p46_oopEseem_expData.csv';
fitName = 'data/digitized/zech_p46_oopEseem_fit.csv';

importedData = readtable(expName);
x = importedData{:, 1}*pi/180;
y = importedData{:, 2};
y = rescaledata(y, 'maxabs');

importedData = readtable(fitName);
xx = importedData{:, 1}*pi/180;
yy = importedData{:, 2};
yy = rescaledata(yy, 'maxabs');

crystalEseemZech = @(p, beta) -( ...
    0.5*sin(2*beta) .* sin(2*p).^4 + ...
    sin(beta) .* cos(2*p).^2 .* sin(2*p).^2);  % p = alpha

p0Z = 31.5*pi/180;
p0Timmel = 41*pi/180;
initFit = rescaledata(crystalEseemZech(p0Z, xx), 'maxabs');
expectFit = rescaledata(crystalEseemZech(p0Timmel, xx), 'maxabs');

p0B = pi/2 - 2*p0Z;  % xi_Bittl = pi/2 - 2*alpha_Zech
initFit2 = rescaledata(crystalEseem(p0B, xx), 'maxabs');

clf
plot(x*180/pi, y, 'o-')
hold on
plot(xx*180/pi, yy)
plot(xx*180/pi, initFit)
plot(xx*180/pi, initFit2)
plot(xx*180/pi, expectFit)
% esfit()

%% System definition

Sys0.S = [1/2 1/2];
Sys0.g = [2.0030 2.0026 2.0023; ... % P700+
         2.0062 2.0051 2.0022];    % A1-
Sys0.gFrame = [-10 -128 -83; ...
                0    0   0]*pi/180;
Sys0.eeFrame = [0 90 0]*pi/180;  % zD direction is along -x of A1-

Sys.g = Sys0.g;
Sys.gFrame = Sys0.gFrame;
Sys.dip = mt2mhz(-0.170);  % MHz
Sys.J = mt2mhz(1e-3);  % MHz
Sys.eeFrame = Sys0.eeFrame;
Sys.nNuc = 0;
Sys.A = [0, 0, 0, 9, 9.4, 12.8];  % MHz
Sys.AFrame = [0 0 0 60 -90 0]*pi/180;

Sys.mwFreq = 9.7;
Sys.x = linspace(344.5, 347, 301);  % mT
Sys.firstPulseTurningAngles = xx;
Sys.nTheta = 90; %, 50};
Sys.nPhi = 45; %, 20};
Sys.lw1 = 0.5;  % mT
Sys.lw2 = 0.5;  % mT
Sys.nSamplegInh1 = 100;
Sys.nSamplegInh2 = 100;

%%

[a0, b0] = oopeseemvsbeta(Sys, 0);
tic
[a2, b2] = oopeseemvsbeta(Sys, 1);
toc

%%

aa0 = averageoversolidangle(a0, Sys.nTheta, Sys.nPhi, 2);
aa1 = averageoversolidangle(a1, Sys.nTheta, Sys.nPhi, 2);
% aa2 = averageoversolidangle(a2, Sys.nTheta, Sys.nPhi, 2);
bb0 = averageoversolidangle(abs(b0), Sys.nTheta, Sys.nPhi, 2);
bb1 = averageoversolidangle(abs(b1), Sys.nTheta, Sys.nPhi, 2);
% bb2 = averageoversolidangle(abs(b2), Sys.nTheta, Sys.nPhi, 2);
% sum(abs(a1(:) - a2(:)))
% sum(abs(b1 - b2))
clf
plot(xx*180/pi, rescaledata(aa0, 'maxabs'))
hold on
plot(xx*180/pi, rescaledata(aa1, 'maxabs'), 'o')
% plot(xx*180/pi, rescaledata(aa2, 'maxabs'), 'x')
disp("xiMean0 = " + string(bb0*180/pi) + newline + ...
    "xiMean1 = " + string(bb1*180/pi) + newline);% + ...
%     "xiMean2 = " + string(bb2*180/pi))

%%

wb = waitbar(0, 'Hi', 'Name', string(datetime));
nSim = 1001;
Sys.nTheta = 40;
Sys.nPhi = 20;
samplings = linspace(10, 1000, nSim);
samplings = round(samplings);
tic
clear oop xi oopPowder xiPlot
xiMean = zeros(1, nSim);
for isim = 1:nSim
    waitbar(isim/nSim, wb, sprintf('%d out of %d', isim, nSim))
    Sys.nSamplegInh1 = samplings(isim);
    Sys.nSamplegInh2 = samplings(isim);
    [oop{isim}, xi{isim}] = oopeseemvsbeta(Sys);
    xiMean(isim) = averageoversolidangle(...
        abs(xi{isim}), Sys.nTheta, Sys.nPhi, 2);
    oopPowder{isim} = averageoversolidangle(...
        oop{isim}, Sys.nTheta, Sys.nPhi, 2);
    xiPlot{isim} = reshape(xi{isim}, [Sys.nTheta, Sys.nPhi])*180/pi;
%     if isim == 1  % Pre-allocation
%         oop = repmat(oop, 1, nSim);
%         xi = repmat(xi, 1, nSim);
%         oopPowder = repmat(oopPowder, 1, nSim);
%         xiPlot = repmat(xiPlot, 1, nSim);
%     end  
end
toc
delete(wb); clear('wb');

%%

clf
plot(x*180/pi, y, 'o-')
% plot(nThetas, xiMeana*180/pi, 'o')
hold on
plot(samplings, xiMean*180/pi, 'o')
clf
plot(xx*180/pi, oopPowder{1} - oopPowder{end})
hold on

%%

% tau = 0.;  % us

% Assume finite bandwidth of the excitation pulse
pulseBw = 0.50;  % mT
pulseFieldPosition = 346.1; % mean(Sys.x);
Sys.excitRange = [-1/2, +1/2]*pulseBw + pulseFieldPosition;

Sys.nTheta = 40;
Sys.nPhi = 30;
[oop, xixi, g1, g2, idxAngle, idxAngle1, idxAngle2] = oopeseemvsbeta(Sys, 0, 1);

figure(1)
g1 = reshape(g1, [Sys.nTheta, Sys.nPhi]);
clf
imagesc(g1)
colorbar
colormap(viridis())
hold on
for ii = 1:numel(idxAngle1)
    [xSquare, ySquare] = drawsquaresexcitrange(idxAngle1(ii), Sys.nTheta);
    plot(xSquare, ySquare, 'k', 'LineWidth', 2)
end
for ii = 1:numel(idxAngle)
    [xSquare, ySquare] = drawsquaresexcitrange(idxAngle(ii), Sys.nTheta);
    plot(xSquare, ySquare, 'r', 'LineWidth', 2, 'LineStyle', '--')
end

figure(2)
g2 = reshape(g2, [Sys.nTheta, Sys.nPhi]);
clf
imagesc(g2)
colorbar
colormap(viridis())
hold on
for ii = 1:numel(idxAngle2)
    [xSquare, ySquare] = drawsquaresexcitrange(idxAngle2(ii), Sys.nTheta);
    plot(xSquare, ySquare, 'k', 'LineWidth', 2)
end
for ii = 1:numel(idxAngle)
    [xSquare, ySquare] = drawsquaresexcitrange(idxAngle(ii), Sys.nTheta);
    plot(xSquare, ySquare, 'r', 'LineWidth', 2, 'LineStyle', '--')
end

figure(3)
clf
xixiplot = reshape(abs(xixi), [Sys.nTheta, Sys.nPhi]);
imagesc(xixiplot)
colorbar
colormap(viridis())
hold on
for ii = 1:numel(idxAngle)
    [xSquare, ySquare] = drawsquaresexcitrange(idxAngle(ii), Sys.nTheta);
    plot(xSquare, ySquare, 'r', 'LineWidth', 2, 'LineStyle', '--')
end

%%

xiMeanInfBw = averageoversolidanglemask(abs(xixi), Sys.nTheta, Sys.nPhi, 2, []);
xiMeanBw = averageoversolidanglemask(abs(xixi), Sys.nTheta, Sys.nPhi, 2, idxAngle);

oopInfBw = averageoversolidanglemask(oop, Sys.nTheta, Sys.nPhi, 2, []);
oopBw = averageoversolidanglemask(oop, Sys.nTheta, Sys.nPhi, 2, idxAngle);

figure(4)
clf
plot(x*180/pi, y, 'o')
hold on
plot(xx*180/pi, yy, '--')
plot(xx*180/pi, rescaledata(oopInfBw, 'maxabs'))
plot(xx*180/pi, rescaledata(oopBw, 'maxabs'))
yline(0, 'Color', 'yellow')
xline(90, 'Color', 'yellow')

yyInterp = interp1(xx, yy, x);
oopInterp = interp1(xx, rescaledata(oopBw, 'maxabs'), x);
plot(x*180/pi, yyInterp, 'x')
plot(x*180/pi, oopInterp, 'x')
rmsdZech = sum( (y - yyInterp).^2)/numel(y);
rmsdMe = sum( (y - oopInterp).^2)/numel(y);

%%

pulseBw = 10;  % mT
pulseFieldPosition = 346.1; % mean(Sys.x);
Sys.excitRange = [-1/2, +1/2]*pulseBw + pulseFieldPosition;
gExcitRange = planck*Sys.mwFreq/bmagn./Sys.excitRange*1e12;
gExcitRange = flip(gExcitRange);

Sys.nTheta = 40;
Sys.nPhi = 30;
tic
[oopNoLw, xiNoLw, ~] = oopeseemvsbeta(Sys, 0, 0);
toc
tic
[oopLw, xiLw, ~] = oopeseemvsbeta(Sys, 1, 0);
toc
tic
[oopBw, xiBw, idxAngle0] = oopeseemvsbeta(Sys, 0, 1);
toc
tic
[oopLwBw, xiLwBw, ~] = oopeseemvsbeta(Sys, 1, 1);
toc

xiMeanNoLw = averageoversolidanglemask(abs(xiNoLw), Sys.nTheta, Sys.nPhi, 2, []);
xiMeanLw = averageoversolidanglemask(abs(xiLw), Sys.nTheta, Sys.nPhi, 2, []);
xiMeanBw = averageoversolidanglemask(abs(xiBw), Sys.nTheta, Sys.nPhi, 2, idxAngle0);
xiMeanLwBw = averageoversolidanglemask(abs(xiLwBw), Sys.nTheta, Sys.nPhi, 2, []);
oopPowderNoLw = averageoversolidanglemask(oopNoLw, Sys.nTheta, Sys.nPhi, 2, []);
oopPowderLw = averageoversolidanglemask(oopBw, Sys.nTheta, Sys.nPhi, 2, idxAngle);
oopPowderBw = averageoversolidanglemask(oopLw, Sys.nTheta, Sys.nPhi, 2, []);
oopPowderLwBw = averageoversolidanglemask(oopLwBw, Sys.nTheta, Sys.nPhi, 2, idxAngle);

xiNoLw = reshape(xiNoLw, [Sys.nTheta, Sys.nPhi]);
xiBw = reshape(xiBw, [Sys.nTheta, Sys.nPhi]);
xiLw = reshape(xiLw, [Sys.nTheta, Sys.nPhi]);
xiLwBw = reshape(xiLwBw, [Sys.nTheta, Sys.nPhi]);
figure(1)
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
imagesc(abs(xiNoLw)); colorbar; colormap(viridis())
title('xiNoLw')
nexttile
imagesc(abs(xiBw)); colorbar; colormap(viridis())
title('xiBw')
nexttile
imagesc(abs(xiLw)); colorbar; colormap(viridis())
title('xiLw')
nexttile
imagesc(abs(xiLwBw)); colorbar; colormap(viridis())
title('xiBwLw')

%%

% function [gax1, gInh1, gax2, gInh2, gMinusAx, savexi, gaussianFactor, xiMultiplicity] = oopEseemVsBeta(Sys, isInhBroadening)
% function [signal, xi, gax1, gax2, spinPairInRange, gaussianFactor, xiMultiplicity, xiMultiplicity2, gMinusAx] = oopeseemvsbeta(Sys, isInhBroadening, isFiniteExcitBw)
function [signal, xi, idxAngleInExcit] = oopeseemvsbeta(Sys, isInhBroadening, isFiniteExcitBw)
    arguments
        Sys
        isInhBroadening = 0
        isFiniteExcitBw = 0
    end

    Sys.isInhBroadening = isInhBroadening;
    Sys.isFiniteExcitBw = isFiniteExcitBw;
    idxAngleInExcit = [];
    idxAngleInExcit1 = [];
    idxAngleInExcit2 = [];

    errorMsg = checkfieldssys(Sys);
    if ~isempty(errorMsg)
        error(errorMsg)
    end
    
    fieldAxis = Sys.x;
    Bmean = mean(fieldAxis);
    firstPulseTurningAngles = Sys.firstPulseTurningAngles;
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
    if isFiniteExcitBw
        gExcitRange = planck*mwFreq/bmagn./Sys.excitRange*1e12;
        gExcitRange = flip(gExcitRange);  % Increasing value
    end
    
    % Parameters for g broadening
    lw1 = mt2mhz(Sys.lw1);
    lw2 = mt2mhz(Sys.lw2);
    % The g-factor linewidth
    glw1 = freq2gvalue(Bmean, lw1);
    glw2 = freq2gvalue(Bmean, lw2);
    % Spacing between g-values
    dgax1 = glw1 / Sys.nSamplegInh1;
    dgax2 = glw2 / Sys.nSamplegInh2;
    % resLw = min(lw1, lw2)*Sys.coeffResLw*1e-3;  % GHz
    
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
    
    signal = zeros(numel(firstPulseTurningAngles), nSolidAngle);
    xi = zeros(1, nSolidAngle);

    if ~isInhBroadening
        for ii = 1:nSolidAngle
            if isFiniteExcitBw
                g1flag = isginexcitrange(g1(ii), gExcitRange);
                g2flag = isginexcitrange(g2(ii), gExcitRange);
                if g1flag
                    idxAngleInExcit1(end + 1) = ii;
                end
                if g2flag
                    idxAngleInExcit2(end + 1) = ii;
                end
%                 if ~g1flag || ~g2flag
%                     continue
%                 end
                if g1flag && g2flag
                    idxAngleInExcit(end + 1) = ii;
                end
            end
            gPlus = 1/2*(g1(ii) + g2(ii));
            gMinus = 1/2*(g1(ii) - g2(ii));
    
            w0_ = gvalue2freq(Bmean, gPlus);
            deltaw_ = gvalue2freq(Bmean, gMinus);
    
            quantNumNuc = -nNuc/2:1:nNuc/2;
            if nNuc > 0
                pascalMatrix = pascal(nHfiLine);
                % Antidiag
                pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
                pascalFactor = pascalFactor'/sum(pascalFactor);
            else
                pascalFactor = 1;
            end
    
            signal__ = zeros(numel(firstPulseTurningAngles), nHfiLine);
            xi__ = zeros(1, nHfiLine);
            for ihfi = 1:nHfiLine
                    % w0__ = w0_ + quantNumNuc(ihfi)*APlus(ii);
                    deltaw__ = deltaw_ + quantNumNuc(ihfi)*AMinus(ii);
                    % savedeltaw(ii, ihfi) = deltaw__;
    
                    xi__(ihfi) = mixinganglexi(deltaw__, JJ, dd(ii));
                    % savexi(ii, ihfi) = xi__(ihfi);
                    signal__(:, ihfi) = ...
                        crystalEseem(xi__(ihfi), firstPulseTurningAngles);
            end
            
            % Average over hfi probability
            xi(ii) = xi__*pascalFactor;
            signal(:, ii) = signal__*pascalFactor;
        end
    else
        for ii = 1:nSolidAngle
            nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
            
            [gax1, gInh1] = creategaxisinhbroadening(g1(ii), glw1, dgax1);
            [gax2, gInh2] = creategaxisinhbroadening(g2(ii), glw2, dgax2);
            ngax1 = numel(gax1);
            ngax2 = numel(gax2);
            
            % It is not necessary to calculate every possible combination
            % but 
            gMinusAx = [flip((gax1 - gax2(1))/2), ...
                (gax1(1) - gax2(2:end))/2];
            ngMinusAx = ngax1 + ngax2 - 1;
            
%             w0_ = gvalue2freq(Bmean, gPlus);
            deltaw_ = gvalue2freq(Bmean, gMinusAx);
            % sintheta_ = sqrt( nVers(1, ii)^2 + nVers(2, ii)^2 );
    
            gaussianFactor = gInh1'*gInh2;
            gaussianFactor = gaussianFactor/sum(gaussianFactor, 'all');
            
            if isFiniteExcitBw
                flag1 = isginexcitrange(gax1, gExcitRange);
                flag2 = isginexcitrange(gax2, gExcitRange);
                spinPairInRange = flag1'*flag2;
            else
                spinPairInRange = ones(ngax1, ngax2);
            end
            diagsToExtract = -(ngax1 - 1):(ngax2 - 1);  % All the diagonals
            xiMultiplicity = sum(spdiags(...
                gaussianFactor.*spinPairInRange, diagsToExtract), 1);
            xiMultiplici
            
            quantNumNuc = -nNuc/2:1:nNuc/2;
            if nNuc > 0
                pascalMatrix = pascal(nHfiLine);
                % Antidiag
                pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
            else
                pascalFactor = 1;
            end
            
            signal__ = zeros(numel(firstPulseTurningAngles), nHfiLine, ngMinusAx);
            xi__ = zeros(1, nHfiLine, ngMinusAx);
            for ihfi = 1:nHfiLine
                % w0__ = w0_ + quantNumNuc(ihfi)*APlus(ii);
                deltaw__ = deltaw_ + quantNumNuc(ihfi)*AMinus(ii);

                xi__(1, ihfi, :) = mixinganglexi(deltaw__, JJ, dd(ii));
%                 savexi = xi__(ihfi, :);
                signal__(:, ihfi, :) = ...
                    crystalEseem(xi__(1, ihfi, :), firstPulseTurningAngles);
            end
            xi(ii) = squeeze(sum(...
                xi__.*reshape(xiMultiplicity, [1, 1, ngMinusAx]), ...
                3 ...
                ))*pascalFactor;
            signal(:, ii) = squeeze(sum(...
                signal__.*reshape(xiMultiplicity, [1, 1, ngMinusAx]), ...
                3 ...
                ))*pascalFactor;
        end
    end
end

function signal = crystalEseem(p, beta)
% p = xi = 2alpha
% beta = xx is the turning angle axis
signal = - (1/2 * sin(2*beta) .* cos(p).^4 + ...
            1/4 * sin(beta) .* sin(2*p).^2);
end

function errorMsg = checkfieldssys(Sys)

% excitRange
errorMsgContent = "excitRange should be a vector of " + ...
    "two elements of increasing value.";
errorMsg = '';
if Sys.isFiniteExcitBw
    if isempty(Sys.excitRange)
        errorMsg = errorMsgContent;
    elseif ~isa(Sys.excitRange, 'double')
        errorMsg = errorMsgContent;
    else
        if size(Sys.excitRange) ~= 2
            errorMsg = errorMsgContent;
        elseif Sys.excitRange(1) >= Sys.excitRange(2)
            errorMsg = errorMsgContent;
        end
    end
end

end

function flag = isginexcitrange(g, gExcitRange)
flag = ones(size(g));
idxOutRange = g < gExcitRange(1) | g > gExcitRange(2);
flag(idxOutRange) = 0;
% if g < gExcitRange(1) || g > gExcitRange(2)
%     flag = 0;  % g outside the range
% else
%     flag = 1;
% end
end

function [xSquare, ySquare] = drawsquaresexcitrange(idxAngle, nTheta)
xSquare0 = floor(idxAngle/nTheta) + 0.5;
if rem(idxAngle, nTheta) == 0
    xSquare0 = xSquare0 - 1;
end
ySquare0 = mod(idxAngle, nTheta) - 0.5;
if ySquare0 == -0.5
    ySquare0 = ySquare0 + nTheta;
end
xSquare = [xSquare0, xSquare0, xSquare0 + 1, xSquare0 + 1, xSquare0];
ySquare = [ySquare0, ySquare0 + 1, ySquare0 + 1, ySquare0, ySquare0];
end