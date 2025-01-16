%
clearvars

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
Sys.trlwpp = 0.35;

Sys.rho = zeros(4, 4);
% Up-down and down-up
% Sys.rho(2, 2) = 1;
% Sys.rho(3, 3) = 1;
% Sys.rho = 1/2*Sys.rho;
% Singlet
Sys.rho(2, 2) = 1;
Sys.rho(2, 3) = -1;
Sys.rho(3, 2) = -1;
Sys.rho(3, 3) = 1;
Sys.rho = 1/2*Sys.rho;

Sys.mwFreq = 9.6;
Sys.x = linspace(340.5, 344, 301);  % mT
Sys.nTheta = 50;
Sys.nPhi = 20;

%%

test = treprstickspectrum(Sys);
aa0 = averageoversolidangle(test, Sys.nTheta, Sys.nPhi, 2);
aa0 = real(aa0);
% figure()
clf
hold on
box on
plot(Sys.x, aa0)

%%
function signal = treprstickspectrum(Sys)
    
    fieldAxis = Sys.x;
    nField = numel(fieldAxis);
    dip = Sys.dip;
    JJ = Sys.J;
    trlwpp = Sys.trlwpp;
    mwFreq = Sys.mwFreq;
    nTheta = Sys.nTheta;
    nPhi = Sys.nPhi;
    rhoInit = Sys.rho;  % This should be in the T+, S, T0, T- basis
    if isfield(Sys, 'nNuc')
        nNuc = Sys.nNuc;
    else
        nNuc = 0;
    end
    nHfiLine = nNuc + 1;

    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi);
    nSolidAngle = nTheta*nPhi;


    % Direction of B0
    nVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
    size(nVers)
    % Effective g-values
    g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
    g1 = sqrt( sum( (g1Tensor*nVers).^2, 1));
    g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
    g2 = sqrt( sum( (g2Tensor*nVers).^2, 1));
    % Dipolar interaction
    zD = erot(Sys.eeFrame)*[0, 0, 1]';
    dd = dipinteraction(dip, nVers, zD);

    % Hyperfine
    if nNuc > 0
        [APlus, AMinus] = calculateaplusaminusnohfi(...
            Sys.A, Sys.AFrame, nVers);
    else
        APlus = zeros(nSolidAngle, 1);
        AMinus = zeros(nSolidAngle, 1);
    end

    signal = zeros(nField, nSolidAngle);
    for ii = 1:nSolidAngle
        disp(floor(ii/nSolidAngle*100))
        gPlus = 1/2*(g1(ii) + g2(ii));
        gMinus = 1/2*(g1(ii) - g2(ii));

        w0_ = gvalue2freq(fieldAxis, gPlus);
        deltaw_ = gvalue2freq(fieldAxis, gMinus);

        quantNumNuc = -nNuc/2:1:nNuc/2;
        if nNuc > 0
            pascalMatrix = pascal(nHfiLine);
            % Antidiag
            pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
            pascalFactor = pascalFactor'/sum(pascalFactor);
        else
            pascalFactor = 1;
        end

        for itrans = 1:4
            for ihfi = 1:nHfiLine
                w0__ = w0_ + quantNumNuc(ihfi)*APlus(ihfi);
                deltaw__ = deltaw_ + quantNumNuc(ihfi)*AMinus(ihfi);
                wReson = myeigenenergies(w0__, deltaw__, JJ, dd(ii), itrans);
                intensityReson = intensityresonance(rhoInit, nVers(:, ii), deltaw__, JJ, dd(ii), itrans);
                signal__ = gaussianresonancebsweep( ...
                    wReson*1e-3, mwFreq, mt2mhz(trlwpp)*1e-3, "lwpp");
                signal(:, ii) = signal(:, ii) + intensityReson'.*signal__'*pascalFactor(ihfi);

            end 
        end
    end

end

