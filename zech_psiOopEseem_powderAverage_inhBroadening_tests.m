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

%% Huge for loop - zero optimization

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
dgax = 0.00001;

thetas = 0:0.5:180;
nTheta = numel(thetas);
phis = 0:0.5:360;
nPhi = numel(phis);
[g1, g2, dd, valueOfCosThetaD] = deal(zeros(nTheta, nPhi));
[nVers, nVersInFrame1, xi, signal_] = deal({});
signalPowder = zeros(numel(xx), 1);


% wbarPh = waitbar(0, '1', 'Name', 'Phi', ...
%     'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');

tic
for ith = 1:1
    % Update waitbar and message
    if ~exist('wbarTh', 'var')
        wbarTh = waitbar(0, '1', 'Name', 'Theta', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    end
    waitbar(ith/nTheta, wbarTh, sprintf('%d out of %d', ith, nTheta))

    theta_ = thetas(ith)*pi/180;
    % theta_ = theta_/multTheta;
    % theta_ = theta_ * 10;
    for iph = 1:1
        % waitbar(iph/nPhi, wbarPh, sprintf('%d out of %d', iph, nPhi))
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

        % Consider inhomogenous broadening
        % Create g-axis in order to 'sample' the inh. broadened g-values
        gax = creategaxis(g1(ith, iph), g2(ith, iph), glw1, glw2, dgax);
        ngax = numel(gax);
        % Inhomogeneously broadened gaussian distributions
        gInh1 = gaussian(gax, g1(ith, iph), glw1);
        gInh2 = gaussian(gax, g2(ith, iph), glw2);
        % Normalization (?)
        gInh1 = gInh1/sum(gInh1);
        gInh2 = gInh2/sum(gInh2);
        
        xi_ = zeros(ngax, ngax);
        % xi__ = zeros(1, ngax);
        dd_ = dd(ith, iph);
        signalInh = zeros(numel(xx), ngax, ngax);
        signalInh_test = zeros(numel(xx), ngax, ngax);
        
        gaussianFactors = gInh1'*gInh2;
        for ig1 = 1:ngax
            %  g1_ = gax(ig1);
            for ig2 = 1:ngax
                %  g2_ = gax(ig2);
                xi_(ig1, ig2) = atan( (dd_ + 2*JJ)*1e6 * ...
                    planck/(2*pi*(bmagn*0.35)) / (gax(ig1) - gax(ig2)) );
%                 
                signalInh(:, ig1, ig2) = ...
                    gInh1(ig1)*gInh2(ig2)*crystalEseem(xi_(ig1, ig2), xx);
%                 signalInh(:, ig1, ig2) = ...
%                     crystalEseem(xi_(ig1, ig2), xx);
%               
                signalInh_test(:, ig1, ig2) = ...
                    gaussianFactors(ig1, ig2)*crystalEseem(xi_(ig1, ig2), xx);
            end
        end
        % xi(ith, iph) = atan( (dd(ith, iph) + 2*JJ)*1e6 * ...
        %     planck/(2*pi*(bmagn*0.35)) / (g1(ith, iph) - g2(ith, iph)) );
        
        signalSum = zeros(numel(xx), 1);
        for ig1 = 1:ngax - 2
            signalSum = signalSum + signalInh(:, ig1 + 2, ig1);
            signalSum = signalSum + signalInh(:, ig1, ig1 + 2);
        end
%         
%         xi{ith, iph} = xi_;
        signal_{ith, iph} = sum(sum(signalInh, 2), 3);
% 
    end

%     signalPowder = signalPowder + signal_{ith, iph}*sin(theta_);

end
toc

% delete(wbarTh);
% delete(wbarg1)

%% Some plots

ith = 1;
iph = 1;
fig = figure(6);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'tight')
nexttile
contour(xi{ith, iph})
colorbar

nexttile
plot(1:numel(xi{ith, iph}(:, 100)), xi{ith, iph}(:, 100))
textStr = sprintf('Max = %.3f, Min = %.3f', ...
    max(xi{ith, iph}(:, 100)), min(xi{ith, iph}(:, 100)));
text(10, -1.5, textStr)

nexttile
gax = creategaxis(g1(ith, iph), g2(ith, iph), glw1, glw2, dgax);
gInh1 = gaussian(gax, g1(ith, iph), glw1);
gInh2 = gaussian(gax, g2(ith, iph), glw2);
gInh1 = gInh1/sum(gInh1);
gInh2 = gInh2/sum(gInh2);
contour(gInh1'*gInh2)
colorbar

nexttile
plot(gax, gInh1, '-o')
hold on
plot(gax, gInh2, '-o')
yyaxis right
plot(gax, gInh1.*gInh2, '-o')

%{
nexttile(5, [2, 2])
plot(x*180/pi, rescaledata(y, 'maxabs'), 'o-')
hold on
plot(xx*180/pi, rescaledata(yy, 'maxabs'))
plot(xx*180/pi, initFit)
plot(xx*180/pi, initFit2)
plot(xx*180/pi, rescaledata(signalPowder, 'maxabs'), 'ko-')
textStr = sprintf('Max = %.3f, Min = %.3f', ...
    max(rescaledata(signalPowder, 'maxabs')), ...
    min(rescaledata(signalPowder, 'maxabs')));
text(10, 0.5, textStr)

textStr = "thetas = 0:0.5:180; phis = 0:0.5:360;" + newline + ...
    "dgax = 0.00005; computational time 5.16 h";
text(110, -0.5, textStr)
%}

%% Optimized for one-dimensional for loop

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
dgax = 0.00001;

thetas = 0:0.5:180;
nTheta = numel(thetas);
phis = 0:0.5:360;
nPhi = numel(phis);
[g1, g2, dd, valueOfCosThetaD] = deal(zeros(nTheta, nPhi));
[nVers, nVersInFrame1, xi, signal_] = deal({});
signalPowder = zeros(numel(xx), 1);


% wbarPh = waitbar(0, '1', 'Name', 'Phi', ...
%     'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');

tic
for ith = 1:1
    % Update waitbar and message
    if ~exist('wbarTh', 'var')
        wbarTh = waitbar(0, '1', 'Name', 'Theta', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    end
    waitbar(ith/nTheta, wbarTh, sprintf('%d out of %d', ith, nTheta))

    theta_ = thetas(ith)*pi/180;
    % theta_ = theta_/multTheta;
    % theta_ = theta_ * 10;
    for iph = 1:1
        % waitbar(iph/nPhi, wbarPh, sprintf('%d out of %d', iph, nPhi))
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

        % Consider inhomogenous broadening
        % Create g-axis in order to 'sample' the inh. broadened g-values
        gax = creategaxis(g1(ith, iph), g2(ith, iph), glw1, glw2, dgax);
        ngax = numel(gax);
        % Inhomogeneously broadened gaussian distributions
        gInh1 = gaussian(gax, g1(ith, iph), glw1);
        gInh2 = gaussian(gax, g2(ith, iph), glw2);
        % Normalization (?)
        gInh1 = gInh1/sum(gInh1);
        gInh2 = gInh2/sum(gInh2);
        
        xi_2 = zeros(ngax, ngax);
        % xi__ = zeros(1, ngax);
        dd_ = dd(ith, iph);
        signalInh_2 = zeros(numel(xx), ngax, ngax);
        
        % Calculate xi at fixed g2 value
        ig2 = 1;
        for ig1 = 1:ngax
            %  g2_ = gax(ig2);
            xi_2(ig1, ig2) = atan( (dd_ + 2*JJ)*1e6 * ...
                planck/(2*pi*(bmagn*0.35)) / (gax(ig1) - gax(ig2)) );
            signalInh_2(:, ig1, ig2) = crystalEseem(xi_2(ig1, ig2), xx);
        end
        
        for ig2 = 2:ngax
            xi_2(:, ig2) = [-flip(xi_2(2:ig2, 1)); ...
                xi_2(1:end + 1 - ig2, 1)];
            signalInh_2(:, :, ig2) = [flip(signalInh_2(:, 2:ig2, 1), 2), ...
                signalInh_2(:, 1:end - ig2 + 1, 1)];
        end
        signal_2{ith, iph} = sum(sum(signalInh_2, 2), 3);

    end

%     signalPowder = signalPowder + signal_{ith, iph}*sin(theta_);

end
toc

%% Check if signalInh is the same

i1 = 7;
i2 = 6;

yOld = signalInh(:, i1, i2);
y_2 =  signalInh_2(:, i1, i2);
% y_2 = signalInh(:, 3, 2);

clf
plot(xx, yOld, '-o')
hold on
plot(xx, y_2, '-o')
plot(xx, yOld - y_2, 'o')

% iax = 26;
% clf
% plot(1:ngax, xi_(:, iax), '-o')
% hold on
% plot(1:ngax, xi_2(:, iax), '-o')
% % plot(1:ngax, xi_(:, iax) - xi_2(:, iax))
% % plot(xx, flip(xi_))
% xlim([0, 40])

for i1 = 1:ngax - 1
    for i2 = 1:ngax - 1
        if xi_(i1, i2) - xi_(i1 + 1, i2 + 1) > 1e-5
            warning('%d %d', i1, i2)
        end
    end
end
    

%%

tic
for ith = 1:1
    % Update waitbar and message
    if ~exist('wbarTh', 'var')
        wbarTh = waitbar(0, '1', 'Name', 'Theta', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    end
    waitbar(ith/nTheta, wbarTh, sprintf('%d out of %d', ith, nTheta))

    theta_ = thetas(ith)*pi/180;
    % theta_ = theta_/multTheta;
    % theta_ = theta_ * 10;
    for iph = 1:1
        % waitbar(iph/nPhi, wbarPh, sprintf('%d out of %d', iph, nPhi))
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

        % Consider inhomogenous broadening
        % Create g-axis in order to 'sample' the inh. broadened g-values
        gax = creategaxis(g1(ith, iph), g2(ith, iph), glw1, glw2, dgax);
        ngax = numel(gax);
        % Inhomogeneously broadened gaussian distributions
        gInh1 = gaussian(gax, g1(ith, iph), glw1);
        gInh2 = gaussian(gax, g2(ith, iph), glw2);
        % Normalization (?)
        gInh1 = gInh1/sum(gInh1);
        gInh2 = gInh2/sum(gInh2);
        
        xi_ = zeros(ngax, 1);
        % xi__ = zeros(1, ngax);
        dd_ = dd(ith, iph);
        signalInh_3 = zeros(numel(xx), ngax);
        
        % Calculate xi at fixed g2 value
        for ig = 1:ngax
            %  g2_ = gax(ig2);
            xi_(ig) = atan( (dd_ + 2*JJ)*1e6 * ...
                planck/(2*pi*(bmagn*0.35)) / (gax(ig) - gax(1)) );
            signalInh_3(:, ig) = crystalEseem(xi_(ig), xx);
        end
        
        gaussianFactors = gInh1'*gInh2;
        signalInhFactors = zeros(1, ngax);
        signalInhFactors(1) = sum(spdiags(gaussianFactors, 0));
        for ig = 1:ngax - 1
            signalInhFactors(ig + 1) = sum(spdiags(gaussianFactors, -ig)) + ...
                sum(spdiags(gaussianFactors, ig));
        %     signalInhFactor2(ig) = sum(spdiags(gaussianFactor, -ig + 1))*2;
        end
        
        signalInh_3 = signalInh_3.*signalInhFactors;
        % Replicate for different values of g1
        
        % xi(ith, iph) = atan( (dd(ith, iph) + 2*JJ)*1e6 * ...
        %     planck/(2*pi*(bmagn*0.35)) / (g1(ith, iph) - g2(ith, iph)) );
        
%         
%         xi{ith, iph} = xi_;
        signal_3{ith, iph} = sum(signalInh_3, 2);
% 
    end

%     signalPowder = signalPowder + signal_{ith, iph}*sin(theta_);

end
toc

%%
matrixTest = ones(10, 10);
sumFactors = 0;
for ii = 1:ngax - 2
    sumFactors = sumFactors + gaussianFactors(ii, ii + 2);
    sumFactors = sumFactors + gaussianFactors(ii + 2, ii);
%     sumFactors = sumFactors + gaussianFactors(ii, ii);
%     sumFactors = sumFactors + matrixTest(ii, ii);
%     sumFactors = sumFactors + matrixTest(ii, ii + 1);
%     sumFactors = sumFactors + matrixTest(ii + 1, ii);
end
sumFactors
% signalInhFactors(1) = sum(spdiags(gaussianFactors, 1));
% signalInhFactors(2)
% signalInhFactors(2)
% signalInhFactors(1) = sum(spdiags(gaussianFactors, 1));
% sum(spdiags(matrixTest, 1)) + sum(spdiags(matrixTest, -1))
sum(spdiags(gaussianFactors, 2)) + sum(spdiags(gaussianFactors, -2))

%%
% delete(wbarTh);
% delete(wbarg1)

xiCalc = xi_(:, 1);
% figure(3)
% clf
% plot(1:ngax, xiCalc)

tic
xiSynth = repmat(xiCalc, 1, ngax);
for ii = 2:fix(ngax/2)
    xiSynth(:, ii) = [-flip(xiSynth(1:ii - 1, 1)); ...
        xiSynth(1:end + 1 - ii, 1)];
end
toc

%%

iCol = 200;
clf
plot(1:ngax, xi_(:, iCol), 1:ngax, xi_2(:, iCol))
% yyaxis right
% plot(1:ngax, xi_(iRow, :) - xiSynth(iRow, :)) 

%%
gaussianFactor = gInh1'*gInh2;
% size(gaussianFactor)

signalInhFactor = zeros(1, ngax);
for ig = 1:ngax
%     if xi_2(ig, 1) ~= -xi_2(ig, 2*ig)
%         disp('Error')
%     end
    signalInhFactor(ig) = sum(spdiags(gaussianFactor, -ig + 1)) + ...
        sum(spdiags(gaussianFactor, ig));
%     signalInhFactor2(ig) = sum(spdiags(gaussianFactor, -ig + 1))*2;
end

signalInh11 = sum(signalInh);

ig1 = 75;
ig2 = 1;
clf
plot(xx, signalInh_2Norm(:, ig1, ig2))
hold on
% plot(xx, crystalEseem(xi_2(ig1, ig2), xx))

%%


%%

% 
% clf
% plot(1:numel(xi{ith, iph}(:, 100)), xi{ith, iph}(:, 100))
% % 
% func = @(xx) atan( (dd_ + 2*JJ)*1e6 * ...
%                     planck/(2*pi*(bmagn*0.35)) ./ (xx - gax(100)) );
% hold on
% plot(1:numel(gax), func(gax), '--')

%%

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

%% 

%{
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
        xi_ = xi{ith, iph};
        [ngax, ~] = size(xi_);
        signalSummed = zeros(numel(xx), ngax);
        for ig1 = 1:ngax
            for ig2 = 1:ngax
                signalSummed(:, ig1) = signalSummed(:, ig1) + ...
                    crystalEseem(xi_(ig1, ig2), xx);
            end
            signal_{ith, iph} = sum(signalSummed, 2);
        end
        % Solid angle normalization
%         theta_ = thetas(ith)*pi/180;
%         signal_{ith, iph} = sum(signalSingleSumg2, 2)*sin(theta_);
%         signalPowder = signalPowder + signal_{ith, iph};
    end
    theta_ = thetas(ith)*pi/180;
    signalPowder = signalPowder + signal_{ith, iph}*sin(theta_);
end

delete(wbarTh)
% delete(wbarPh)
%}

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

clf
aa = creategaxis(g1(140, 1), g2(140, 1), glw1, glw2, dgax);
plot(aa, gaussian(aa, g1(140, 1), glw1), ...
    aa, gaussian(aa, g2(140, 1), glw2))
yyaxis right
plot(aa, gaussian(aa, g1(140, 1), glw1) .* gaussian(aa, g2(140, 1), glw2))

%%

function xx = creategaxis(g1, g2, glw1, glw2, dg)
    xxFactor = 1;
    xmin = max([g1 - xxFactor*glw1, g2 - xxFactor*glw2]);
    xmax = min([g1 + xxFactor*glw1, g2 + xxFactor*glw2]);
    xx = xmin:dg:xmax;
end















