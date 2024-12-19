%% Plots of ESEEM calculations
clearvars

% Definitions
% 1/2*sin(2(d-J)tau) neglected in front of every formula
sigxSinglet = ...
    @(beta, p) - (sin(beta).*sin(2.*p).^2 + sin(2.*beta).*cos(p).^4);

sigxUd = @(beta, p) cos(p).^2.*(...
    sin(beta).*sin(2.*p) - sin(2.*beta).*(2 + sin(2.*p))./2);
sigxDu = @(beta, p) cos(p).^2.*(...
    - sin(beta).*sin(2.*p) - sin(2.*beta).*(2 - sin(2.*p))./2);

sigxUdRot = @(beta, p) sin(beta).*(-sin(p(2)).^2.*cos(p(1)).^2.*sin(p(1)).^2 ...
    + 4.*cos(p(2)).*cos(p(1)).^3.*sin(p(1)))./2 ...
    + sin(2.*beta).*(-cos(p(2)).^2.*cos(p(1)).^2 ...
    + sin(p(2)).^2.*cos(p(1)).^2.*sin(p(1)).^2./4 ...
    + sin(p(2)).^2.*(cos(p(3)).^2 - sin(p(3)).^2).*cos(p(1)).^2 ...
    - cos(p(2)).*cos(p(1)).^3.*sin(p(1)));
sigxDuRot = @(beta, p) sin(beta).*(-sin(p(2)).^2.*cos(p(1)).^2.*sin(p(1)).^2 ...
    - 4.*cos(p(2)).*cos(p(1)).^3.*sin(p(1)))./2 ...
    + sin(2.*beta).*(-cos(p(2)).^2.*cos(p(1)).^2 ...
    + sin(p(2)).^2.*cos(p(1)).^2.*sin(p(1)).^2./4 ...
    + sin(p(2)).^2.*(cos(p(3)).^2 - sin(p(3)).^2).*cos(p(1)).^2 ...
    + cos(p(2)).*cos(p(1)).^3.*sin(p(1)));

bb = linspace(0, pi, 1000);

% Eckvahl
ddip = 52/(2.48^3);  % Mhz
ddip = -3*ddip;  % Bittl defines D(th) = ddip(cos(th)^2 - 1/3)
Jexch = 0.1;  % MHz
Jexch = -2*Jexch;  % J_bittl = -2Jeck;
gA = [2.0031, 2.0044, 2.0046];
gB = [2.0034, 2.0041, 2.0043];

% PSI values from PSI example easyspin
% Jexch = unitconvert(1e-3,'mT->MHz'); % MHz
% ddip = unitconvert(-0.170,'mT->MHz'); % MHz
% gA = [2.0033, 2.0024, 2.0020];
% gB = [2.0065, 2.0053, 2.0022];

theta = pi/2;
phi = pi/2;
xiTheta0 = calculatexi(gA, gB, ddip, Jexch, 0, 0);
xiRot = calculatexi(gA, gB, ddip, Jexch, theta, phi);

sigSingletTheta0 = sigxSinglet(bb, xiTheta0);
sigSingletRotYAxis = sigxSinglet(bb, xiRot);
sigCissTheta0 = sigxUd(bb, xiTheta0) + sigxDu(bb, xiTheta0);
sigCissRotYAxis = sigxUdRot(bb, [xiRot, theta, phi]) ...
    + sigxDuRot(bb, [xiRot, theta, phi]);

% Plot
load("plotColors.mat")
plotLw = 1.5;

ff = figure();
clf
plot(bb, sigSingletTheta0, 'LineWidth', plotLw, 'Color', plotColors(1))
hold on
plot(bb, sigCissTheta0, 'LineWidth', plotLw, 'Color', plotColors(2))
plot(bb, sigSingletRotYAxis, 'LineWidth', plotLw, 'Color', plotColors(3))
plot(bb, sigCissRotYAxis, 'LineWidth', plotLw, 'Color', plotColors(4))

xlim([min(bb), max(bb)])

ax = gca;
grid("on")
ax.XGrid = 'off';
% chH = get(gca, 'Children');
% set(gca, 'Children', flipud(chH))

labelSingletTheta0 = append("Singlet, ", char(952), "=0");
labelCissTheta0 = append("CISS, ", char(952), "=0");
labelSingletRot = append("Singlet, ", char(952), "=", char(960), ...
    "/2, ", char(966), "=", char(960), "/2");
labelCissRot = append("CISS, ", char(952), "=", char(960), ...
    "/2, ", char(966), "=", char(960), "/2");
legend({labelSingletTheta0, labelCissTheta0, labelSingletRot, labelCissRot}, ...
    "Location", "southeast")

xAxisLabel = append("Turning angle ", char(946));
yAxisLabel = append("S_x(T + 2", char(964), ")");
yticks(0)
labelaxesfig(gca, xAxisLabel, yAxisLabel)

text(.05, .9, "D, J, g_A and g_B taken from [3]", "Units", "normalized", ...
    "FontSize", 14)
ff.Position(3) = 1.3*ff.Position(3);
savefigas(ff, ...
    '/home/gianlum33/files/projects/oop_ciss_calculations/images/posterRsc_oopEseem_eckvahl.svg')

%% 
function xi = calculatexi(gA, gB, ddip, Jexch, theta, phi)
    
    % Effective g-value
    nVers = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    gAeff = sqrt(gA.^2*nVers'.^2);
    gBeff = sqrt(gB.^2*nVers'.^2);
    
    dOmega = (gAeff - gBeff)*bmagn*0.35*2*pi/planck*1e-6;
    ddipTheta = ddip*(cos(theta)^2 - 1/3);
    
    xi = atan((ddipTheta + 2*Jexch)/dOmega);
end
