%% 1. Exhange coupling
clearvars

Sys.S = [1/2 1/2];
Sys.g = [2, 2.0075];

Sys.initState = 'singlet';

% Exchange coupling
nj = 100;
Js = linspace(1e-2, 100, nj);  % MHz

gPlus = (Sys.g(1) + Sys.g(2))/2;
mwFreq = 9.7;  % GHz
B0 = planck/bmagn/gPlus*mwFreq*1e12;  % mT
nPoint = 5000;
msAx = linspace(9.4, 10, nPoint);  % GHz, mw sweep axis
bsAx = mhz2mt(msAx*1e3, gPlus);  % mT, b field sweep

% 
% Frequency sweep
%
mhzLw = 10;
Sys.lw = mhzLw;  % MHz
 
Exp.mwRange = [min(msAx), max(msAx)];
Exp.Field = B0;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

msSim = zeros(nj, nPoint);

for ij = 1:nj
    Sys.J = Js(ij);
    msSim(ij, :) = pepper(Sys, Exp);
    % msSim(ij, :) = msSim(ij, :)/max(abs(msSim(ij, :)));
end

msTotIntegral = sum(msSim, 2)/nPoint;

%
% Field sweep
%
clear("Exp")

Sys.lw = mhz2mt(mhzLw);  % mT

Exp.Range = [min(bsAx), max(bsAx)];
Exp.mwFreq = mwFreq;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

bsSim = zeros(nj, nPoint);

for ij = 1:nj
    Sys.J = Js(ij);
    bsSim(ij, :) = pepper(Sys, Exp);
    % bsSim(ij, :) = bsSim(ij, :)/max(abs(bsSim(ij, :)));
end

bsTotIntegral = sum(bsSim, 2)/nPoint;

clf
plot(Js, msTotIntegral, 'o-', 'DisplayName', 'Mw sweep')
yyaxis right
plot(Js, bsTotIntegral, 'o-', 'DisplayName', 'B0 sweep')
xlabel('J / MHz')
yyaxis left
ylabel('Integral of trEPR spectrum')
legend()

%% Plot

ij = 100;
clf
plot(msAx, msSim(ij, :))
% yyaxis right
plot(msAx, flip(bsSim(ij, :)))
yyaxis right
plot(msAx, msSim(ij, :) - flip(bsSim(ij, :)))


%% Dipolar coupling
clearvars

Sys.S = [1/2 1/2];
Sys.g = [2, 2.0075];

Sys.initState = 'singlet';

% Dipolar coupling
ndip = 100;
dips = linspace(1, 100, ndip);  % MHz

gPlus = (Sys.g(1) + Sys.g(2))/2;
mwFreq = 9.7;  % GHz
B0 = planck/bmagn/gPlus*mwFreq*1e12;  % mT
nPoint = 5000;
msAx = linspace(9.4, 10, nPoint);  % GHz, mw sweep axis
bsAx = mhz2mt(msAx*1e3, gPlus);  % mT, b field sweep

% 
% Frequency sweep
%
mhzLw = 10;
Sys.lw = mhzLw;  % MHz
 
Exp.mwRange = [min(msAx), max(msAx)];
Exp.Field = B0;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

msSim = zeros(ndip, nPoint);

for ij = 1:ndip
    Sys.dip = dips(ij);
    msSim(ij, :) = pepper(Sys, Exp);
    % msSim(ij, :) = msSim(ij, :)/max(abs(msSim(ij, :)));
end

msTotIntegral = sum(msSim, 2)/nPoint;

%
% Field sweep
%
clear("Exp")

Sys.lw = mhz2mt(mhzLw);  % mT

Exp.Range = [min(bsAx), max(bsAx)];
Exp.mwFreq = mwFreq;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

bsSim = zeros(ndip, nPoint);

for ij = 1:ndip
    Sys.dip = dips(ij);
    bsSim(ij, :) = pepper(Sys, Exp);
    % bsSim(ij, :) = bsSim(ij, :)/max(abs(bsSim(ij, :)));
end

bsTotIntegral = sum(bsSim, 2)/nPoint;

clf
plot(dips, msTotIntegral, 'o-', 'DisplayName', 'Mw sweep')
yyaxis right
plot(dips, bsTotIntegral, 'o-', 'DisplayName', 'B0 sweep')
xlabel('J / MHz')
yyaxis left
ylabel('Integral of trEPR spectrum')
legend()

%% Plot

ij = 4;
clf
plot(msAx, msSim(ij, :))
yyaxis right
plot(msAx, flip(bsSim(ij, :)))

%% Crystal rotation
clearvars

Sys.S = [1/2 1/2];
Sys.g = [2, 2.0075];
Sys.dip = 5;  % MHz

Sys.initState = 'singlet';

gPlus = (Sys.g(1) + Sys.g(2))/2;
mwFreq = 9.7;  % GHz
B0 = planck/bmagn/gPlus*mwFreq*1e12;  % mT
nPoint = 5000;
msAx = linspace(9.4, 10, nPoint);  % GHz, mw sweep axis
bsAx = mhz2mt(msAx*1e3, gPlus);  % mT, b field sweep

% 
% Frequency sweep
%
mhzLw = 10;
Sys.lw = mhzLw;  % MHz
 
Exp.mwRange = [min(msAx), max(msAx)];
Exp.Field = B0;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

Opt.GridSize = [91 0];  % 181 orientations
Opt.GridSymmetry = 'C1';

% Powder average
msSim = pepper(Sys, Exp, Opt);

clear("Opt")

Exp.CrystalSymmetry = 1;
Exp.MolFrame = [0 0 0]*pi/180;
Exp.SampleFrame = [0 0 0]*pi/180; 

Opt.separate = 'orientations';
rho = deg2rad(linspace(0, 180, 181));              % 10 degree increments
solidAngleWeight = sin(rho')/sum(sin(rho), 'all');

% Rotation around x
Exp.SampleRotation = {'x', rho};
msSimRotX = pepper(Sys, Exp, Opt);
msPowderRotX  = sum(solidAngleWeight.*msSimRotX);

% Rotation around y
Exp.SampleRotation = {'y', rho};
msSimRotY = pepper(Sys, Exp, Opt);
msPowderRotY  = sum(solidAngleWeight.*msSimRotY);

%
% Field sweep
%
clear("Exp")

Sys.lw = mhz2mt(mhzLw);  % mT

Exp.Range = [min(bsAx), max(bsAx)];
Exp.mwFreq = mwFreq;
Exp.nPoints = nPoint;
Exp.Harmonic = 0;

clear("Opt")
Opt.GridSize = [91 0];  % 181 orientations
Opt.GridSymmetry = 'C1';

% Powder average
bsSim = pepper(Sys, Exp, Opt);

clear("Opt")

Exp.CrystalSymmetry = 1;
Exp.MolFrame = [0 0 0]*pi/180;
Exp.SampleFrame = [0 0 0]*pi/180; 

Opt.separate = 'orientations';
rho = deg2rad(linspace(0, 180, 181));
solidAngleWeight = sin(rho')/sum(sin(rho), 'all');

% Rotation around x
Exp.SampleRotation = {'x', rho};
bsSimRotX = pepper(Sys, Exp, Opt);
bsPowderRotX  = sum(solidAngleWeight.*bsSimRotX);

% Rotation around y
Exp.SampleRotation = {'y', rho};
bsSimRotY = pepper(Sys, Exp, Opt);
bsPowderRotY  = sum(solidAngleWeight.*bsSimRotY);

%% Plot

clf
plot(msAx, msPowderRotX, msAx, msPowderRotY)
% clf
% plot(bsAx, bsPowderRotX, bsAx, bsPowderRotY)

