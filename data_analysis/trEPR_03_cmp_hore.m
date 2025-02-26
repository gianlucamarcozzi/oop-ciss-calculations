% clearvars
addpath("util/")

%% CHI = 0, SINGLET INIT STATE

% Initial density matrix: singlet
[Sys, Exp] = importparam("trEPR_03_01_param.txt");
Exp.gridType = "zech";
Exp.gInhBroadening = 1;

tic
[test, th, ph] = mytrepr(Sys, Exp);
toc
aa = sum(test, [2, 3]);
aa = aa/max(abs(aa));
aa = aa;

hhFolder = "../data/digitized/";
hh = load(hhFolder + "hore_psiXband_01.mat");

clf
tL = tiledlayout(1, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile(1)
% hold on
plot(Exp.x, aa, "DisplayName", "mytrepr")
hold on
% plot(Exp.x, aa, "DisplayName", "mytrepr")
plot(hh.x*9.4945/9.6, hh.y, "DisplayName", "Hore NO CISS")
% plot(Exp.x, a2, "DisplayName", "mytrepr2")
% plot(Exp.x, a3, "DisplayName", "mytrepr3")
legend()

%%

% Initial density matrix: ud
Sys.rho = zeros(4, 4);
Sys.rho(3, 3) = 1;
Exp.gInhBroadening = 1;

% Initial density matrix: singlet
[Sys, Exp] = importparam("trEPR_03_01_param.txt");
Exp.gridType = "zech";
Exp.gInhBroadening = 0;

tic
[test, th, ph] = mytrepr(Sys, Exp);
toc
aa = sum(test, [2, 3]);
aa = aa/max(abs(aa));

hhFolder = "../data/digitized/";
hh = load(hhFolder + "hore_psiXband_05.mat");

figure(2)
clf
tL = tiledlayout(1, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile(1)
% hold on
plot(Exp.x, aa, "DisplayName", "mytrepr")
hold on
plot(hh.x*9.4945/9.6, hh.y, "DisplayName", "Hore CISS")
legend()
