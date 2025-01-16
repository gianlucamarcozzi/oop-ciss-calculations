% clearvars
addpath("util/")

%%

[Sys, Exp] = importparam("trEPRstickSpectra_02_01_param.txt");
Exp.gridType = "zech";

[test, th] = mytrepr(Sys, Exp);
aa = sum(test, 2);
aa = aa/max(aa);
aa = real(aa);

%%

hhFolder = "/home/gianluca/files/work/projects/oop-ciss-calculations/data/digitized/";
hh = load(hhFolder + "hore_psiXband_05.mat");

clf
tL = tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile(1)
% hold on
plot(Exp.x, aa, "DisplayName", "Hore")
legend()

nexttile
plot(hh.x, hh.y, "DisplayName", "Hore")
xlim(setaxlim(hh.x, 2))
legend()
