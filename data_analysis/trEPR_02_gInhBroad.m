%%
clearvars
addpath("../util/")

%% NORMAL TREPR VS GINHBROAD

% Normal trEPR
[Sys, Exp] = importparam("trEPR_02_01_param.txt");
Exp.gridType = "zech";
Exp.gInhBroadening = 0;

[test, th, ph] = mytrepr(Sys, Exp);
aa = sum(test, [2, 3]);  % Sum over orientations
aa = aa/max(abs(aa));

% With gInhBroadening
Exp.gInhBroadening = 1;

tic
[test, th, ph] = mytrepr(Sys, Exp);
toc
bb = sum(test, [2, 3]);  % Sum over orientations
bb = bb/max(abs(bb));

% Import result of p2 script
p2folder = "../from_bittl/PSI_trEPR/";
p2T = readtable(p2folder + "test2_spec_01_cmp.asc", 'FileType', 'text');
p2x = p2T{:, 1}/10;
p2y = p2T{:, 2}/max(abs(p2T{:, 2}));
hold on

figure(1)
clf
nexttile
plot(Exp.x, aa)
hold on
plot(Exp.x, bb)
plot(p2x, p2y)
nexttile
plot(Exp.x, aa - bb)

