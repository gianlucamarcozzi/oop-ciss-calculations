%%
clearvars
addpath("../util/")

%% COMPARISON BETWEEN TWO DIFFERENT WAYS OF CHOOSING POWDER AVG ANGLES

% General import
[Sys, Exp] = importparam("trEPR_01_01_param.txt");

% Zech (constant delta theta)
Exp.gridType = "zech";
[testz, thz] = mytrepr(Sys, Exp);
aaz = sum(testz, 2)/sum(sin(repmat(thz, [Exp.nPhi, 1])), 'all');
aaz = squeeze(aaz);
aaz = sum(aaz, 2);

% My way (constant theta)
Exp.gridType = "other";
[test, th] = mytrepr(Sys, Exp);
aaa = averageoversolidangle(test, th, Exp.nPhi, 2);
aaa = sum(aaa, 2);

% PLOT
% figure()
clf
tL = tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");
ax1 = nexttile;
hold on
box on
% plot(Exp.x, aa0/max(abs(aa0)), "DisplayName", "mytrepr")
plot(Exp.x, aaa/max(aaa), ...
    "DisplayName", "Const. " + char(0x0394) + char(0x03B8))
plot(Exp.x, aaz/max(aaz), ...
    "DisplayName", "Const. " + char(0x0394) + "cos(" + char(0x03B8) + ")")
legend("Location", "northwest")
ylim(setaxlim(aaz/max(aaz), 1.05))

ax2 = nexttile;
hold on
box on
% plot(Exp.x, aa1/max(abs(aa1)) - p2{:, 2}/max(abs(p2{:, 2})))
% plot(Exp.x, aa2/max(abs(aa2)) - p2{:, 2}/max(abs(p2{:, 2})))
plot(Exp.x, aaa/max(aaa) - aaz/max(aaz), 'k', ...
    "DisplayName", "Const. " + char(0x0394) + char(0x03B8) + " " + ...
    char(0x2014) + " const. " + char(0x0394) + "cos(" + char(0x03B8) + ")")
legend("Location", "northwest")
linkaxes([ax1, ax2], 'x');
xticklabels(ax1, {});
labelaxesfig(tL, "Magnetic field / mT", "trEPR signal / a.u.")
ylim(setaxlim(aaa/max(aaa) - aaz/max(aaz), 1.05))

%% COMPARISON WITH EASYSPIN

[Sys, Exp] = importparam("trEPR_01_02_easyspin_param.txt");
if Sys.nNuc == 0
    Sys = nucspinrmv(Sys, 1);
end

[B, ySim] = pepper(Sys, Exp);
clf
plot(B, ySim/max(abs(ySim)), 'DisplayName', 'Easyspin')
hold on
plot(Exp.x, aaa/max(aaa), 'DisplayName', 'Gianluca')
legend()

%% MYTREPR vs p2 SCRIPT

% My trepr
[Sys, Exp] = importparam("trEPR_01_03_param.txt");
Exp.gridType = "zech";

[test, th] = mytrepr(Sys, Exp);
aa = sum(test, 2)/sum(sin(repmat(th, [Exp.nPhi, 1])), 'all');
aa = squeeze(aa);
aa = sum(aa, 2);
% aa = aa(:, 3);

% Import result of p2 script
p2folder = "../from_bittl/PSI_trEPR/";
p2T = readtable(p2folder + "test2_spec_01_cmp.asc", 'FileType', 'text');
p2x = p2T{:, 1}/10;
p2y = p2T{:, 2}/max(abs(p2T{:, 2}));
hold on
% plot(p2{:, 1}/10, p2{:, 2}/max(abs(p2{:, 2})), 'DisplayName', 'p2 script')

% figure()
clf
tL = tiledlayout(2, 1, "TileSpacing", "compact", "Padding", "compact");
ax1 = nexttile;
hold on
box on
% plot(Exp.x, aa0/max(abs(aa0)), "DisplayName", "mytrepr")
% plot(Exp.x, aaa/max(aaa), ...
%     "DisplayName", "Const. " + char(0x0394) + char(0x03B8))
plot(Exp.x, aa/max(abs(aa)), ...
    "DisplayName", "Const. " + char(0x0394) + "cos(" + char(0x03B8) + ")")
plot(p2x, p2y, "DisplayName", "p2 script")
legend("Location", "northwest")
% ylim(setaxlim(aa/max(aa), 1.05))

ax2 = nexttile;
hold on
box on
% plot(Exp.x, aa1/max(abs(aa1)) - p2{:, 2}/max(abs(p2{:, 2})))
% plot(Exp.x, aa2/max(abs(aa2)) - p2{:, 2}/max(abs(p2{:, 2})))
plot(Exp.x, aa/max(abs(aa)) - p2y, 'k', ...
    "DisplayName", "Difference")
legend("Location", "northwest")
linkaxes([ax1, ax2], 'x');
xticklabels(ax1, {});
labelaxesfig(tL, "Magnetic field / mT", "trEPR signal / a.u.")
ylim(setaxlim(aa/max(aa) - p2y, 1.05))
