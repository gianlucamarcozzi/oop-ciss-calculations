%%
clearvars
addpath("../util/")

%% CHI = 0, SINGLET INIT STATE

% Initial density matrix: singlet
[Sys, Exp] = importparam("trEPR_03_01_param.txt");
Exp.gridType = "zech";
Exp.gInhBroadening = 0;

%
tic
aout = mytrepr(Sys, Exp);
toc
% Average over solid angle (no sin(theta) weight needed if gridType zech)
a1 = sum(aout, 1);
aP = sum(a1(:, :, [1, 2]) , 3);
aQ = sum(a1(:, :, [3, 4]) , 3);
aa = sum(a1 , 3);
aNorm = max(abs(aa));
aP = aP/aNorm;
aQ = aQ/aNorm;
aa = aa/aNorm;
aP = real(aP);
aQ = real(aQ);
aa = real(aa);

hhFolder = "../data/digitized/";
hh = load(hhFolder + "hore_psiXband_01.mat");

figure(1)
clf
tL = tiledlayout(1, 1, "TileSpacing", "compact", "Padding", "compact");
nexttile(1)
% hold on
plot(Exp.x, aa, "DisplayName", "mytrepr")
hold on
% plot(Exp.x, aa, "DisplayName", "mytrepr")
plot(hh.x*9.4945/9.6, hh.y, "DisplayName", "Hore NO CISS")
% plot(Exp.x, aP, "DisplayName", "aP")
% plot(Exp.x, aQ, "DisplayName", "aQ")
legend()

iimax = 0;
for ii = round(linspace(1, 200, 200))
    iimax = iimax + 1;
    [~, imax1(iimax)] = max(abs(test(ii, :, 4)));
end

%%

% Initial density matrix: singlet
[Sys, Exp] = importparam("trEPR_03_01_param.txt");
Xs = [0, 0.05, 0.1, 0.15];
Xs = linspace(0, 1, 11);
nSim = numel(Xs);
Exp.gridType = "zech";
Exp.gInhBroadening = 0;
for ii = 1:nSim
    Sys.chi = acos(1 - Xs(ii))*180/pi;
    Sys = getrhofromparamfile(Sys);

    tic
    aouts{ii} = mytrepr(Sys, Exp);
    toc
end
%%
figure(2)
clf
tiledlayout("flow", "TileSpacing", "compact", "Padding", "compact")

for ii = 1:nSim
    % Average over solid angle (no sin(theta) weight needed if gridType zech)
    a1 = sum(aouts{ii}, 1);
    aP = sum(a1(:, :, [1, 2]) , 3);
    aQ = sum(a1(:, :, [3, 4]) , 3);
    aa = sum(a1 , 3);
    % aNorm = max(abs(aa));
    aNorm = 1000;
    aP = aP/aNorm;
    aQ = aQ/aNorm;
    aa = aa/aNorm;
    aP = real(aP);
    aQ = real(aQ);
    aa = real(aa);

    nexttile
    plot(Exp.x, aa, '-', "DisplayName", "tot")
    hold on

    plot(Exp.x, aP, '--', "DisplayName", "P")
    hold on
    plot(Exp.x, aQ, '--', "DisplayName", "Q")
    ylim([-1, 1])
    title(string(Xs(ii)))
end

%%
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
plot(Exp.x, aP, '--', "DisplayName", "P")
plot(Exp.x, aQ, '--', "DisplayName", "Q")

legend()

iimax = 0;
for ii = round(linspace(1, 200, 200))
    iimax = iimax + 1;
    [~, imax2(iimax)] = max(abs(test(ii, :, 4)));
end

