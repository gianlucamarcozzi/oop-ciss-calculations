function [signal, thetas, phis] = mytrepr(Sys, Exp)

    dip = Sys.dip;
    if ~isfield(Sys, 'nNucs')
        Sys.nNucs = 0;
    end
    if ~isfield(Exp, 'gInhBroadening')
        Exp.gInhBroadening = 0;
    end

    nThetas = Exp.nThetas;
    nPhis = Exp.nPhis;
    Exp.nSolidAngles = nThetas*nPhis;

    % Grid
    [thetas, phis] = createthetaphigrid(nThetas, nPhis, Exp.gridType);

    % Direction of B0
    rVers = ang2vec(ones(nThetas, 1)*phis, thetas*ones(1, nPhis));
    % Effective g-values
    g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
    g1s = sqrt( sum( (g1Tensor*rVers).^2, 1));
    g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
    g2s = sqrt( sum( (g2Tensor*rVers).^2, 1));
    % Dipolar interaction
    zD = erot(Sys.eeFrame)*[0, 0, 1]';
    dipolars = dipinteraction(dip, rVers, zD);

    % Hyperfine: A plus (Apl) array (Apls) and A minus (Amin) array (Amins) 
    if Sys.nNucs > 0
        [Apls, Amins] = calculateApllusAmininusnohfi(...
            Sys.A, Sys.AFrAmine, rVers);
    else
        Apls = zeros(Exp.nSolidAngles, 1);
        Amins = zeros(Exp.nSolidAngles, 1);
    end

    Exp.rVers = rVers;
    Sys.g1s = g1s;
    Sys.g2s = g2s;
    Sys.dipolars = dipolars;
    Sys.zD = zD;
    Sys.Apls = Apls;
    Sys.Amins = Amins;
    
    signal = signaltreprginhbroadening(Sys, Exp);
end

%% SIGNALTREPR
function signal = signaltrepr(Sys, Exp)
    %% Variables
    g1s = Sys.g1s;
    g2s = Sys.g2s;
    nNucs = Sys.nNucs;
    nHfiLines = nNucs + 1;
    Apls = Sys.Apls;
    Amins = Sys.Amins;
    dipolars = Sys.dipolars;
    JJ = Sys.J;
    lwpp = Sys.lwpp;
    resonLwpp = mt2mhz(lwpp)*1e-3;
    rhoInit = Sys.rho;
    
    fieldAxis = Exp.x;
    nFields = Exp.nPoints;
    nSolidAngles = Exp.nSolidAngles;
    rVers = Exp.rVers;
    mwFreq = Exp.mwFreq;

    quantNumNucs = -nNucs/2:1:nNucs/2;

    %% Determine Pascal factor of hfi lines
    if nNucs > 0
        pascalMatrix = pascal(nHfiLines);
        % Antidiag
        pascs = pascalMatrix(nHfiLines:nHfiLines - 1:end - 1);
        pascs = pascs'/sum(pascs);
    else
        pascs = 1;
    end

    %% Spectrum
    % The general convention of dimensions for matrices: the indeces are
    % (iSolidAngle, iField, itrans)
    signal = zeros(nSolidAngles, nFields, 4);
    for ii = 1:nSolidAngles
        % disp(floor(ii/nSolidAngles*100))
        gPlus = 1/2*(g1s(ii) + g2s(ii));
        gMinus = 1/2*(g1s(ii) - g2s(ii));

        w0s = gvalue2freq(fieldAxis, gPlus);
        deltaws = gvalue2freq(fieldAxis, gMinus);
        
        % Variables that are scalars at the given iSolidAngle
        Apl = Apls(ii);
        Amin = Amins(ii);
        dipolar = dipolars(ii);
        for ihfi = 1:nHfiLines
            % Variables that are scalars at the given ihfi
            quantNumNuc = quantNumNucs(ihfi);
            pasc = pascs(ihfi);

            % In the for loop, the variables that change are called with
            % the convention "om" instead of "w"
            om0s = w0s + quantNumNuc*Apl;
            deltaoms = deltaws + quantNumNuc*Amin;  % (1, nField)
            omResons = myeigenenergies(om0s, deltaoms, JJ, dipolar); % (1, nField, 4)
            intensityResons = intensityresonance(rhoInit, rVers(:, ii), ...
                    deltaoms, JJ, dipolar);
            for itrans = 1:4
                %% Calculate new resonance transition
                newSig = gaussianresonancebsweep( ...
                    omResons(:, :, itrans)*1e-3, mwFreq, ...
                    resonLwpp, "lwpp");  % (1, nField)
                % Multiply by scaling factors
                newSig = intensityResons(:, :, itrans).*newSig*pasc;
                % Sum to previous signal
                signal(ii, :, itrans) = signal(ii, :, itrans) + newSig;    
            end
        end 
    end
end

%% SIGNALTREPRGINHBROADENING
function signal = signaltreprginhbroadening(Sys, Exp)
    
    gInhBroadening = Exp.gInhBroadening;
    if gInhBroadening
        N_GINH = 100;
        lw1 = mt2mhz(Sys.lw1);
        lw2 = mt2mhz(Sys.lw2);
        glw1 = freq2gvalue(mean(Exp.x), lw1);
        glw2 = freq2gvalue(mean(Exp.x), lw2);
        resonLwpp = mt2mhz(Sys.trlwpp)*1e-3;
        rng(0, "twister")  % Set random seed to 0 for reproducibility
    else
        resonLwpp = mt2mhz(Sys.lwpp)*1e-3;
        N_GINH = 1;
        glw1 = 0;
        glw2 = 0;
    end
    g1s = Sys.g1s;
    g2s = Sys.g2s;
    nNucs = Sys.nNucs;
    nHfiLines = nNucs + 1;
    Apls = Sys.Apls;
    Amins = Sys.Amins;
    dipolars = Sys.dipolars;
    zD = Sys.zD;
    JJ = Sys.J;
    rhoInit = Sys.rho;
    
    fieldAxis = Exp.x;
    nFields = Exp.nPoints;
    nSolidAngles = Exp.nSolidAngles;
    rVers = Exp.rVers;
    mwFreq = Exp.mwFreq;

    quantNumNucs = -nNucs/2:1:nNucs/2;
    %% Determine Pascal factor of hfi lines
    if nNucs > 0
        pascalMatrix = pascal(nHfiLines);
        % Antidiag
        pascs = pascalMatrix(nHfiLines:nHfiLines - 1:end - 1);
        pascs = pascs'/sum(pascs);
    else
        pascs = 1;
    end

    %% Spectrum
    signal = zeros(nSolidAngles, nFields, 4);
    ts = datetime("now", 'Format', 'HH:mm:ss');
    fprintf("%s: Start parfor: ", ts)
    for ii = 1:nSolidAngles
        % The general convention of dimensions for matrices: the indeces are
        % (igInh, iField, itrans) (except for signal)
        if gInhBroadening
            % Get numbers from gaussian distribution centered at g1 and g2
            g1rand = normrnd(g1s(ii), glw1, [N_GINH, 1]);
            g2rand = normrnd(g2s(ii), glw2, [N_GINH, 1]);
            gpls = 1/2*(g1rand + g2rand);   % (N_GINH, 1)
            gmins = 1/2*(g1rand - g2rand);  % (N_GINH, 1)
        else
            gpls = 1/2*(g1s(ii) + g2s(ii));
            gmins = 1/2*(g1s(ii) - g2s(ii));
        end
        
        w0s = gvalue2freq(fieldAxis, gpls);  % (N_GINH, nField) 
        deltaws = gvalue2freq(fieldAxis, gmins); % (N_GINH, nField)

        % Variables that are scalars at the given iSolidAngle
        Apl = Apls(ii);
        Amin = Amins(ii);
        dipolar = dipolars(ii);

        %%
        for ihfi = 1:nHfiLines
            % Variables that are scalars at the given ihfi
            quantNumNuc = quantNumNucs(ihfi);
            pasc = pascs(ihfi);
            
            % In the for loop, the variables that change are called with
            % the convention "om" instead of "w"
            om0s = w0s + quantNumNuc*Apl;
            deltaoms = deltaws + quantNumNuc*Amin;
            omResons = myeigenenergies(om0s, deltaoms, JJ, dipolar);
            intensityResons = intensityresonance(rhoInit, rVers(:, ii), ...
                    deltaoms, JJ, dipolar, zD);

            for itrans = 1:4
                %% Calculate new resonance transition
                newSig = gaussianresonancebsweep( ...
                    omResons(:, :, itrans)*1e-3, mwFreq, ...
                    resonLwpp, "lwpp");

                % Multiply by scaling factors
                newSig = squeeze(intensityResons(:, :, itrans)).*newSig*pasc;
                % Sum all (eventual) gInh contributions
                newSig = sum(newSig, 1);

                % Sum to previous signal
                signal(ii, :, itrans) = signal(ii, :, itrans) + newSig;    
            end
            % fprintf("%f\n", sum(signal(ii, :, :), "all"));
            % figure(2)
            % clf
            % plot(1:301, sum(signal(ii, :, :), 3))
            % hold on
            % plot(1:301, signal(ii, :, 1), '--')
            % plot(1:301, signal(ii, :, 2), '--')
            % plot(1:301, signal(ii, :, 3), '--')
            % plot(1:301, signal(ii, :, 4), '--')

        end 
    end
end
