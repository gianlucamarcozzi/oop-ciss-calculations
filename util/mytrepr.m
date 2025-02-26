function [signal, thetas, phis] = mytrepr(Sys, Exp)

    dip = Sys.dip;
    lwpp = Sys.lwpp;

    if ~isfield(Sys, 'nNuc')
        Sys.nNuc = 0;
    end

    nTheta = Exp.nTheta;
    nPhi = Exp.nPhi;
    Exp.nSolidAngle = nTheta*nPhi;

    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi, Exp.gridType);

    % Direction of B0
    rVers = ang2vec(ones(nTheta, 1)*phis, thetas*ones(1, nPhi));
    % Effective g-values
    g1Tensor = rotatematrixeuangles(diag(Sys.g(1, :)), Sys.gFrame(1, :));
    g1 = sqrt( sum( (g1Tensor*rVers).^2, 1));
    g2Tensor = rotatematrixeuangles(diag(Sys.g(2, :)), Sys.gFrame(2, :));    
    g2 = sqrt( sum( (g2Tensor*rVers).^2, 1));
    % Dipolar interaction
    zD = erot(Sys.eeFrame)*[0, 0, 1]';
    dd = dipinteraction(dip, rVers, zD);

    % Hyperfine
    if Sys.nNuc > 0
        [APlus, AMinus] = calculateaplusaminusnohfi(...
            Sys.A, Sys.AFrame, rVers);
    else
        APlus = zeros(Exp.nSolidAngle, 1);
        AMinus = zeros(Exp.nSolidAngle, 1);
    end

    Exp.rVers = rVers;
    Sys.g1 = g1;
    Sys.g2 = g2;
    Sys.dd = dd;
    Sys.APlus = APlus;
    Sys.AMinus = AMinus;
    
    if isfield(Exp, "gInhBroadening") && Exp.gInhBroadening
        signal = signaltreprginhbroadening(Sys, Exp);
    else
        signal = signaltrepr(Sys, Exp);
    end
end

%% SIGNALTREPR
function signal = signaltrepr(Sys, Exp)
    
    g1 = Sys.g1;
    g2 = Sys.g2;
    nNuc = Sys.nNuc;
    nHfiLine = nNuc + 1;
    APlus = Sys.APlus;
    AMinus = Sys.AMinus;
    dd = Sys.dd;
    JJ = Sys.J;
    lwpp = Sys.lwpp;
    rhoInit = Sys.rho;
    
    fieldAxis = Exp.x;
    nField = Exp.nPoints;
    nSolidAngle = Exp.nSolidAngle;
    rVers = Exp.rVers;
    mwFreq = Exp.mwFreq;

    quantNumNuc = -nNuc/2:1:nNuc/2;
    % Determine Pascal factor of hfi lines
    if nNuc > 0
        pascalMatrix = pascal(nHfiLine);
        % Antidiag
        pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
        pascalFactor = pascalFactor'/sum(pascalFactor);
    else
        pascalFactor = 1;
    end
    
    signal = zeros(nField, nSolidAngle, 4);
    for ii = 1:nSolidAngle
        % disp(floor(ii/nSolidAngle*100))
        gPlus = 1/2*(g1(ii) + g2(ii));
        gMinus = 1/2*(g1(ii) - g2(ii));

        w0 = gvalue2freq(fieldAxis, gPlus);
        deltaw = gvalue2freq(fieldAxis, gMinus);

        for ihfi = 1:nHfiLine
            % In the for loop, the variables that change are called with
            % the convention "om" instead of "w"
            om0 = w0 + quantNumNuc(ihfi)*APlus(ii);
            deltaom = deltaw + quantNumNuc(ihfi)*AMinus(ii);
            omReson = myeigenenergies(om0, deltaom, JJ, dd(ii));
            intensityReson = intensityresonance(rhoInit, rVers(:, ii), ...
                    deltaom, JJ, dd(ii));
            for itrans = 1:4
                % Calculate new resonance transition
                newSig = gaussianresonancebsweep( ...
                    omReson(itrans, :)*1e-3, mwFreq, ...
                    mt2mhz(lwpp)*1e-3, "lwpp");
                % Multiply by scaling factors
                newSig = intensityReson(itrans, :)'.*newSig'*pascalFactor(ihfi);
                % Sum to previous signal
                signal(:, ii, itrans) = signal(:, ii, itrans) + newSig;
                % signal(:, ii) = signal(:, ii) + ...
                %     intensityReson'.*signal__'*pascalFactor(ihfi);
    
            end
        end 
    end
end

%% SIGNALTREPRGINHBROADENING
function signal = signaltreprginhbroadening(Sys, Exp)
    
    N_RAND_POINTS = 2000;
    g1 = Sys.g1;
    g2 = Sys.g2;
    nNuc = Sys.nNuc;
    nHfiLine = nNuc + 1;
    APlus = Sys.APlus;
    AMinus = Sys.AMinus;
    dd = Sys.dd;
    JJ = Sys.J;
    trlwpp = Sys.trlwpp;
    rhoInit = Sys.rho;
    lw1 = mt2mhz(Sys.lw1);
    lw2 = mt2mhz(Sys.lw2);
    glw1 = freq2gvalue(mean(Exp.x), lw1);
    glw2 = freq2gvalue(mean(Exp.x), lw2);

    fieldAxis = Exp.x;
    nField = Exp.nPoints;
    nSolidAngle = Exp.nSolidAngle;
    rVers = Exp.rVers;
    mwFreq = Exp.mwFreq;

    quantNumNuc = -nNuc/2:1:nNuc/2;
    % Determine Pascal factor of hfi lines
    if nNuc > 0
        pascalMatrix = pascal(nHfiLine);
        % Antidiag
        pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
        pascalFactor = pascalFactor'/sum(pascalFactor);
    else
        pascalFactor = 1;
    end
    
    % Set random seed to 0 for reproducibility
    rng(0, "twister")

    signal = zeros(nField, nSolidAngle, 4);
    ts = datetime("now", 'Format', 'HH:mm:ss');
    fprintf("%s: Start parfor: ", ts)
    parfor ii = 1:nSolidAngle
        % fprintf("%i, ", floor(ii/nSolidAngle*100))
        % Generate numbers in a gaussian distribution centered at g1 and g2
        g1rand = normrnd(g1(ii), glw1, [N_RAND_POINTS, 1]);
        g2rand = normrnd(g2(ii), glw2, [N_RAND_POINTS, 1]);
        gPlus = 1/2*(g1rand + g2rand);
        gMinus = 1/2*(g1rand - g2rand);
        % gPlus = gPlus(:);
        % gMinus = gMinus(:);

        w0 = gvalue2freq(fieldAxis, gPlus);
        deltaw = gvalue2freq(fieldAxis, gMinus);

        for ihfi = 1:nHfiLine
            % In the for loop, the variables that change are called with
            % the convention "om" instead of "w"
            om0 = w0 + quantNumNuc(ihfi)*APlus(ii);
            deltaom = deltaw + quantNumNuc(ihfi)*AMinus(ii);
            om0Resh = reshape(om0, [1, size(om0)]);
            deltaomResh = reshape(deltaom, [1, size(deltaom)]);
            omReson = myeigenenergies(om0Resh, deltaomResh, JJ, dd(ii));
            intensityReson = intensityresonance(rhoInit, rVers(:, ii), ...
                    deltaom, JJ, dd(ii));
            for itrans = 1:4
                % Calculate new resonance transition
                newSig = gaussianresonancebsweep( ...
                    omReson(itrans, :, :)*1e-3, mwFreq, ...
                    mt2mhz(trlwpp)*1e-3, "lwpp");
                newSig = squeeze(newSig);

                % Multiply by scaling factors
                newSig = squeeze(intensityReson(itrans, :, :))'.*newSig'*pascalFactor(ihfi);
                newSig = sum(newSig, 2);
                % Sum to previous signal
                signal(:, ii, itrans) = signal(:, ii, itrans) + newSig;
                % signal(:, ii) = signal(:, ii) + ...
                %     intensityReson'.*signal__'*pascalFactor(ihfi);
    
            end
        end 
    end
end
