function [signal, thetas] = mytrepr(Sys, Exp)
    
    dip = Sys.dip;
    JJ = Sys.J;
    lwpp = Sys.lwpp;

    rhoInit = Sys.rho;  % This should be in the T+, S, T0, T- basis
    if isfield(Sys, 'nNuc')
        nNuc = Sys.nNuc;
    else
        nNuc = 0;
    end
    nHfiLine = nNuc + 1;

    fieldAxis = Exp.x;
    nField = Exp.nPoints;
    mwFreq = Exp.mwFreq;
    nTheta = Exp.nTheta;
    nPhi = Exp.nPhi;

    % Grid
    [thetas, phis] = createthetaphigrid(nTheta, nPhi, Exp.gridType);
    nSolidAngle = nTheta*nPhi;

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
    if nNuc > 0
        [APlus, AMinus] = calculateaplusaminusnohfi(...
            Sys.A, Sys.AFrame, rVers);
    else
        APlus = zeros(nSolidAngle, 1);
        AMinus = zeros(nSolidAngle, 1);
    end

    signal = zeros(nField, nSolidAngle);
    for ii = 1:nSolidAngle
        disp(floor(ii/nSolidAngle*100))
        gPlus = 1/2*(g1(ii) + g2(ii));
        gMinus = 1/2*(g1(ii) - g2(ii));

        w0_ = gvalue2freq(fieldAxis, gPlus);
        deltaw_ = gvalue2freq(fieldAxis, gMinus);

        quantNumNuc = -nNuc/2:1:nNuc/2;
        if nNuc > 0
            pascalMatrix = pascal(nHfiLine);
            % Antidiag
            pascalFactor = pascalMatrix(nHfiLine:nHfiLine - 1:end - 1);
            pascalFactor = pascalFactor'/sum(pascalFactor);
        else
            pascalFactor = 1;
        end

        for itrans = 1:4
            for ihfi = 1:nHfiLine
                w0__ = w0_ + quantNumNuc(ihfi)*APlus(ii);
                deltaw__ = deltaw_ + quantNumNuc(ihfi)*AMinus(ii);
                wReson = myeigenenergies(w0__, deltaw__, JJ, dd(ii), itrans);
                intensityReson = intensityresonance(...
                    rhoInit, rVers(:, ii), deltaw__, JJ, dd(ii), itrans);
                signal__ = gaussianresonancebsweep( ...
                    wReson*1e-3, mwFreq, mt2mhz(lwpp)*1e-3, "lwpp");
                signal(:, ii) = signal(:, ii) + intensityReson'.*signal__'*pascalFactor(ihfi);

            end 
        end
    end

end

