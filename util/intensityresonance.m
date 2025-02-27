function intensityResonance = intensityresonance(rhoInit, rVers, deltaw, J, d, zD)
    arguments
        rhoInit     % (4, 4) double
        rVers       % (3, 1) double
        deltaw      % (1, nField) or double
        J           % (1, 1) double
        d          % (1, 1) double
        zD
    end
    %% Mixing angle and transition probabilities in eigenbasis of spinH
    alpha = 1/2*atan(deltaw./(J + d/2));  % (1, nField)
    szalpha = size(alpha);
    alpha = reshape(alpha, [1, numel(alpha)]);
    
    transProb = [sin(alpha).^2/2; ...    % 1 to 2
                cos(alpha).^2/2; ...     % 3 to 4
                cos(alpha).^2/2; ...     % 1 to 3
                sin(alpha).^2/2];        % 2 to 4

    %% Diagonalize rho to get population differences
    % rhoInit is assumed to be in the product basis
    % Rotate according to angles theta and phi
    rho0 = changeframe0toB(rhoInit, rVers, zD);

    % Diagonalization matrix from product basis to eigenbasis H
    UMatrix = diagonalizingmatrix(alpha, "product");  % (4, 4, nField)
    rho1 = pagemtimes(UMatrix, ...
        pagemtimes(rho0, pagetranspose(UMatrix)));  % (4, 4, nField)
    
    % Transform from product to coupled NO NEED
    % rho = rotateproduct2coupled(rho1);
    rho = rho1;

    % Diagonal elements
    % theta = acos(rVers(3)/norm(rVers));
    populs = mydiag3d(rho);  % (4, nField)
    % populs(2, :) = 1/2*(sin(alpha(1))*cos(theta) + cos(alpha(1)))^2;
    % populs(3, :) = 1/2*(cos(alpha(1))*cos(theta) - sin(alpha(1)))^2;
    populDiff = [populs(1, :) - populs(2, :); ...
                populs(3, :) - populs(4, :); ...
                populs(1, :) - populs(3, :); ...
                populs(2, :) - populs(4, :)];
    
    % p = 1/4*sin(theta)^2 - populs(1);
    % diff2 = - populs(2); 
    % diff3 = - populs(3); 
    % fprintf("%f, %f, %f\n", diff1, diff2, diff3)

    %% Intensity of resonances
    intensityResonance = -populDiff.*transProb;  % (4, nField)
    intensityResonance = reshape(intensityResonance, [4, szalpha]);  % (4, 1, nField)
    intensityResonance = permute(intensityResonance, [2, 3, 1]);  % (1, nField, 4)
end

function rho = changeframe0toB(rhoInit, rVers, zD)
    %% Go from frame 0 (x, y, z) to frame D (x', y', zD). 
    % The transformation matrix is M_0toD = erot(phiD, thetaD, 0) or 
    % erot(phiD, thetaD, -phi) if we want to mantain one of the axes fixed.
    % To check that this is the right transformation:
    % M_0toD * zD_0 should give zD_D = [0; 0; 1], where zD_0 are the
    % coordinates of zD in frame 0, equivalent to the parameter zD passed
    % to this function, and zD_D are the coordinates of zD in frame D. 
    % Note also that, if rhoInit is up-down for example, the representation
    % rho = 0 everywhere except for rho(2, 2) = 1 is only correct in
    % the (x', y', zD) and not in the initial frame (unless zD is along the
    % initial z).
    
    % Polar angles of zD
    [phiD, thetaD] = vec2ang(zD);
    % Rotation matrix
    M_0toD = erot(phiD, thetaD, 0);
    % Direction of B0 in frame D
    rVersD = M_0toD*rVers;

    %% Change of frame from D (x', y', zD) to B (x'', y'', B0Vers)
    % Here the procedure is analogous. rhoInit is the density matrix in the
    % D frame. To express it in the B frame, we calculate the appropriate
    % angles and then apply a rotation. rhoInit is a density matrix
    % therefore the rotation will be represented as expm[-i(Sza+Szb)phi]
    % and expm[-i(Sya+Syb)theta] with the appropriate phi and theta angles.
    % Note: the convention R_y = expm[-i(Sya+Syb)theta] (with the minus) is
    % such that rho_new = R_y(theta) * rho * R_y(theta)^(-1)

    % Polar angles of rVersD (hence in the D frame)
    [phiB, thetaB] = vec2ang(rVersD);
    
    % fprintf("%f, %f\n", theta*180/pi, phi)
    % Sza + Szb
    szapszb = 1/2*diag([1, 0, 0, -1]);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sya + Syb
    % syapsyb = 
    %   0   -i/2    -i/2    0
    %   i/2   0       0     -i/2
    %   i/2   0       0     -i/2
    %   0   i/2     i/2     0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    syapsyb = zeros(4, 4);
    syapsyb(1, 2:3) = [-1, -1];
    syapsyb(2, [1, 4]) = [1, -1];
    syapsyb(3, [1, 4]) = [1, -1];
    syapsyb(4, 2:3) = [1, 1];
    syapsyb = syapsyb*1i/2;
    
    % The rotation around z is not changing anything in case of S and/or T0
    % born rhoInit because szapszb commutes with any density matrix arising
    % from S and T0 states.
    rotz = expm(-1i*szapszb*(phiB));  % Rot 1: around z of minus phi (phi!)
    rotzT = expm(1i*szapszb*(phiB));
    roty = expm(-1i*syapsyb*(-thetaB)); % Rot 2: around y of minus theta
    rotyT = expm(1i*syapsyb*(-thetaB));

    rho1 = rotz*rhoInit*rotzT;
    % fprintf("%f", sum(rhoInit - rho1, 'all'))
    rho = roty*rho1*rotyT;
end

function rho = rotateproduct2coupled(rhoInit)
    %% Transform rho to the T+, S, T0, T- basis
    % This is the correct rotation matrix because:
    % rotProd2Coupled*[uu, ud, du, dd] = 
    % [uu, (ud - du)/sqrt(2), (ud + du)/sqrt(2), dd] = [Tp, S, T0, Tm]
    rotProd2Coupled = diag([1, 1/sqrt(2), 1/sqrt(2), 1]);
    rotProd2Coupled(2, 3) = -1/sqrt(2);
    rotProd2Coupled(3, 2) = 1/sqrt(2);

    % Matrices transform according to M' = R*M*R'
    rho = pagemtimes(rotProd2Coupled, ...
                pagemtimes(rhoInit, rotProd2Coupled'));
end

function U = diagonalizingmatrix(alpha, initBasis)
    arguments
        alpha
        initBasis = "product"
    end
    U = zeros(4, 4, numel(alpha));
    if initBasis == "coupled"
        % U = [1      0           0       0;
        %     0   cos(alpha)   sin(alpha)  0;
        %     0   -sin(alpha)  cos(alpha)  0;
        %     0       0          0         1];    
        U(1, 1, :) = 1;
        U(2, 2, :) = cos(alpha);
        U(2, 3, :) = sin(alpha);
        U(3, 2, :) = -sin(alpha);
        U(3, 3, :) = cos(alpha);
        U(4, 4, :) = 1;
    elseif initBasis == "product"
        % U = [1      0                         0               0;
        %     0   (cos + sin)/sqrt(2)   (sin - cos)/sqrt(2)     0;
        %     0   (cos - sin)/sqrt(2)   (cos + sin)/sqrt(2)     0;
        %     0       0                         0               1];    
        U(1, 1, :) = 1;
        U(2, 2, :) = (cos(alpha) + sin(alpha))/sqrt(2);
        U(2, 3, :) = (sin(alpha) - cos(alpha))/sqrt(2);
        U(3, 2, :) = -U(2, 3, :);
        U(3, 3, :) = U(2, 2, :);
        U(4, 4, :) = 1;
    end
end

function diag = mydiag3d(matr)
    
    % matr should be a (4, 4, nField) matrix.
    szmatr = size(matr);
    diag = zeros(szmatr(2:end));
    for ii = 1:4
        diag(ii, :) = matr(ii, ii, :);
    end
end