function intensityResonance = intensityresonance(rhoInit, rVers, deltaw, J, d, itrans)
    
    % Rotate rhoInit depending on theta, phi
    rho0 = rotatezfirstysecond(rhoInit, rVers);

    % rho should be in the product basis
    rho1 = rotateproduct2coupled(rho0);
    % rho0 = rhoInit;
    
    alpha = 1/2*atan(deltaw./(J + d/2));
    
    transProb = [sin(alpha).^2/2; ...    % 1 to 2
                cos(alpha).^2/2; ...    ss % 3 to 4
                cos(alpha).^2/2; ...     % 1 to 3
                sin(alpha).^2/2];        % 2 to 4

    UMatrix = diagonalizingmatrix(alpha);
    rho = pagemtimes(UMatrix, pagemtimes(rho1, pagetranspose(UMatrix)));
    populs = mydiag3d(rho);
    populDiff = [populs(1, :) - populs(2, :); ...
                populs(3, :) - populs(4, :); ...
                populs(1, :) - populs(3, :); ...
                populs(2, :) - populs(4, :)];
    
    if abs(rVers(3)) < 0.1 && itrans == 1
%         nVers
%         rho0
%         rho1
%         trace(rho(:, :, 1))
%         populDiff(:, 1)
    end
    intensityResonance = -populDiff(itrans, :).*transProb(itrans, :);
end

function rho = rotatezfirstysecond(rhoInit, rVers)
    
    theta = atan(hypot(rVers(1), rVers(2))/rVers(3));
    if rVers(1) ~= 0 || rVers(2) ~= 0
        phi = atan(rVers(2)/rVers(1));
    else
        phi = 0;  % This prevents phi = atan(0/0) = NaN
    end
    
    szapszb = 1/2*diag([1, 0, 0, -1]);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    rotz = expm(1i*szapszb*(-phi));  % Rot 1: around z of minus phi (phi!)
    rotzT = expm(-1i*szapszb*(-phi));
    roty = expm(1i*syapsyb*(-theta)); % Rot 2: around y of minus theta
    rotyT = expm(-1i*syapsyb*(-theta));
    
    rho = rotyT*rhoInit*roty;
    rho = rotzT*rho*rotz;

end

function rho = rotateproduct2coupled(rhoInit)
    % Transform rho to the T+, S, T0, T- basis
    rotProd2Coupled = diag([1, 1/sqrt(2), 1/sqrt(2), 1]);
    rotProd2Coupled(2, 3) = 1/sqrt(2);
    rotProd2Coupled(3, 2) = -1/sqrt(2);

    rho = rotProd2Coupled'*rhoInit*rotProd2Coupled;
end

function U = diagonalizingmatrix(alpha)
    
    % U = [1      0           0       0;
    %     0   cos(alpha)   sin(alpha)  0;
    %     0   -sin(alpha)  cos(alpha)  0;
    %     0       0          0         1];

    U = zeros(4, 4, numel(alpha));
    U(1, 1, :) = 1;
    U(2, 2, :) = cos(alpha);
    U(2, 3, :) = sin(alpha);
    U(3, 2, :) = -sin(alpha);
    U(3, 3, :) = cos(alpha);
    U(4, 4, :) = 1;
     
end

function diag = mydiag3d(matr)
    
    diag = squeeze(matr(1, :, :));
    for ii = 1:4  % 4 = first and second dimension of matr
        diag(ii, :) = matr(ii, ii, :);
    end
end