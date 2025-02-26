function [APlus, AMinus] = calculateaplusaminusnohfi(A, AFrame, rVers)
    %% For spin 1
    % Calculate the hyperfine tensor for the given A values and A frame
    Ahfi1Tensor = rotatematrixeuangles(diag(A(1:3)), AFrame(1:3));
    % Project along rVers
    Ahfi1 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi1Tensor, rVers)).^2)));

    %% For spin 2
    % Calculate the hyperfine tensor for the given A values and A frame
    Ahfi2Tensor = rotatematrixeuangles(diag(A(4:6)), AFrame(4:6));
    % Project along rVers
    Ahfi2 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi2Tensor, rVers)).^2)));

    %% Aplus and Aminus
    APlus = (Ahfi1 + Ahfi2)/2;  % MHz
    AMinus = (Ahfi1 - Ahfi2)/2;  % MHz
end