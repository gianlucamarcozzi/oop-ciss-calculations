function [APlus, AMinus] = calculateaplusaminusnohfi(A, AFrame, nVers)
    % For spin 1
    Ahfi1Tensor = rotatematrixeuangles(diag(A(1, 1:3)), AFrame(1, 1:3));
    Ahfi1 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi1Tensor', nVers)).^2)));
    % For spin 2
    Ahfi2Tensor = rotatematrixeuangles(diag(A(1, 4:6)), AFrame(1, 4:6));
    Ahfi2 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi2Tensor', nVers)).^2)));

    APlus = (Ahfi1 + Ahfi2)/2;  % MHz
    AMinus = (Ahfi1 - Ahfi2)/2;  % MHz
end