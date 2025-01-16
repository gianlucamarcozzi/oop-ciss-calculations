function [APlus, AMinus] = calculateaplusaminusnohfi(A, AFrame, rVers)
    % For spin 1
    Ahfi1Tensor = rotatematrixeuangles(diag(A(1:3)), AFrame(1:3));
    Ahfi1 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi1Tensor', rVers)).^2)));
    % For spin 2
    Ahfi2Tensor = rotatematrixeuangles(diag(A(4:6)), AFrame(4:6));
    fprintf("Size Ahfi2Tensor: %d %d\n", size(Ahfi2Tensor))
    Ahfi2 = squeeze( ...
        sqrt( sum( (pagemtimes(Ahfi2Tensor', rVers)).^2)));
    fprintf("Size Ahfi2: %d %d\n", size(Ahfi2))

    APlus = (Ahfi1 + Ahfi2)/2;  % MHz
    AMinus = (Ahfi1 - Ahfi2)/2;  % MHz
end