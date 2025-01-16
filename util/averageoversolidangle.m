function outV = averageoversolidangle(inV, thetas, nPhi, dim)
arguments
    inV
    thetas
    nPhi
    dim
end

% [thetas, ~] = createthetaphigrid(nTheta, nPhi);
nTheta = numel(thetas);

if isscalar(dim)
    thetas = repmat(thetas, [nPhi, 1]);
    weight = sin(thetas);
    weight = weight/sum(weight, 'all');

    weight = reshape(weight, [ones(1, dim - 1), nTheta*nPhi]);
    % size(weight)

    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
else
    weight = sin(thetas)/sum(sin(thetas)*ones(1, nPhi), 'all');

    weight = reshape(weight, [ones(1, min(dim(1) - 1, 1)), nTheta]);
    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
end
end

