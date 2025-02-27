function outV = averageoversolidangle(inV, thetas, nPhis, dim)
arguments
    inV
    thetas
    nPhis
    dim
end

% [thetas, ~] = createthetaphigrid(nTheta, nPhi);
nTheta = numel(thetas);

if isscalar(dim)
    thetas = repmat(thetas, [nPhis, 1]);
    weight = sin(thetas);
    weight = weight/sum(weight, 'all');

    if dim ~= 1
        weight = reshape(weight, [ones(1, dim - 1), nTheta*nPhis]);
    end

    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
else
    weight = sin(thetas)/sum(sin(thetas)*ones(1, nPhis), 'all');

    weight = reshape(weight, [ones(1, min(dim(1) - 1, 1)), nTheta]);
    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
end
end

