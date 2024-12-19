function outV = averageoversolidangle(inV, nTheta, nPhi, dim)
arguments
    inV
    nTheta
    nPhi
    dim
end

[thetas, ~] = createthetaphigrid(nTheta, nPhi);

if numel(dim) == 1
    thetas = repmat(thetas, [nPhi, 1]);
    weight = sin(thetas);
    weight = weight/sum(weight);

    weight = reshape(weight, [ones(1, dim - 1), nTheta*nPhi]);

    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
else
    weight = sin(thetas)/sum(sin(thetas)*ones(1, nPhi), 'all');

    weight = reshape(weight, [ones(1, min(dim(1) - 1, 1)), nTheta]);
    outV = sum(weight.*inV, dim);
    outV = squeeze(outV);
end
end

