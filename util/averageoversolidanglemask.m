function outV = averageoversolidanglemask(inV, nTheta, nPhi, dim, mask)
arguments
    inV
    nTheta
    nPhi
    dim
    mask = []
end

mask = checkparams(inV, nTheta, nPhi, dim, mask);
sizeV = size(inV);
[thetas, ~] = createthetaphigrid(nTheta, nPhi);

if numel(dim) == 1    
    % Put dim as first dimension
    idxNonAvg = [1:(dim - 1), (dim + 1):numel(sizeV)];
    inV = permute(inV, [dim, idxNonAvg]);  % Size is (nTheta*nPhi, ...)
    
    % Create logical mask
    logicMask = zeros(1, nTheta*nPhi);
    logicMask(mask) = 1;
    logicMask = boolean(logicMask);
    
    % Create weight
    thetas = repmat(thetas, [nPhi, 1]);
    weight = sin(thetas);  % Size is (nTheta*nPhi, 1)
    weight = weight(logicMask)/sum(weight(logicMask));
    
    % Average
    outV = sum(weight.*inV(logicMask, :), 1);
    
    % Reshape
    if numel(size(outV)) > 2
        outV = reshape(outV, sizeV(idxNonAvg));
    end
else
%     weight = sin(thetas)/sum(sin(thetas)*ones(1, nPhi), 'all');
% 
%     weight = reshape(weight, [ones(1, min(dim(1) - 1, 1)), nTheta]);
%     outV = sum(weight(mask).*inV(mask), dim);
%     outV = squeeze(outV);
    error('Still needs to be implemented')
end
end

function mask = checkparams(inV, nTheta, nPhi, dim, mask)

% Check numel(dim)
switch numel(dim)
    case 1
        if size(inV, dim) ~= nTheta*nPhi
            error("size(inV, dim) must equal nTheta*nPhi.")
        end
    case 2
        if size(inV, dim(1)) ~= nTheta
            error("size(inV, dim(1)) must equal nTheta.")
        elseif size(inV, dim(2)) ~= nPhi
            error("size(inV, dim(2)) must equal nPhi.")
        end
    otherwise
        error("numel(dim) must be equal to 1 or 2.")
end

% Check mask
if isempty(mask)
    if numel(dim) == 1
        mask = boolean(ones(1, nTheta*nPhi));
    else
        mask = boolean(ones(nTheta, nPhi));
    end
end

end