function [newMatr, idxs1, idxs2] = restrictgaussianfactor(matr, xx1, xx2)
% [newMatr, idxs1, idxs2] = restrictgaussianfactor(matr, xx1, xx2)
    [~, idxMax1] = max(xx1);
    [~, idxMax2] = max(xx2);
    [~, idxFwhm1] = min(abs(matr(:, idxMax2) - max(matr(:))/2));
    idxDxFwhm1 = abs(idxMax1 - idxFwhm1);
    [~, idxFwhm2] = min(abs(matr(idxMax1, :) - max(matr(:))/2));
    idxDxFwhm2 = abs(idxMax2 - idxFwhm2);
    
    xxFactor = 2;
    dx1 = round(xxFactor*idxDxFwhm1);
    dx2 = round(xxFactor*idxDxFwhm2);
    idxs1 = max(idxMax1 - dx1, 1):...
            min(idxMax1 + dx1, numel(xx1));
    idxs2 = max(idxMax2 - dx2, 1):...
            min(idxMax2 + dx2, numel(xx1));
    newMatr = matr(idxs1, idxs2);
end
