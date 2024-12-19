function xi = mixinganglexi(deltaw, JJ, dd)
%MIXINGANGLEXI Summary of this function goes here
%   Detailed explanation goes here
% if size(deltaw) ~= size(dd)

    
xi = atan( (dd/2 + JJ) ./ deltaw);
end

