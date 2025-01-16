function Sys = addnucspin(Sys, nucStr, A0, A0Frame, nNuc)
% Sys = addnucspin(Sys, nucStr, A, nNuc)
%
% Input
%   Sys
%   nucStr
%   A0      (1, nElectrons x 3)
%   A0Frame (1, nElectrons x 3)
%   nNuc
%
% Output
%   Sys
%     .Nucs         of the form [nucStr, nucStr, nucStr, ...]
%     .Sys.A        (nNuc, nElectrons x 3)
%     .Sys.AFrame   (nNuc, nElectrons x 3)

Sys.Nucs = nucStr;
Sys.Nucs = repmat(append(Sys.Nucs, ', '), [1, nNuc]);
Sys.Nucs = Sys.Nucs(1:end - 2);

Sys.A = repmat(A0, [nNuc, 1]);
Sys.AFrame = repmat(A0Frame, [nNuc, 1]);

end