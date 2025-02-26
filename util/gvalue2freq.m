function freq = gvalue2freq(B0, g)
% freq = gvalue2freq(B0, g)
%
% B0 in mT, output freq in MHz
freq = bmagn*B0/planck.*g*1e-9;
end