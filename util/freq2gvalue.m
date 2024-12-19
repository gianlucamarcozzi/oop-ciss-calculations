function g = freq2gvalue(B0, freq)
% g = freq2gvalue(B0, freq)
%
% B0 in mT, freq in MHz
g = planck/bmagn/B0*freq*1e9;
end
