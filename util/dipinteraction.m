function dd = dipinteraction(dip, nVers, zD)
% dd = dipinteraction(dip, nVers, zD)
% 
% Dipolar interaction
% Makes use of the fact that cos(thetaD) = dotProduct(B0, zD)
dd = dip*((sum(nVers.*zD)).^2 - 1/3);
end

