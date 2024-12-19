function [gax, gInh] = creategaxisinhbroadening(g, glw, dgax)
% [gax, gInh] = creategaxisinhbroadening(g, glw, dgax)
% Here glw1 and glw2 are the lwpp, therefore the variance will be
% lwpp*2
gax = (g - glw):dgax:(g + glw);
    
% Inhomogeneously broadened gaussian distributions
normTerm = sqrt(2/pi)/glw;
gInh = normTerm*exp(-2 * (gax - g).^2 ./ glw^2);

end