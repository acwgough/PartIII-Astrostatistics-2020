function d=lumdist(z,om,ol,w)

%compute luminosity distance in flat LCDM universe 
%by integration

%need to multiply by c/Ho to get distance in physical units

integrand = @(x) 1./(x.*sqrt(om./x - (om+ol-1) + ol*x.^(-(1+3*w)) ));

d = (1+z)*quad(integrand,1/(1+z),1);

end

