function dLs=lumdist_vec(zs,om,ol,w)

%compute luminosity distance in flat LCDM universe 
%by integration

%need to multiply by c/Ho to get distance in physical units

% z is allowed to be a sorted vector of redshifts
if ~issorted(zs)
    error('z are not sorted!');
end

if size(zs,2) > 1
    error('z must be a column vector of redshifts');
end

n_z = length(zs);

integrand = @(x) 1./(x.*sqrt(om./x - (om+ol-1) + ol*x.^(-(1+3*w)) ));

% lum dist = (1+z) x comoving transverse distance Dm

% d = (1+z)*quad(integrand,1/(1+z),1);

% Dm = quad(integrand,1/(1+z),1)

segments = nan(n_z,1);

segments(1) = quad(integrand,1/(1+zs(1)),1);

for i=2:n_z
   
    segments(i) = quad(integrand,1/(1+zs(i)),1/(1+zs(i-1)) );
    
end

Dms = cumsum(segments);
dLs = (1+zs).*Dms;

end

