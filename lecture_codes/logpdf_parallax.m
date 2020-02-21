function y = logpdf_parallax(rs,omega,somega,Lgal)

logpdf = rs*0 -500;

q = (rs > 0);

logpdf(q) = -0.5*log(2*pi*somega) -0.5*(omega - 1./rs).^2 /somega^2 + 2*log(rs) - rs/Lgal;
%normpdf(1./x,omega,somega).* x.^2 .* exp(-x/Lgal);

y = logpdf;

end