function fxy = pdfHsTp(hs, tz)
% pdfHsTp returns the density values of a the joint PDF of hs and tz.
%	Inputs: 
%   hs - Significant wave height.
%   tz - Zero-up-crossing period.

alpha = 2.776;
beta = 1.471;
gamma = 0.8888;
a1 = 0.1;
a2 = 1.489;
a3 = 0.1901;
b1 = 0.04;
b2 = 0.1748;
b3 = -0.2243;
fx = wblpdf(hs - gamma, alpha, beta);
mu = a1 + a2 .* hs.^a3;
sigma = b1 + b2 .* exp(b3 .* hs);
fybarx = lognpdf(tz, mu, sigma);
fxy = fx .* fybarx;
end