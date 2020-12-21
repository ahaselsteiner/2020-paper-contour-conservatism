function [hs, tz] = drawSampleHsTz(n)
% drawSampleHsTp returns a random sample of the joint distrib. of hs and tz.
%	Inputs: 
%   n - length of sample.

r1 = rand(n, 1);
r2 = rand(n, 1);

alpha = 2.776;
beta = 1.471;
gamma = 0.8888;
a1 = 0.1;
a2 = 1.489;
a3 = 0.1901;
b1 = 0.04;
b2 = 0.1748;
b3 = -0.2243;

hs = wblinv(r1, alpha, beta) + gamma;

mu = a1 + a2 .* hs.^a3;
sigma = b1 + b2 .* exp(b3 .* hs);

tz = logninv(r2, mu, sigma);
end