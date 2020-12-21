function [u, hs] = drawSampleU10Hs(n)
% drawSampleU10Hs returns a random sample of U10 and Hs.
%	Inputs: 
%   n - length of sample.

r1 = rand(n, 1);
r2 = rand(n, 1);

alpha = 10.0;
beta = 2.42;
delta = 0.761;
c6 = 0.488;
c7 = 0.0114;
c8 = 2.03;
c9 = 0.714;
c10 = 1.70;
c11 = 0.304;
c12 = 8.77;

U = ExponentiatedWeibull(alpha, beta, delta);
u = U.icdf(r1);
betahs = c9 + c10 ./ (1 + exp(-1 * c11 * (u - c12)));
alphahs = (c6 + c7 * u.^c8) ./ 2.0445.^(1 ./ betahs);

hs = nan(n, 1);
for i = 1 : n
        HS = ExponentiatedWeibull(alphahs(i), betahs(i), 5);
        hs(i) = HS.icdf(r2(i));
end
end