function fxy = pdfU10Hs(u, hs)
% pdfHsTp returns the density values of a the joint PDF of u and hs.
%	Inputs: 
%   u - 10-minute mean wind speed.
%   hs - Significant wave height.
%
% This model was presented in the OMAE 2020 paper "Global hierarchical 
% models for  wind and wave contours". It represents wind and wave at 
% the North Sea, at the FINO1 platform.

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

fxy = nan(size(u));
for i = 1 : size(u, 1)
    for j = 1 : size(u, 2)
        U = ExponentiatedWeibull(alpha, beta, delta);
        fu = U.pdf(u(i, j));
        betahs = c9 + c10 ./ (1 + exp(-1 * c11 * (u(i, j) - c12)));
        alphahs = (c6 + c7 * u(i, j).^c8) ./ 2.0445.^(1 ./ betahs);
        HS = ExponentiatedWeibull(alphahs, betahs, 5);
        fhsbaru = HS.pdf(hs(i, j));
        fxy(i, j) = fu * fhsbaru;
        %disp(['fxy(i,j) = ' num2str(fxy(i, j)) ' for i=' num2str(i) ', j=' num2str(j)]);
    end
end
% ExponentiatedWeibull.pdf(0) = NaN, but we want it to be 0.
fxy(isnan(fxy)) = 0;
end