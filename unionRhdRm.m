function [P, xInside, yInside] = unionRhdRm(x, y, fxy, fc, Rmx, Rmy)
% HDR2D computes a highest density contour.
%	Inputs: 
%   x - x-values where the density function is evaluated as a meshgrid.
%   y - y-values where the density function is evaluated as a meshgrid.
%   fxy - Density function at x,y.
%   Fxy - Cumulative distribution function at x, y.
%   fc - Density threshold.
%   Rmx - Mild region x-value boundary of polygon 
%   Rmy - Mild region y-value boundary of polygon 

% Calculate included probability as a function of density level.
dx = x(1, 2) - x(1, 1);
dy = y(2, 1) - y(1, 1);
isInHdr = fxy > fc;
isInRm = inpolygon(x, y, Rmx, Rmy);
xInside = x(isInHdr | isInRm);
yInside = y(isInHdr | isInRm);
fxyInside = fxy(isInHdr | isInRm);
P = sum(fxyInside) * dx * dy;
end