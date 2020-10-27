function [fxybar, x, y] = computefxybar(x, y, funfxy)
% Computes fxybar as outlined in XXX
%	Inputs: 
%   x - x-values where the density function is evaluated as a meshgrid.
%   y - y-values where the density function is evaluated as a meshgrid.
%   fxy - Density function at x,y.
%	Outputs: 
%   fxybar - Average density value in the cell centered at x,y.
%   x - x-values of the cell's center.
%   y - y-values of the cell's center.

if (x(1,1)  ~= 0) || (y(1,1)  ~= 0)
    error('Grid must start at (0,0)');
end

xx = x(1, :);
yy = y(:, 1);
for i = 2 : size(x, 2)
    for j = 2 : size(x, 1)
        Fji = integral2(funfxy, 0, x(1, i), 0, y(j), 'Method', 'iterated', 'RelTol', 1e-16);
        disp(['i,j=' num2str(i) ',' num2str(j) ', Fji=' num2str(Fji)]);
        F(j, i) = Fji;
    end
end

disp(['1 - the full integral should be 0 and was ' num2str(1 - max(max(F)))]);

cellArea = ((xx(1, 2) - xx(1, 1)) * (yy(2, 1) - yy(1, 1)));
for i = 1 : length(xx) - 1
    for j = 1 : length(yy) - 1
        cellProb(j, i) = F(j + 1, i + 1) - F(j, i + 1) - F(j + 1, i) + F(j, i);
    end
end

fxybar = cellProb / cellArea;

end