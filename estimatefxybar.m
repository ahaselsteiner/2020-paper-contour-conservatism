function [fxybar, x, y] = estimatefxybar(x, y, Fx, Fygivenx)
% Estimates fxybar as outlined in XXX
%	Inputs: 
%   x - x-values where the density function is evaluated as a meshgrid.
%   y - y-values where the density function is evaluated as a meshgrid.

%	Outputs: 
%   fxybar - Average density value in the cell centered at x,y.
%   x - x-values of the cell's center.
%   y - y-values of the cell's center.

xx = x(1, :);
yy = y(:, 1);
dx = xx(2) - xx(1);
dy = yy(2) - yy(1);
cellArea = dx * dy;

for i = 1 : length(xx) - 1
    for j = 1 : length(yy) - 1
        Px = Fx(xx(i) + 0.5 * dx) - Fx(xx(i) - 0.5 * dx);
        Py = Fygivenx(yy(j) + 0.5 * dy, xx(i)) - Fygivenx(yy(j) - 0.5 * dy, xx(i));
        cellProb(j, i) = Px * Py;
        %disp(['i,j=' num2str(i) ',' num2str(j) ', cellProbji=' num2str(cellProb(j, i))]);
    end
end

fxybar = cellProb / cellArea;

end