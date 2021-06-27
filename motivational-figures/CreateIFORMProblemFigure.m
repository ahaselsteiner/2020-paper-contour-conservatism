% Figure showing cases where IFORM underpredicts alpha.

thisFolderName = '2020-paper-contour-conservatism';
addpath([thisFolderName '/compute-hdc'])
addpath([thisFolderName '/exponentiated-weibull'])
addpath([thisFolderName '/example1-subfunctions'])
lw = 1;


fig = figure('position', [100, 100, 1100, 340]);
layout = tiledlayout(1, 3);

% Sea state contour.
nexttile;

PM.modelType = 'CMA';
PM.distributions = {'weibull'; 'lognormal'};
PM.isConditionals = {[0 0 0]; [1 1]};
PM.coeffs = {{2.776 1.471 0.8888}; 
             { @(x1)0.1000 + 1.489 * x1^0.1901;
               @(x1)0.0400 + 0.1748 * exp(-0.2243*x1)}
            };
PM.labels = {'Significant wave height (m)';
              'Zero-upcrossing period (s)'};
alpha = 1 / (50 * 365.25 * 24 / 6);
[x1Array, x2Array] = computeIformContour(PM, alpha, 180, 0);
tzToTpFactor = 1.2796;


h_contour = plot([x2Array x2Array(1)] * tzToTpFactor, [x1Array x1Array(1)], ...
    '-b', 'linewidth', lw, 'DisplayName', 'IFORM contour');
hold on
hs = 0:0.2:20;
tp = 2:0.2:20;
[hgrid, tgrid] = meshgrid(hs, tp);
rgrid = responseTwoPeaks(hgrid, tgrid);
rcap = 16.6;
c = contour(tgrid, hgrid, rgrid, [rcap rcap], '-r');
c = c(:,2:end);
c = [c [tp(end); hs(end)]];
h_failuresurface = patch(c(1,:), c(2,:), [1 0.5 0.5], 'DisplayName', 'Failure region');
h_approximated = plot([0, 20], [18.3, 14.75], '--k', 'DisplayName', 'Approximated failure surface');
xlabel(PM.labels{2});
ylabel(PM.labels{1});

% Wind-wave contour.
nexttile;

PM.name = 'V-Hs model from X for dataset D (FINO1)';
PM.modelType = 'CMA';
PM.distributions = {'exponentiated-weibull'; 'exponentiated-weibull'};
PM.isConditionals = {[0 0 0 ]; [1 1 0]};
PM.coeffs = {
    {10.0 2.42 0.761}; 
    { 
    @(x1) (0.488 + 0.0114 * x1^2.03) / (2.0445^(1 / (0.714 + 1.70 / (1 + exp(-0.304 * (x1 - 8.77))))));
    @(x1) 0.714 + 1.70 / (1 + exp(-0.304 * (x1 - 8.77)));
    2}
    };
PM.labels = {'Wind speed (m s^{-1})';
    'Significant wave height (m)'};

hubHeight = 90;
jointModelHeight = 10;
alpha = 0.14; % as for normal wind speed in IEC 61400-3-1 p. 35
v_factor = (hubHeight / jointModelHeight)^alpha;
v_factor = v_factor * 0.95; % 61400-3-1:2019 p. 34


alpha = 1 / (50 * 365.25 * 24);
[x1Array, x2Array] = computeIformContour(PM, alpha, 180, 0);

plot([x1Array x1Array(1)] * v_factor, [x2Array x2Array(1)], '-b', 'linewidth', lw);
hold on
v = 0:0.25:40;
hs = 0:0.2:20;

[vgrid, hgrid] = meshgrid(v, hs);
sgrid = 0.012 + 0.021 ./ (1 + exp(-0.3*(vgrid - 10)));
tgrid = sqrt(2 * pi * hgrid ./ (9.81 * sgrid));
R = ResponseEmulator10mWaterDepth;
rgrid = R.ICDF1hr(vgrid, hgrid, tgrid, 0.5);
rcap = 113.3 * 10^6;
c = contour(vgrid, hgrid, rgrid, [rcap rcap], '-r');
c = c(:,2:end);
c = [[v(1); hs(end)] c [v(end); hs(end)]];
patch(c(1,:), c(2,:), [1 0.5 0.5]);
plot([14.9, 40], [20, 10.9], '--k');
xlabel(PM.labels{1});
ylabel(PM.labels{2});

% Hx-Hy contour.
nexttile;
hold on


% Define the parameters of the two considered response functions.
phi = 315 / 180 * pi;
a = 1.4;
b = 5.0;

hlim = 9;
hx = -hlim:0.2:hlim;
hy = -hlim:0.2:hlim;
[hxgrid, hygrid] = meshgrid(hx, hy);
rgrid = directionalResponse(hxgrid, hygrid, phi, a, b);
rcap = 10.31;
c = contour(hxgrid, hygrid, rgrid, [rcap rcap], '-r');
c = c(:,2:end);
patch([min(hx) max(hx) max(hx) min(hx)], [max(hy) max(hy) min(hy) min(hy)], [1 0.5 0.5]);
patch(c(1,:), c(2,:), [1 1 1]);
IFORM_hx_hy = readmatrix('iform_hx_hy.txt');
plot(IFORM_hx_hy(:,1), IFORM_hx_hy(:,2), '-b', 'linewidth', lw);
plot([-6.7, 9], [9, 1.5], '--k');
xlabel('x-component of significant wave height (m)');
ylabel('y-component of significant wave height (m)');

% Create a Legend with the data from multiple axes
lg = legend(nexttile(1), [h_contour, h_failuresurface, h_approximated], ...
     'Orientation', 'Horizontal');
lg.Layout.Tile = 'north';

exportgraphics(layout, 'IFORM_problem.pdf')

function r = directionalResponse(hx, hy, phi, a, b)
    hxPrime = hx * cos(phi) + hy * sin(phi);
    hyPrime = -1 * hx * sin(phi) + hy * cos(phi);
    r = sqrt(a * hxPrime.^2 + b * hyPrime.^2);
end
