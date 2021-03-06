DO_USE_FXBAR = 1;

tr = 50; % Return period in years.
ts = 1; % Environmental state duration in hours.
alpha = 1 / (tr * 365.25 * 24 / ts);

thisFolderName = '2020-paper-contour-conservatism';
addpath([thisFolderName '/compute-hdc'])
addpath([thisFolderName '/exponentiated-weibull'])
addpath([thisFolderName '/example2-subfunctions'])

du = 0.05;
dhs = du;
cellSize = du * dhs;

% Results from different grid resolutions:
% du = [0.5,        0.1,        0.05,         0.01]; dhs = du;
% fm = [1.1159e-07, 1.1272e-07, 1.1332e-07,   1.1339e-07]

u = du/2 : du : 40;
hs = dhs/2 : dhs : 25;
[U, HS] = meshgrid(u, hs);
uCell = u(1 : end - 1) + (u(2) - u(1)) / 2;
hsCell = hs(1 : end - 1) + (hs(2) - hs(1)) / 2;
[UCELL, HSCELL] = meshgrid(uCell, hsCell);

disp('Evaluating fxy in the variable space ...')
fxy = pdfU10Hs(U, HS);
if DO_USE_FXBAR == 1
    Fx = @(x) (x>0) * (1 - exp( -1 .* (x ./ 10.0).^2.42)).^0.761;
    betahs = @(u) 0.714 + 1.70 ./ (1 + exp(-1 * 0.304 * (u - 8.77)));
    alphahs = @(u) (0.488 + 0.0114 * u.^2.03) ./ 2.0445.^(1 ./ betahs(u));
    Fygivenx = @(y, x) (y>0) * (1 - exp( -1 .* (y ./ alphahs(x)).^betahs(x))).^5;
    E = ExponentiatedWeibull;
    expWblPdf = @(x, alpha, beta, delta) E.pdf(x, alpha, beta, delta);
    funfxy = @(x, y) expWblPdf(x, 10.0, 2.42, 0.761) .* expWblPdf(y, alphahs(x), betahs(x), 5);
    
    Ftot = integral2(funfxy, 0, max(u), 0, max(hs), 'RelTol', 1e-16);
    disp(['In Example 2: 1 - the full integral2 should be 0 and was ' num2str(1 - Ftot)]);
    fxybar = estimatefxybar(U, HS, Fx, Fygivenx);
    disp(['In Example 2: 1 - the sum of all fxybar should be 0 and was ' num2str(1 - sum(sum(fxybar .* cellSize)))]);

    % To be compatible with the rest of the code
    fxy = fxybar;
    u = uCell;
    hs = hsCell;
    U = UCELL;
    HS = HSCELL;
end

%figure
%contour(u, hs, fxy, [10^-8 10^-7 10^-6 10^-5 10^-4 10^-4 10^-3 10^-2])

v_in = 3;
CC = contourc(u, hs, fxy, [10^-9 10^-9]);
CCpiece = [CC(1,CC(1,:) < v_in); CC(2,CC(1,:) < v_in)];
CCpiece = CCpiece(:, 2:end); % First value is an artifact
V50 = 28.6;
Hs50 = 13.87;
hsMedian = @(v) 0.488 + 0.0114 * v.^2.03;
polygonU_Mid  = [V50, floor(V50) : -1 : v_in];
polygonHs_Mid = hsMedian(polygonU_Mid);
polygonU =  [0 V50  polygonU_Mid   CCpiece(1, :)];
polygonHs = [0 0    polygonHs_Mid  CCpiece(2, :)];

[P, uInside, hsInside] = unionRhdRm(U, HS, fxy, 0, polygonU, polygonHs);
disp(['1 - integrating over the whole variable space should yield 0 and was ' num2str(1 - P)])

disp('Calculating the adjusted HD contour ...')
fun = @(f) unionRhdRm(U, HS, fxy, f, polygonU, polygonHs) - (1 - alpha);
startValue = 0.001;
fm = fzero(fun, startValue);
[P, uInside, hsInside] = unionRhdRm(U, HS, fxy, fm, polygonU, polygonHs);

% Compute the adjusted contour using the found fm value.
C = contourc(u, hs, fxy, [fm, fm]);
Cu = C(1,2:end);
Chs = C(2,2:end);

isInRm = inpolygon(Cu, Chs, polygonU, polygonHs);
CuRel = Cu(~isInRm);
ChsRel = Chs(~isInRm);
CuRel = [polygonU(1) polygonU(2) CuRel CCpiece(1, :) polygonU(1)];
ChsRel = [polygonHs(1) polygonHs(2) ChsRel CCpiece(2, :) polygonHs(1)];

% Now calculate a "normal HD contour" using the compute-hdc software
% package.
disp('Calculating the normal HD contour ...')
PM.name = 'V-Hs model from X for dataset D (FINO1)';
PM.modelType = 'CMA';
PM.distributions = {'exponentiated-weibull'; 'exponentiated-weibull'};
PM.isConditionals = {[0 0 0 ]; [1 1 0]};
PM.coeffs = {
    {10.0 2.42 0.761}; 
    { 
    @(x1) (0.488 + 0.0114 * x1^2.03) / (2.0445^(1 / (0.714 + 1.70 / (1 + exp(-0.304 * (x1 - 8.77))))));
    @(x1) 0.714 + 1.70 / (1 + exp(-0.304 * (x1 - 8.77)));
    5}
    };
PM.labels = {'Wind speed (m/s)';
    'Significant wave height (m)'};
PM.gridCenterPoints = {0.05:0.1:50; 0.05:0.1:30};
[fmnormal, uHd, hsHd] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);

% To achieve smoother contour coordinates, use Matlab's contour function.
% instead of using computeHDC's hsHd{1}, tzHd{1}.
C = contourc(u, hs, fxy, [fmnormal, fmnormal]);
uHd = C(1,2:end);
hsHd = C(2,2:end);

% Plot the mild region with a Monte Carlo sample.
figRM = figure();
hold on
[uMC, hsMC] = drawSampleU10Hs(100000);
scatter(uMC, hsMC, '.k')
pgon = polyshape(polygonU, polygonHs);
plot(pgon);
xlim([0 30]);
ylim([0 14]);
xlabel('Wind speed (m s^{-1})');
ylabel('Significant wave height (m)');
set(figRM, 'renderer', 'painter')


% Plot the two contours, the mild regions, V50 and Hs50.
fig = figure();
plot(pgon);
hold on
contour(u, hs, fxy, [fmnormal, fmnormal], '-b', 'linewidth', 2);
plot(CuRel, ChsRel, '-r', 'linewidth', 2)
plot(Cu, Chs, '--r')
UWbl = ExponentiatedWeibull(10.0, 2.42, 0.761);
V50 = UWbl.icdf(1 - alpha);
plot([V50, V50,], [0, max(hs)], '--k');
fun = @(hs) integral2(funfxy, 0, 50, 0, hs, 'RelTol', 1e-16) - (1 - alpha);
startValue = 20;
Hs50 = fzero(fun, startValue);
plot([0 max(u)], [Hs50, Hs50], '--k');
xlabel('Wind speed (m s^{-1})');
ylabel('Significant wave height (m)');
legend({'Mild region', 'Normal HD contour', ...
    'Adjusted HD contour', 'Boarder of HD region within the mild region', ...
    'Marginal return values V_{50} and H_{s50}'}, ...
    'location', 'southoutside', 'NumColumns', 2);
legend box off
xlim([0 35]);
ylim([0 18]);
title([num2str(tr) '-yr environmental contours (alpha=' num2str(alpha) ')']);
