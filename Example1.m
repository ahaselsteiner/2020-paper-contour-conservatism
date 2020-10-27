DO_USE_FXBAR = 1;
DO_FULL_SEA_STATE = 0;

hsThresholds = [2 : 2 : 14]; % Mild region's upper Hs thresholds.
tr = 50; % Return period in years.
ts = 6; % Sea state duration in hours.
alpha = 1 / (tr * 365.25 * 24 / ts);
tzToTpFactor = 1.2796;

thisFolderName = '2020-paper-contour-conservatism';
addpath([thisFolderName '/compute-hdc'])
addpath([thisFolderName '/example1-subfunctions'])

dhs = 0.1;
dtz = dhs;
% Results from different grid resolutions:
% with polygonHs = [0 0 2 2]; polygonTz = [0 20 20 0];
% dhs = [0.5, 0.1, 0.05, 0.01, 0.005, 0.002, 0.0001]; dhs = du;
% Results using fxy (and not fxybar):
% fc = [NaN, 2.7825e-05, NaN, 2.5129e-06, 9.5934e-07, 2.0294e-06, script throws error]
% Results using computefxybar:
% fc = [1.701e-06, 4.653370166672485e-06, X, X, X, X, X]
% Results using estimatefxybar:
% fc = [1.825e-06, 1.831e-06, 1.825e-06, X, X, X, X]

hs = dhs/2 : dhs : 30;
tz = dhs/2 : dtz : 25;
[HS, TZ] = meshgrid(hs, tz);
hsCell = hs(1 : end - 1) + (hs(2) - hs(1)) / 2;
tzCell = tz(1 : end - 1) + (tz(2) - tz(1)) / 2;
[HSCELL, TZCELL] = meshgrid(hsCell, tzCell);

fxy = pdfHsTp(HS, TZ);

if DO_USE_FXBAR == 1
    Fx = @(x) wblcdf(x - 0.8888, 2.776, 1.471);
    Fygivenx = @(y, x) logncdf(y, 0.1 + 1.489 .* x.^0.1901, 0.04 + 0.1748 .* exp(-0.2243 .* x));
    funfxy = @(x, y) wblpdf(x - 0.8888, 2.776, 1.471) .* ...
        lognpdf(y, 0.1 + 1.489 .* x.^0.1901, 0.04 + 0.1748 .* exp(-0.2243 .* x));
    
    Ftot = integral2(funfxy, 0, max(hs), 0, max(tz), 'Method', 'iterated', 'RelTol', 1e-16);
    disp(['In Example 1: 1 - the full integral should be 0 and was ' num2str(1 - Ftot)]);
    %fxybar = computefxybar(HS, TZ, funfxy);
    fxybar = estimatefxybar(HS, TZ, Fx, Fygivenx);

    % To be compatible with the rest of the code
    fxy = fxybar;
    hs = hsCell;
    tz = tzCell;
    HS = HSCELL;
    TZ = TZCELL;

    %figure
    %surf(TZCELL', HSCELL', fxybar')
    %xlabel('Zero-up-crossiong period (s)');
    %ylabel('Significant wave height (m)');
end

for i = 1 : length(hsThresholds)
    hsThreshold = hsThresholds(i); 
    polygonHs = [0 0 hsThreshold hsThreshold];
    polygonTz = [0 15 15 0];

    [P, hsInside, tzInside] = unionRhdRm(HS, TZ, fxy, 0, polygonHs, polygonTz);
    disp(['1 - integrating over the whole variable space should yield 0 and was ' num2str(1 - P)]);

    fun = @(f) unionRhdRm(HS, TZ, fxy, f, polygonHs, polygonTz) - (1 - alpha);
    fc = fzero(fun, 0.001);
    [P, hsInside, tzInside] = unionRhdRm(HS, TZ, fxy, fc, polygonHs, polygonTz);

    % Compute the adjusted contour using the found fc value.
    C = contourc(hs, tz, fxy, [fc, fc]);
    Chs = C(1,2:end);
    Ctz = C(2,2:end);
    isInRm = inpolygon(Ctz, Chs, polygonTz, polygonHs);
    ChsRel = Chs(~isInRm);
    CtzRel = Ctz(~isInRm);
    
    % Add the mild region's boundary.
    TzMildRUpper = [CtzRel(end) : -0.5 : polygonTz(4) polygonTz(4)];
    HsMildRUpper = zeros(1, length(TzMildRUpper)) + polygonHs(4);
    adjustedCtz = [polygonTz(1) polygonTz(2) polygonTz(3) CtzRel TzMildRUpper polygonTz(1)];
    adjustedChs = [polygonHs(1) polygonHs(2) polygonHs(3) ChsRel HsMildRUpper polygonHs(1)];

    % Now calculate a "normal HD contour" using the compute-hdc software
    % package.
    %
    % Define a joint distribution for the sea state. The used model was 
    % proposed by Vanem and Bitner-Gregersen, DOI: 10.1016/j.apor.2012.05.006 .
    PM.name = 'Vanem and Bitner-Gregerse (2012), DOI: 10.1016/j.apor.2012.05.006';
    PM.modelType = 'CMA';
    PM.distributions = {'weibull'; 'lognormal'};
    PM.isConditionals = {[0 0 0]; [1 1]};
    PM.coeffs = {{2.776 1.471 0.8888}; 
                                 { @(x1)0.1000 + 1.489 * x1.^0.1901;
                                   @(x1)0.0400 + 0.1748 * exp(-0.2243*x1)}
                                };
    PM.labels = {'Significant wave height (m)';
                                 'Zero-up-crossing period (s)'};
    PM.gridCenterPoints = {0:0.05:20; 0:0.05:18};
    [fcnormal, hsHd, tzHd] = computeHdc(PM, alpha, PM.gridCenterPoints, 0);
    %hsHd = hsHd{1};
    %tzHd = tzHd{1};
    C = contourc(hs, tz, fxy, [fcnormal, fcnormal]);
    hsHd = C(1,2:end);
    tzHd = C(2,2:end);
    
    % Computing the response for the contours.
    tpHd = tzToTpFactor * tzHd;
    rHdNormal = responseTwoPeaks(hsHd, tpHd);
    [maxRNormal, iMaxNormal] = max(rHdNormal);

    tpAdjusted = tzToTpFactor * adjustedCtz;
    rHdAdjusted = responseTwoPeaks(adjustedChs, tpAdjusted);
    [maxRAdjusted, iMaxAdjusted] = max(rHdAdjusted);
    

    fig = figure();
    %scatter(tzInside, hsInside, '.k');
    pgon = polyshape(polygonTz, polygonHs);
    plot(pgon);
    hold on
    xlabel('Zero-up-crossing period (s)');
    ylabel('Significant wave height (m)');

    % Plot the normal HD contour. To have a smoother plot, use Matlab's contour
    % function.
    %plot(tzHd, hsHd, '-b', 'linewidth', 2)
    %contour(tz, hs, fxy', [fcnormal, fcnormal], '-b', 'linewidth', 2);
    plot(tzHd, hsHd, '-b', 'linewidth', 2)

    % Plot the adjusted HD contour.

    plot(adjustedCtz, adjustedChs, '-r', 'linewidth', 2)
    plot(Ctz, Chs, '--r')
    
    plot(adjustedCtz(iMaxAdjusted), adjustedChs(iMaxAdjusted), 'xr', ...
    'markersize', 10, 'linewidth', 2);

    % Plot the univariate return value
    Hsn = wblinv(1 - alpha, 2.776, 1.471) + 0.8888;
    plot([0 max(tz)], [Hsn, Hsn], '--k');

    r = responseTwoPeaks(HS, TZ * tzToTpFactor);
    [C, h] = contour(tz, hs, r', [5 10 15 17 20 25 30 35 40], 'color', [0.5 0.5 0.5]);
    clabel(C, h, 'FontSize', 6, 'Color', [0.5 0.5 0.5])

    legend({'Mild region', 'Normal HD contour', ...
        'Adjusted HD contour', 'border of HD region within the mild region', ...
        'Maximum response', 'H_{s50}', 'Constant response'}, ...
        'location', 'southoutside', 'NumColumns', 2);
    legend box off

    xlim([0 20]);
    ylim([0 20]);
    title([num2str(tr) '-yr environmental contours (alpha=' num2str(alpha) ')']);

    disp('Evaluating the long term response CDF ...');
    pfNormal = 1 - longTermResponseCdfE1TwoPeaks(maxRNormal);
    pfAdjusted = 1 - longTermResponseCdfE1TwoPeaks(maxRAdjusted);

    if DO_FULL_SEA_STATE == 1
        disp('Computing the all sea states approach ...');
        responseAllSeaStates = allSeaStateApproachE1TwoPeaks(alpha);
    else
        responseAllSeaStates = 16.9967; % Previously calculated.
    end

    % Create a table with the main results.
    if i == 1
        method = {'Normal HD contour'; 'Adjusted HD contour'};
        mildRegionHsMax = [0; max(polygonHs)];
        maxHs = [max(hsHd); max(Chs)];
        maxResponse = [maxRNormal; maxRAdjusted];
        maxResponseHs = [hsHd(iMaxNormal); adjustedChs(iMaxAdjusted)];
        maxResponseTp = [tpHd(iMaxNormal); tpAdjusted(iMaxAdjusted)];
        consevatismFactorR = [maxRNormal / responseAllSeaStates; maxRAdjusted / responseAllSeaStates];
        consevatismFactorPf = [alpha / pfNormal; alpha / pfAdjusted;];
    else
        c = 1;
        method{i + c} = 'All sea states';
        mildRegionHsMax(i + c) = max(polygonHs);
        maxHs(i + c) = max(Chs);
        maxResponse(i + c)=  maxRAdjusted;
        maxResponseHs(i + c) = adjustedChs(iMaxAdjusted);
        maxResponseTp(i + c) = tpAdjusted(iMaxAdjusted);
        consevatismFactorR(i + c) = maxRAdjusted / responseAllSeaStates;
        consevatismFactorPf(i + c) = alpha / pfAdjusted;
    end
end
c = 2;
method{i + c} = 'All sea states';
mildRegionHsMax(i + c) = NaN;
maxHs(i + c) = NaN;
maxResponse(i + c)=  responseAllSeaStates;
maxResponseHs(i + c) = NaN;
maxResponseTp(i + c) = NaN;
consevatismFactorR(i + c) = 1;
consevatismFactorPf(i + c) = 1;

Table = table(method, mildRegionHsMax, maxHs, maxResponse, maxResponseHs, maxResponseTp, consevatismFactorR, consevatismFactorPf)

