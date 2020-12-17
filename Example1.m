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

dhs = 0.05;
dtz = dhs;

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
    fxybar = estimatefxybar(HS, TZ, Fx, Fygivenx);

    % To be compatible with the rest of the code
    fxy = fxybar;
    hs = hsCell;
    tz = tzCell;
    HS = HSCELL;
    TZ = TZCELL;
end

%figure
%contour(tz, hs, fxy', [10^-10 10^-9 10^-8 10^-7 10^-6 10^-5 10^-4 10^-4 10^-3 10^-2])

CC = contourc(tz, hs, fxy', [10^-9 10^-9]);
for i = 1 : length(hsThresholds)
    hsThreshold = hsThresholds(i);  
    RM = [CC(1, CC(2,:) < hsThreshold); CC(2, CC(2,:) < hsThreshold)];
    polygonHs = RM(2, :);
    polygonTz = RM(1, :);
    [mx, idx] = max(polygonHs);
    polygonHs = flip(circshift(polygonHs, -1 * idx));
    polygonTz = flip(circshift(polygonTz, -1 * idx));
    if polygonTz(1) > 2 * polygonTz(2)
        polygonTz = circshift(polygonTz, -1);
        polygonHs = circshift(polygonHs, -1);
    end
    
    [P, hsInside, tzInside] = unionRhdRm(HS, TZ, fxy, 0, polygonHs, polygonTz);
    disp(['1 - integrating over the whole variable space should yield 0 and was ' num2str(1 - P)]);

    fun = @(f) unionRhdRm(HS, TZ, fxy, f, polygonHs, polygonTz) - (1 - alpha);
    fm = fzero(fun, 0.001);
    [P, hsInside, tzInside] = unionRhdRm(HS, TZ, fxy, fm, polygonHs, polygonTz);

    % Compute the adjusted contour using the found fm value.
    C = contourc(hs, tz, fxy, [fm, fm]);
    Chs = C(1,2:end);
    Ctz = C(2,2:end);
    isInRm = inpolygon(Ctz, Chs, polygonTz, polygonHs);
    ChsRel = Chs(~isInRm);
    CtzRel = Ctz(~isInRm);
    
    % Add the mild region's boundary.
    adjustedCtz = [CtzRel polygonTz CtzRel(1)];
    adjustedChs = [ChsRel polygonHs ChsRel(1)];

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
    
    % To achieve smoother contour coordinates, use Matlab's contour function.
    % instead of using computeHDC's hsHd{1}, tzHd{1}.
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
    
    % Plot the two contours, the mild regions and Hs50.
    fig = figure();
    pgon = polyshape(polygonTz, polygonHs);
    plot(pgon);
    hold on
    plot(tzHd, hsHd, '-b', 'linewidth', 2)
    plot(adjustedCtz, adjustedChs, '-r', 'linewidth', 2)
    plot(Ctz, Chs, '--r')
    plot(adjustedCtz(iMaxAdjusted), adjustedChs(iMaxAdjusted), 'xr', ...
    'markersize', 10, 'linewidth', 2);
    Hsn = wblinv(1 - alpha, 2.776, 1.471) + 0.8888;
    plot([0 max(tz)], [Hsn, Hsn], '--k');
    r = responseTwoPeaks(HS, TZ * tzToTpFactor);
    [C, h] = contour(tz, hs, r', [5 10 15 17 20 25 30 35 40], 'color', [0.5 0.5 0.5]);
    clabel(C, h, 'FontSize', 6, 'Color', [0.5 0.5 0.5])
    xlabel('Zero-up-crossing period (s)');
    ylabel('Significant wave height (m)');
    legend({'Mild region', 'Normal HD contour', ...
        'Adjusted HD contour', 'Boundary of HD region within mild region', ...
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
        mildRegionHsMax = [0; hsThreshold];
        maxHs = [max(hsHd); max(Chs)];
        maxResponse = [maxRNormal; maxRAdjusted];
        maxResponseHs = [hsHd(iMaxNormal); adjustedChs(iMaxAdjusted)];
        maxResponseTp = [tpHd(iMaxNormal); tpAdjusted(iMaxAdjusted)];
        consevatismFactorR = [maxRNormal / responseAllSeaStates; ...
            maxRAdjusted / responseAllSeaStates];
        consevatismFactorPf = [alpha / pfNormal; alpha / pfAdjusted;];
    else
        c = 1;
        method{i + c} = 'Adjusted HD contour';
        mildRegionHsMax(i + c) = hsThreshold;
        maxHs(i + c) = max(Chs);
        maxResponse(i + c)=  maxRAdjusted;
        maxResponseHs(i + c) = adjustedChs(iMaxAdjusted);
        maxResponseTp(i + c) = tpAdjusted(iMaxAdjusted);
        consevatismFactorR(i + c) = maxRAdjusted / responseAllSeaStates;
        consevatismFactorPf(i + c) = alpha / pfAdjusted;
    end
end
c = 2;
method{i + c} = 'IFORM contour';
mildRegionHsMax(i + c) = NaN;
maxHs(i + c) = 15.2; % see 10.1016/j.marstruc.2020.102863;
maxResponse(i + c)=  16.57; % see 10.1016/j.marstruc.2020.102863
maxResponseHs(i + c) = NaN;
maxResponseTp(i + c) = NaN;
consevatismFactorR(i + c) = 0.97; % see 10.1016/j.marstruc.2020.102863
consevatismFactorPf(i + c) = 0.61; % see 10.1016/j.marstruc.2020.102863

c = 3;
method{i + c} = 'All sea states';
mildRegionHsMax(i + c) = NaN;
maxHs(i + c) = NaN;
maxResponse(i + c)=  responseAllSeaStates;
maxResponseHs(i + c) = NaN;
maxResponseTp(i + c) = NaN;
consevatismFactorR(i + c) = 1;
consevatismFactorPf(i + c) = 1;

Table = table(method, mildRegionHsMax, maxHs, maxResponse, ...
    maxResponseHs, maxResponseTp, consevatismFactorR, consevatismFactorPf)

