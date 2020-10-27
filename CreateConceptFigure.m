x = A.Tz;
y = A.Hs;

%[xnew, ynew] = meshgrid(linspace(0, 20, 100));
%znew = interp2(x, y, z, xnew, ynew, 'spline');

figure
scatter(x, y, 'k.')
hold on
bins = 50;
[n, c] = hist3([x, y], [bins bins]);
contour(c{1}, c{2}, n', [25, 50, 100, 200])
xlim([0 1.1 * max(x)]);
ylim([0 1.1 * max(y)]);
xlabel('Zero-up-crossing period (s)');
ylabel('Significant wave height (m)');
box off
