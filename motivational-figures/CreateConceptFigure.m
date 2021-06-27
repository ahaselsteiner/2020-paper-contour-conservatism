x = A.Tz(1:1:end);
y = A.Hs(1:1:end);

figure
binscatter(x, y, [200 200])
xlim([0 1.1 * max(x)]);
ylim([0 1.1 * max(y)]);
xlabel('Zero-up-crossing period (s)');
ylabel('Significant wave height (m)');
colormap(gca,'parula')
box off
colorbar off
