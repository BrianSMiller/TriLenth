% True location of source;
x = 0;
y = 0;

% Locations of observations (i.e. sonobuoy positions)
n = 2
xi = randn(n,1) + x;
yi = randn(n,1);
xi = [-1 0]';
yi = [0 1]';

% Simulate the observed bearings (in radians)
sigma = 0.1; % Standard deviation of bearings 

% Simulate the observed bearings (in radians)
theta = atan2(y-yi, x-xi) + randn(n,1) * sigma;

A = exp(-0.5*(sigma).^2);
kappaInv = 2*(1-A) + ((1-A).^2*(0.48794 - 0.82905 * A - 1.3915 * A.^2))/A;
kappa = 1/kappaInv;

[xHat yHat qHat] = lenth2(xi, yi, theta, kappa);
hMatlabMl = plot(xHat,yHat,'ks');
% hold on;
hError = error_ellipse(qHat,[xHat yHat],'conf',.95,'style','b--');

hBuoy = plot(xi, yi, 'o');
hold on;
hBearing = plot([xi xi + 10*cos(theta)]',[yi yi+10*sin(theta)]','k')
hTrue = plot(x, y,'rx');
[xHat yHat qHat kappa] = lenth(xi, yi, theta);
hEstimated = plot(xHat,yHat,'g^');
error_ellipse(qHat,[xHat yHat],'conf',.95,'style','g');
hCb = colorbar;
% ylabel','-Log-Likelihood');
legend([hTrue, hBearing(1), hEstimated, hError],'True Location','Measured Bearings (including error)','Estimated Location','Error Ellipse');


%%
