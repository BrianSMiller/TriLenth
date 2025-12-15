% Test script for functions triangulation functions lenth and lenth2. The
% script will simulate a source at location (x,y) and then simulate
% bearings to that source from receivers (sonobuoys) located at (xi,yi).
% The bearings will have simulated error of (sigma) radians. Functions
% lenth and lenth 2 take the sonobuoy locations and bearings+error and
% estimate the source location using "Lenths method" (1981) On Finding the
% Source of a Signal. Technometrics Vol. 23(2):149-154.
% 
% The difference between lenth and lenth2 is that lenth2 can be used when
% you have only 2 bearings and an can estimate their errors. Otherwise the
% function "lenth" can be used and will provide an estimate of kappa
% along with the estimated source location.

% True location of source;
x = 3;
y = 3;

% Locations of observations (i.e. sonobuoy positions)
xi = [-1 0 1]';
yi = [0 1 0]';
n = length(xi);
% xi = randn(n,1) + x;
% yi = randn(n,1);
% Simulate the observed bearings (in radians)
sigma = 0.1; % Standard deviation of bearings 

% Simulate the observed bearings (in radians)
theta = atan2(y-yi, x-xi) + randn(n,1) * sigma;

% The next three lines calculate kappa, which is required for the function
% lenth2. 
A = exp(-0.5*(sigma).^2);
kappaInv = 2*(1-A) + ((1-A).^2*(0.48794 - 0.82905 * A - 1.3915 * A.^2))/A;
kappa = 1/kappaInv;

% Triangulate using both methods
[xHat yHat qHat kap] = lenth(xi,yi,theta);
[xHat2 yHat2 qHat2] = lenth2(xi, yi, theta, kappa);


%% Plot the results, including error ellipse;

% Buoys and bearings, and source location;
hBuoy = plot(xi, yi, 'ko','markerFaceColor','k');
hold on;
hBearing = plot([xi xi + 10*cos(theta)]',[yi yi+10*sin(theta)]','k');
hTrue = plot(x, y,'ro','markerFaceColor','r');

% Lenth2 location and error ellipse
hLenth2 = plot(xHat2,yHat2,'bs');
hErrorLenth2 = error_ellipse(qHat2,[xHat2 yHat2],'conf',.95,'style','b--');

% Lenth location and error ellipse
hLenth = plot(xHat,yHat,'g^');
hErrorLenth = error_ellipse(qHat,[xHat yHat],'conf',.95,'style','g');

% hCb = colorbar;
% ylabel','-Log-Likelihood');
legend([hTrue, hBuoy, hBearing(1), hLenth, hErrorLenth, hLenth2, hErrorLenth2],...
   {'True Location','Buoys','Bearings (including error)',...
    'Lenth Location','Lenth Error Ellipse','Lenth2',' Lenth2 error'},...
    'location','northwest');
axis equal;
xlabel('X location (km)');
ylabel('Y location (km)');
hold off;


