function [result, x, y] = triangulationError(xi,yi, xLim, yLim, step, sigma, n)
% [result, x, y] = triangulationError(xi,yi, xLim, yLim, step, sigma, n)
% Estiamate the error surface around a sonobuoy array. The locations of the
% sonobuoys are (xi,yi) in kilometers in cartesian coordinates. The limits of 
% the simulation are xLim and yLim, also in km in cartesian coordinates.
% The simulation is conducted over a grid with a step size of 'step'. Sigma
% is the standard deviation of bearing measurements (usually 5 degrees for
% sonobuoys). The number of simulated localisations at each grid step is
% specified as n. 
%
% Several 'preset' array types can be specified by including the following
% strings as the only input to this function: 
% 'pair9nmi' - a pair of sonobuoys 9 nmi apart from each other, 
% 'enrich', - a triplet of sonobuoys spaced similar to that deployed during
%             the 2019 ENRICH Voyage
% 'equilateral' - a triplet of sonobuoys deployed in an equilateral triangle
%
% Brian S Miller - Australian Antarctic Division - brian.miller@aad.gov.au
%
nmi2meter = 1854; % Convert nmi to m;
arrayType = 'customArray';
plotLocalizations = false;

if nargin == 1;
    arrayType = xi;
end

% Setup receiver locations (xi, yi)
% Locations of observations (i.e. sonobuoy positions)
if nargin < 2
    if nargin < 1; 
        arrayType = 'enrich';
    end
    switch arrayType
        case 'pair9nmi' % two sonobuoys 9 nmi apart
            xi = [-4.5 4.5]' * nmi2meter; % receiver x coordinates
            yi = [0 0]' * nmi2meter; % receiver y coordinates
        case 'enrich' % ENRICH style triplet (6 nmi isocoles right triangle)
%             xi = [0 0 -6]' * nmi2meter;
%             yi = [0 -6 0]' * nmi2meter;
            xi = [-4.2 0 4.2]' * nmi2meter;
            yi = [0 4.4 0]' * nmi2meter;
        case 'equilateral'
            xi = [-3 0 3]' * nmi2meter;
            yi = [0 5.2 0]' * nmi2meter;
        otherwise
            xi = [0 9]' * nmi2meter; % receiver x coordinates
            yi = [0 0]' * nmi2meter; % receiver y coordinates
    end
end

numBuoys = length(xi); % One bearing for each receiver location

cenx = sum(xi)/(numBuoys);
ceny = sum(yi)/(numBuoys);

xi = xi - cenx;
yi = yi - ceny;


if nargin < 5
    % Loop over a grid of locations and simulate bearings to the source
    xLim = [-5e4 5e4]; % X boundaries of area
    yLim = [-5e4 5e4]; % Y boundaries of area (in m)
    step = 1e3;  % 1 km steps
end
if nargin < 6
    sigma = 5 * pi/180; % Standard deviation of bearings (in radians)
end
if nargin < 7
    n = 100; % 100 observations per grid location
end



nRow = length(xLim(1):step:xLim(2));
nCol = length(yLim(1):step:yLim(2));

% Precalculate angular concentration based on nominal sonobuoy angular
% accuracy
A = exp(-0.5*(sigma).^2);
kappaInv = 2*(1-A) + ((1-A).^2*(0.48794 - 0.82905 * A - 1.3915 * A.^2))/A;
kappa = 1/kappaInv;

result = nan(nRow,nCol); % preallocate result for speed;
total = length(result(:));
% x and y are true location of source
count = 1; % Total number of locations processed
colCount = 1;
strLen = 0;
for x = xLim(1):step:xLim(2) % Loop over x locations
    rowCount = 1;
    for y = yLim(1):step:yLim(2); % Loop over y locations
        rmsError = nan(n,1);
        for i = 1:n          % Loop n times at each location
          
            % Simulate the observed bearings (in radians)
            theta = atan2(y-yi, x-xi) + randn(numBuoys,1) * sigma;
            
%             [xHat, yHat, qHat] = lenth2(xi, yi, theta, kappa);
            [xHat, yHat] = lenth(xi, yi, theta);
            if isempty(xHat) | isempty(yHat); % Lenths version failed, so use the matlab nonlinear least squares optimizer instead
                [xHat, yHat, qHat] = lenth2(xi, yi, theta, kappa);
            end
            rmsError(i) = abs(xHat - x)+abs(yHat - y);

            if plotLocalizations
                hMatlabMl = plot(xHat,yHat,'ks');
                hold on;
                hError = error_ellipse(qHat,[xHat yHat],'conf',.95,'style','b--');
                                set(gca,'ylimmode','auto','xlimmode','auto')
                hBuoy = plot(xi, yi, 'o');
%                 hold on;
%                 hBearing = plot([xi xi + 10*cos(theta)]',[yi yi+10*sin(theta)]','k')
                hBearing = plot([xi abs(xi-x).*cos(theta)]',[yi abs(yi-y).*sin(theta)]','k')
                hTrue = plot(x, y, 'rx');
%                 [xHat yHat qHat kappa] = lenth(xi, yi, theta);
%                 hEstimated = plot(xHat,yHat,'g^');
                error_ellipse(qHat,[xHat yHat],'conf',.95,'style','g');
%                 hCb = colorbar;

                % ylabel','-Log-Likelihood');
                legend([hTrue, hBearing(1), hMatlabMl, hError],...
                    'True Location','Measured Bearings (including error)','Estimated Location','Error Ellipse',...
                    'location','southoutside');
                pause(0.1);
                hold off;
            end
        end
        fprintf('%s',repmat(char(8),1,strLen));
        txt = sprintf('%g of %g locations calculated',count,total);
        fprintf('%s',txt);
        strLen = length(txt);
        count = count + 1;
        result(rowCount,colCount) = mean(rmsError);
        rowCount = rowCount + 1;
    end

    colCount = colCount + 1;
end
fprintf('\n');
x = xLim(1):step:xLim(2);
y = yLim(1):step:yLim(2);
[X, Y] = meshgrid(x,y);
distance = abs(cenx-X)+abs(cenx-Y);
percentError = (result./distance)*100;
save(sprintf('%sSonobuoyAccuracy.mat',arrayType));
%%
load(sprintf('%sSonobuoyAccuracy.mat',arrayType));
cmap = brewermap(256,'RdBu');

figure('units','centimeters','paperPositionMode','auto','position',[7 8 8.5 6]);
colormap(cmap);
radius = 25;
[cx,cy] = circle(radius);
cenx = sum(xi)/(numBuoys * 1e3);
ceny = sum(yi)/(numBuoys * 1e3);
cx = cx + cenx;
cy = cy + ceny;

% Plot absolute error
imagesc(x/1e3,y/1e3,result/1e3)
set(gca,'fontsize',8);
cb = colorbar('fontsize',8);
set(gca,'clim',[0 100]);
cbarrow('up');
ylabel(cb,'RMS triangulation error (km)')
xlabel('X position (km)');
ylabel('Y position (km)');
hold on;
plot(xi/1e3, yi/1e3,'m.','markerFaceColor','m')

% hCircle = plot(cx,cy,'k','lineWidth',2);
% set(gca,'xtick',[-50:25:50],'ytick',[-50:25:50]);
grid on;
saveAsPng(sprintf('%sSonobuoyAccuracy',arrayType));


% Plot percentage error
figure('units','centimeters','paperPositionMode','auto','position',[7 8 8.5 6]);
colormap(cmap)
imagesc(x/1e3,y/1e3,percentError);
set(gca,'fontsize',8);
cb = colorbar('fontsize',8);
set(gca,'clim',[0 200]);
cbarrow('up');
ylabel(cb,'Percentage error')

hold on;
plot(xi/1e3, yi/1e3,'m.','markerFaceColor','m')
xlabel('X position (km)');
ylabel('Y position (km)');

[cx,cy] = circle(radius);
% plot(cx+cenx,cy+ceny,'k','lineWidth',2);
% set(gca,'xtick',[-50:25:50],'ytick',[-50:25:50]);
grid on;

saveAsPng(sprintf('%sSonobuoyPercentError',arrayType));
% keyboard;



function [x y] = circle(r)
% Return points x,y on the circle of radius r
theta = linspace(0,2*pi,90);
x = r*cos(theta);
y = r*sin(theta);