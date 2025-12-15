function [xHat yHat qHat] = lenth2(xi, yi, theta, kappa)
% [xyMle qHat xSe ySe] = lenth2(xi, yi, theta, kappa)
% Estimate the position of a source given bearings from a known location.
% xi, yi are an n x 1 vectors of positions in eastings and northings
%   respectively.
% theta is a n x 1 vector of angles (not azimuths) to the
%   target. Azimuths in degrees (measured from true North) can be converted
%   to angles such that: theta = (90-az)*(pi/180).
% xHat, yHat are the maximum likelihood x and y positions of the source.
% kappa is a measure of the concentration of angles (inversely proportional
%   to bearing error)
% qHat is the covariance matrix which can be used to obtain the error
%   ellipse
%
% Algorithm
% From Lenth (1981) On Finding the Source of a Signal. Technometrics Vol.
% 23(2):149-154.
%
% Specifically, this algorithm uses matlab's fminsearch to minimise the
% negative log-likelihood (given by equation 2.2).
% 
% This function assumes that the concentration parameter kappa, has been
% determined empirically, thus the covariance matrix can be computed even
% when only two bearings are observed. 
% See lenth.m, if the concentration parameter is unknown, but 3 or more
% bearings are available for triangulation.
%
% Brian S. Miller, Australian Antarctic Division brian.miller@aad.gov.au
%
addpath('c:\analysis\stats\CircStat\')
displayLikelihoodSurface = false;

% Convert bearings into direction cosines and sines.
ci = cos(theta);
si = sin(theta);
xyReceiver = [xi yi];
% Compute initial estimate
[xInitial, yInitial] = lenthImprove([], [], ci, si, xi, yi);
xyInitial = [xInitial yInitial];

[v, fval] = fminsearch(@(v) -logl(v, xyReceiver, theta, kappa),...
    xyInitial,optimset('Display','off','MaxIter',1e5));

xHat = v(1);
yHat = v(2);
[qHat] = lenthCov(xi, yi, theta, xHat, yHat, kappa);

if displayLikelihoodSurface;
    xyMin = min(xyReceiver);
    xyMax = max(xyReceiver);
    xyDist = norm(xyMax-xyMin);
    lowerBound = xyMin - 1 * xyDist;
    upperBound = xyMax + 1 * xyDist;

    n = 100;
    x = linspace(lowerBound(1),upperBound(1),n);
    y = linspace(lowerBound(2),upperBound(2),n);
    for i = 1:n;
        for j = 1:n;
            v = [x(i) y(j)];
            L(i,j) = logl(v, xyReceiver, theta, kappa);
        end
    end
    % Find the peak of the likelihood surface
    [Lmin ix] = max(L(:));
    [i j] = ind2sub(size(L),ix);

    % Show the likelihood surface. pcolor requires a matrix
    if length(x) > 1 && length(y) > 1;
        hSurf = pcolor(x,y,L');
        set(hSurf,'lineStyle','none');
        set(gca,'ydir','normal');
        hold on;
         plot(x(i),y(j),'ks');
%                 plot(phi,r,'bs');
        xlabel('Easting (m)');
        ylabel('Northing (m)');
        %         [C h] = contour(x,y,L',(max(L(:))-6):max(L(:)),'k-');
        %         clabel(C,h);
        colormap(autumn(16));
    end
end

function L = logl(v,  xyReceiver, theta, kappa)
n = length(theta);
xi = xyReceiver(:,1);
yi = xyReceiver(:,2);
mu = expectedBearing(v(1), v(2), xi, yi);
err = circ_dist(theta,mu);           
L = n * log(besseli(0,kappa)) + kappa.*sum(cos(theta-mu));
% L = sum(dnorm(err, 0, kappa, true)); % Use normal distribution instead of bessel function?

function mu = expectedBearing(x, y, xi, yi)
mu = atan2(y-yi,x-xi);

function density = dnorm(x, m, sd, logFlag)
density = 1/(sqrt(2 * pi) * sd) * exp(-((x - m).^2)/(2 * sd.^2));
if logFlag
    density = log(density);
end