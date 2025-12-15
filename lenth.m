function [xHat, yHat, qHat, kappa] = lenth(xi, yi, theta)
% [xyMle qHat xSe ySe] = lenth(xi, yi, theta)
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
% This function applies the algorithm from Lenth (1981) On Finding the
% Source of a Signal. Technometrics Vol. 23(2):149-154.
%
% Specifically, this algorithm iteratively solves equation 2.6 "until none
% of the parameters change by more than a negligible amount."
%
% Brian S. Miller, Australian Antarctic Division brian.miller@aad.gov.au

maxIterations = 1e3; % Abort after this many iterations
threshold = 1e-4; % Define the "negligible amount" of change in parameters

% Convert bearings into direction cosines and sines.
ci = cos(theta);
si = sin(theta);
% cs = [ci si];

% Obtain initial estimate
[xHat yHat] = lenthImprove([], [], ci, si, xi, yi);

% Iterate until convergence or maximum number of iterations.
eps = inf;

while eps >= threshold
    [xUpdated yUpdated] = lenthImprove(xHat, yHat, ci, si, xi, yi);

    eps = norm([xHat-xUpdated yHat-yUpdated],2)/sqrt(norm([xUpdated yUpdated],2));
    xHat = xUpdated; yHat = yUpdated;

    maxIterations = maxIterations - 1;
    if (maxIterations < 0)
%         disp(sprintf('Convergence failed: epsilon=%g',eps));
        xHat = []; yHat = []; qHat = []; kappa = [];
        return
    end 
end %while

[kappa kappaInv] = lenthKappa(xi, yi, theta, xHat, yHat);   % Equation 2.10
qHat = lenthCov(xi, yi, theta, xHat, yHat, kappa);          % Equation 2.9




