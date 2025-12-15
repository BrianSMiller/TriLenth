
function [qHat] = lenthCov(xi, yi, theta, xHat, yHat, kappa)
% Estimate the covariance matrix qHat from Lenth (1981) On Finding the
% Source of a Signal. Technometrics. 23(2):149-154. Equation 2.10.
% Brian S. Miller, Australian Antarctic Division brian.miller@aad.gov.au
ci = cos(theta);
si = sin(theta);

dx = xHat - xi;
dy = yHat - yi;

d = sqrt(dx.^2 + dy.^2);
ciStar = dx./d.^3;
siStar = dy./d.^3;

% Estimate the covariance matrix for the MLE position.
qHat = [sum(siStar.*si),                     -0.5*sum(siStar.*ci+ciStar.*si);
        -0.5 * sum(siStar.*ci + ciStar.*si), sum(ciStar.*ci)                ];

qHat = 1./kappa .* inv(qHat);