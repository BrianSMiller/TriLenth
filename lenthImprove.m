function [xImproved yImproved] = lenthImprove(xEst, yEst, ci, si, xi, yi)
% Compute improved estimate of xHat and yHat by solving equation 2.6 from
% Lenth (1981). Improved estimates will eventually converge upon the MLE
% of xHat and yHat by iteratively calling this function using the improved
% positions as the estimated positions.
%
% Brian S. Miller, Australian Antarctic Division brian.miller@aad.gov.au

if (isempty(xEst)) % Initial estimate
    ciStar = ci;
    siStar = si;
else
    dx = xEst - xi;
    dy = yEst - yi;
    d3 = sqrt(dx.^2 + dy.^2).^3;
    ciStar = dx./d3;
    siStar = dy./d3;
end
z = si.*xi - ci.*yi;

A = [sum(si .*siStar) -sum(ci.*siStar);     % Equation 2.6,left hand side
    -sum(si .*ciStar)  sum(ci.*ciStar)];
b = [sum(siStar .* z); -sum(ciStar .*z)];   % Equation 2.6 right hand side

imp = A \ b;
xImproved = imp(1);
yImproved = imp(2);
end

