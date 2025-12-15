function [kappa kappaInv] = lenthKappa(xi, yi, theta, xHat, yHat)
% [kappa kappaInv] = lenthKappa(xi, yi, theta, xHat, yHat)
% Estimate the concentration parameter, kappa, from Lenth (1981) On Finding
% the Source of a Signal. Technometrics. 23(2):149-154. Equation 2.10.
%
% This function returns a value for kappa estimated from empirical data for
% a single localisation. 
%
% Kappa can alternatively be measured as the concentration parameter of a
% Von Mises distribution if the angular precisions (ie standard deviation)
% of each station is known. See: Mardia (1972) Statistics of 
% Directional Data. Academic Press. New York
%
% Brian S. Miller, Australian Antarctic Division brian.miller@aad.gov.au

dx = xHat - xi;
dy = yHat - yi;
muHat = atan2(dy,dx);
cBar = mean(cos(theta-muHat));
kappaInv = 2*(1-cBar) + (1-cBar)^2 .* (0.48794 - 0.82905 .* cBar - 1.3915.*cBar.^2)./cBar;
kappa = 1/kappaInv;