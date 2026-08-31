function [k, theta] = gamma_coef(mode_val, sd_val)
% PURPOSE
% -------
% Convert a (mode, standard deviation) pair into the (shape k, scale
% theta) parameters of a Gamma distribution.  Direct port of
% GammaCoef.m from the Ferreira / Miranda-Agrippino / Ricco (FMAR)
% replication package; used for the Gamma hyperprior on the global
% tightness lambda in the FMAR marginal-likelihood selection.
%
% MODEL / EQUATIONS
% -----------------
% For X ~ Gamma(k, theta) with density p(x) ~ x^(k-1) exp(-x/theta):
%     mode = (k - 1) * theta   (k > 1),     var = k * theta^2.
% Solving the two equations for (k, theta) given (mode, sd) yields
%     k     = ( 2 + m^2/s^2 + sqrt( (4 + m^2/s^2) * m^2/s^2 ) ) / 2,
%     theta = sqrt( s^2 / k ),
% which is exactly the formula in FMAR's GammaCoef.m (and in the
% Giannone-Lenza-Primiceri 2015 code it derives from).
%
% INPUTS
% ------
% mode_val : scalar > 0, desired mode of the Gamma distribution.
% sd_val   : scalar > 0, desired standard deviation.
%
% OUTPUTS
% -------
% k, theta : shape and scale such that mode and sd are as requested.
%
% DIMENSIONS
% ----------
% Scalars only.
%
% NOTES
% -----
% The log-density used downstream is
%     log p(x) = (k-1) log x - x/theta - k log theta - gammaln(k).

assert(mode_val > 0 && sd_val > 0, 'gamma_coef: mode and sd must be > 0.');

ratio = mode_val^2 / sd_val^2;
k     = (2 + ratio + sqrt((4 + ratio) * ratio)) / 2;
theta = sqrt(sd_val^2 / k);
end
