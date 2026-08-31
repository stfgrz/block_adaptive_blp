function x = draw_ig(shape, rate)
% PURPOSE
% -------
% Draw one variate from the Inverse-Gamma distribution IG(shape, rate),
% parameterised so that
%     p(x) = rate^shape / Gamma(shape) * x^(-shape-1) exp(-rate/x),  x>0.
% This is the parameterisation used in every conditional posterior of
% samplers/gibbs_block_horseshoe.m.
%
% MODEL / EQUATIONS
% -----------------
% If X ~ IG(a, b) then 1/X ~ Gamma(a, rate = b) = Gamma(a, scale = 1/b).
% Hence: draw g ~ Gamma(a, scale = 1), return x = b / g.
%
% INPUTS
% ------
% shape : scalar > 0 (a above).
% rate  : scalar > 0 (b above).
%
% OUTPUTS
% -------
% x : scalar IG(shape, rate) draw. E[x] = rate/(shape-1) for shape > 1.
%
% DIMENSIONS
% ----------
% Scalars only.
%
% NOTES
% -----
% Uses utils/draw_gamma.m; no toolbox needed; reproducible under rng.

assert(shape > 0 && rate > 0, 'draw_ig: shape and rate must be > 0.');
g = draw_gamma(shape, 1);
x = rate / g;
end
