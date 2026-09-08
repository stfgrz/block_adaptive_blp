function x = draw_gamma(shape, scale)
% PURPOSE
% -------
% Draw one random variate from the Gamma(shape, scale) distribution using
% only base-MATLAB rand/randn.  This avoids a dependence on the
% Statistics and Machine Learning Toolbox (gamrnd / randg).
%
% MODEL / EQUATIONS
% -----------------
% Density (shape a > 0, scale s > 0):
%     p(x) = x^(a-1) exp(-x/s) / (Gamma(a) s^a),  x > 0,
% so E[x] = a*s.  The sampler is the Marsaglia-Tsang (2000) squeeze
% method for a >= 1:
%     d = a - 1/3,  c = 1/sqrt(9d),
%     repeat: z ~ N(0,1), v = (1 + c z)^3 (reject if v <= 0), u ~ U(0,1);
%     accept if  u < 1 - 0.0331 z^4
%            or  log(u) < 0.5 z^2 + d - d v + d log(v);
%     return d*v*scale.
% For a < 1, use the boosting identity
%     Gamma(a) =d Gamma(a+1) * U^(1/a),   U ~ Uniform(0,1),
% which follows from the density transformation in Marsaglia-Tsang.
%
% INPUTS
% ------
% shape : scalar > 0.
% scale : scalar > 0 (so rate = 1/scale).
%
% OUTPUTS
% -------
% x : scalar Gamma(shape, scale) draw.
%
% DIMENSIONS
% ----------
% Scalars only. Callers loop; all loops in this project are tiny.
%
% NOTES
% -----
% Fully driven by rand/randn, hence reproducible under rng(seed).

% NOTE ON THE ARGUMENT CHECKS.  These were written with assert(), which
% in Octave is an m-file doing a great deal of work on every call: it
% cost about 100 microseconds per draw, and this function is called
% roughly 400,000 times per block-adaptive estimator run, so the checks
% alone accounted for most of the Monte Carlo's runtime.  The plain
% conditional below is the same check with the same error, at a
% negligible cost.  Nothing about the algorithm or the consumption of
% the random stream changes.
if ~(isscalar(shape) && shape > 0)
    error('draw_gamma: shape must be a positive scalar.');
end
if ~(isscalar(scale) && scale > 0)
    error('draw_gamma: scale must be a positive scalar.');
end

if shape < 1
    % Boost: draw Gamma(shape+1,1), multiply by U^(1/shape).
    x = draw_gamma(shape + 1, 1) * rand()^(1 / shape) * scale;
    return;
end

d = shape - 1/3;
c = 1 / sqrt(9 * d);
while true
    z = randn();
    v = (1 + c * z)^3;
    if v <= 0
        continue;                       % outside the support of v
    end
    u = rand();
    if u < 1 - 0.0331 * z^4             % fast squeeze acceptance
        x = d * v * scale;
        return;
    end
    if log(u) < 0.5 * z^2 + d - d * v + d * log(v)   % exact acceptance
        x = d * v * scale;
        return;
    end
end
end
