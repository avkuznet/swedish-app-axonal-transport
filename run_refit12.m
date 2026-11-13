function run_refit12(k)
%RUN_REFIT12  Corrected refit of condition k (1..6).
%
%   All conditions: pooled targets where the source pools them
%   (p211_data_pooled) and exhaustive whole-track classes.
%   APP (k = 1, 2): in addition, movie initial states drawn from the open-
%   population composition at the in-axon half-life of 3.63 h, so that the
%   movies and the spatial model sample the same population.  An extra
%   start rescales the reference kds and ksd (ratio kept) to 20 lam.
%   Output p211_refit12_<k>.mat.  MATLAB R2026a.

S = load('p211_fitted.mat');  P = S.R.P;
Dp = p211_data_pooled();  d = Dp(k);
opt = struct('nstart',6,'maxfun',1500);
if k <= 2
    lam = log(2)/(3.63*3600);
    opt.lam = lam;
    q = P{k};  f = max(1, 20*lam/q.ksd);
    q.kds = q.kds*f;  q.ksd = q.ksd*f;
    opt.starts = {q};
end
t0 = tic;
R = p211_refit(d, P{k}, opt);
R.elapsed = toc(t0);
save(sprintf('p211_refit12_%d.mat', k), 'R');
fprintf('%s done in %.0f s: objective %.4f (reference vector under the same rules %.4f), rms %.1f%%\n', ...
        d.name, R.elapsed, R.fchk, R.f0, 100*R.rms);
end
