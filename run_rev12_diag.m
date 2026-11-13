function run_rev12_diag(outfile)
%RUN_REV12_DIAG  Diagnose the reference vectors (p211_fitted.mat) under corrected observation rules.
%
%   The reference vectors are the starting points of run_refit12, not the
%   adopted fits; this script shows why they were refitted.
%
%   (1) unclassified moving tracks (|net velocity| <= 0.01 um/s)
%   (2) objective with the original targets, with exhaustive (renormalised)
%       whole-track classes, and with pooled Fig. 6C/6D/7F targets
%   (3) open-population composition of the spatial model against the
%       stationary distribution of the chain, and the movie statistics when
%       the movie samples the open population instead
%   MATLAB R2026a.

if nargin < 1, outfile = 'run_rev12_diag_output.txt'; end
S = load('p211_fitted.mat');  P = S.R.P;  D0 = p211_data();  Dp = p211_data_pooled();
lam = log(2)/(3.63*3600);
fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('P211  reference vectors under corrected observation rules   %s   MATLAB %s\n', ...
   char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
o = struct('N',20000,'seed',99001);
obj = @(t) sum(t.rel(~isnan(t.rel)).^2);
rmsr = @(t) sqrt(mean(t.rel(~isnan(t.rel)).^2));
pr('\n(1,2) condition | unclassified | objective: original, exhaustive, pooled+exhaustive | rms %% (same)\n');
for k = 1:6
    o.vthr = D0(k).vthr;
    m  = p211_sim(P{k}, o);
    me = p211_sim(P{k}, setfield(o,'exhaustive',true)); %#ok<SFLD>
    t0 = p211_compare(D0(k), m);  t1 = p211_compare(D0(k), me);  t2 = p211_compare(Dp(k), me);
    pr('   %-12s  %.3f   %.4f %.4f %.4f   %.1f %.1f %.1f   (pooled model: L %.2f um, pauses %.2f /10 s)\n', ...
       D0(k).name, m.fracUncl, obj(t0), obj(t1), obj(t2), 100*rmsr(t0), 100*rmsr(t1), 100*rmsr(t2), me.Lall, me.fpAll);
end

pr('\n(3) stationary distribution of the chain vs open-population composition (t1/2 = 3.63 h)   [A R P D S]\n');
for k = 1:6
    pi0 = p211_occ(P{k});  c = p211_spatcomp(P{k}, lam);
    pr('   %-12s  chain %s   open %s\n', D0(k).name, sprintf('%.3f ', pi0), sprintf('%.3f ', c));
end
pr('\n   APP movies sampled from the open population (no refit): whole-track classes and objective\n');
for k = 1:2
    o.vthr = D0(k).vthr;  o.exhaustive = true;
    m1 = p211_sim(P{k}, o);
    m2 = p211_sim(P{k}, setfield(o,'init',p211_spatcomp(P{k}, lam))); %#ok<SFLD>
    t1 = p211_compare(Dp(k), m1);  t2 = p211_compare(Dp(k), m2);
    pr('   %-7s chain-sampled  A/S/R %.3f %.3f %.3f  obj %.4f | open-population  A/S/R %.3f %.3f %.3f  obj %.4f  (measured %.2f %.2f %.2f)\n', ...
       D0(k).name, m1.fracA, m1.fracS, m1.fracR, obj(t1), m2.fracA, m2.fracS, m2.fracR, obj(t2), D0(k).fracA, D0(k).fracS, D0(k).fracR);
    o = rmfield(o,'exhaustive');
end
pr('\n   slow rates of the reference vectors relative to lam = %.2e 1/s\n', lam);
for k = 1:6
    p = P{k};
    pr('   %-12s  kpd %.2e kdp %.2e kds %.2e ksd %.2e   (ksd/lam %.3g, kds/lam %.3g)\n', D0(k).name, p.kpd, p.kdp, p.kds, p.ksd, p.ksd/lam, p.kds/lam);
end
fclose(fid);
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
