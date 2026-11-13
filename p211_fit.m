function [p, info] = p211_fit(d, opt)
%P211_FIT  Infer the ten transition constants for one experimental condition.
%
%   [p, info] = p211_fit(d)  with d one element of p211_data().
%
%   ROLE IN THE SHORT PAPER.  This is the random search of Sec. 2.3: a
%   coarse random sweep of the parameter space chooses the starting points
%   of the Nelder-Mead refinements.  It was used, with the original
%   (directional) targets of p211_data, in the searches that produced the
%   reference vectors stored in p211_fitted.mat.  Those vectors are the
%   best found over several such searches with different seeds and
%   budgets, all rescored on one common sample of 20,000 trajectories.
%   run_refit12 starts from them and refits under the corrected
%   observation rules (p211_refit).  run_short_all does not call this
%   function, and calling it does not regenerate p211_fitted.mat.
%
%   PARAMETRISATION.  The chain is a tree, so detailed balance fixes four
%   of the eight rates once the occupancies are known:
%
%       kpa = kap*piA/piP,   kpr = krp*piR/piP,
%       kdp = kpd*piP/piD,   ksd = kds*piD/piS.
%
%   The search therefore runs over the five occupancies (through a
%   softmax, so they stay positive and sum to one), four forward
%   timescales kap, krp, kpd, kds, and the two speeds.  That is the same
%   ten degrees of freedom as before, but every one of them is now a
%   quantity with a scale the data speaks to directly, which is why the
%   search starts close and converges in a few hundred evaluations rather
%   than a few thousand.
%
%   Ten constants are fitted to up to thirteen published statistics, but
%   the two three-category partitions (time budget; whole-track classes)
%   each sum to one, so the statistics carry at most eleven independent
%   constraints for APP and lysosomes and nine for Rab5, whose whole-track
%   classes were not reported.  APP and lysosomes are therefore barely
%   over-determined and Rab5 is under-determined by one: some direction in
%   its parameter space is not fixed by the data at all.  (These counts
%   refer to the directional targets of p211_data; with the pooled targets
%   of p211_data_pooled Rab5 has at most seven independent constraints.)
%
%   The objective is stochastic.  It is made deterministic for the
%   optimiser by common random numbers: p211_sim gets the same seed every
%   call, so the same realisations are reused and the search sees a smooth
%   surface.  The optimum is then re-evaluated on an independent, larger
%   sample to confirm it is not an artefact of that one stream.

if nargin < 2, opt = struct; end
if ~isfield(opt,'N'),      opt.N      = 1500;  end
if ~isfield(opt,'Npre'),   opt.Npre   = 600;   end
if ~isfield(opt,'npre'),   opt.npre   = 400;   end
if ~isfield(opt,'Nchk'),   opt.Nchk   = 20000; end
if ~isfield(opt,'maxfun'), opt.maxfun = 900;   end
if ~isfield(opt,'ntry'),   opt.ntry   = 2;     end
if ~isfield(opt,'quiet'),  opt.quiet  = false; end

so = struct('N',opt.N,'vthr',d.vthr,'seed',20261104);
sp = struct('N',opt.Npre,'vthr',d.vthr,'seed',20261104);

% ---- starting guess straight from the published numbers
piS = d.fracS;  if isnan(piS), piS = 0.10; end
piD = 0.05*(1-piS);
rest = 1 - piS - piD;
pi0 = [d.thetaA*rest, d.thetaR*rest, d.thetaP*rest, piD, piS];
pi0 = max(pi0, 1e-3);  pi0 = pi0/sum(pi0);

% engagement rate implied by the reversal frequency: inside a moving
% track a reversal needs one A->P->R round trip, so krev ~ thetaA*kap/2
kap0 = max(0.3, 2*(d.frev/10)/max(d.thetaA,0.05));

x0 = [ log(2.0*d.vsA), log(2.0*d.vsR), log(pi0), ...
       log(kap0), log(kap0), log(0.05), log(0.005) ];

% ---- global pre-search.  The simplex on its own settles in the slow,
%      long-run corner of the space and never reaches the fast-switching
%      regime that a reversal frequency of several per ten seconds
%      demands, so the starting point is chosen by a coarse random sweep
%      first.  Half the draws perturb the data-informed guess, half are
%      drawn from the full box.
rs = RandStream('mt19937ar','Seed',4242);
lo = [log(0.2) log(0.2) -4*ones(1,5) log(0.05) log(0.05) log(2e-3) log(2e-5)];
hi = [log(5.0) log(5.0)  2*ones(1,5) log(8.00) log(8.00) log(1.00) log(2e-1)];

cand = zeros(opt.npre, numel(x0));
for i = 1:opt.npre
    if mod(i,2) == 1
        s = 0.2 + 1.1*rand(rs);
        cand(i,:) = x0 + s*randn(rs,size(x0));
    else
        cand(i,:) = lo + (hi-lo).*rand(rs,size(x0));
    end
end
cand(1,:) = x0;

fpre = inf(opt.npre,1);
for i = 1:opt.npre
    fpre(i) = p211_resid(cand(i,:), d, sp);
end
[~, ord] = sort(fpre);

o = optimset('Display','off','MaxFunEvals',opt.maxfun, ...
             'MaxIter',opt.maxfun,'TolX',1e-4,'TolFun',1e-5);
obj = @(x) p211_resid(x, d, so);

bestf = inf; bestx = cand(ord(1),:);
for trial = 1:opt.ntry
    [xh, fh] = fminsearch(obj, cand(ord(trial),:), o);
    [xh, fh] = fminsearch(obj, xh, o);
    if fh < bestf, bestf = fh; bestx = xh; end
end

p = p211_unpack(bestx);

% ---- independent check at the optimum, fresh and larger sample
so2 = so; so2.N = opt.Nchk; so2.seed = 99001;
m = p211_sim(p, so2);

info.f      = bestf;
info.model  = m;
info.occ    = m.occ;
info.occAn  = p211_occ(p);
info.veff   = m.veff;
info.table  = p211_compare(d, m);
r = info.table.rel(~isnan(info.table.rel));
info.rmsRel = sqrt(mean(r.^2));
info.nTarget = numel(r);

if ~opt.quiet
    fprintf('\n%-14s  objective %.4g   rms relative residual %.1f%%  (%d targets)\n', ...
            d.name, bestf, 100*info.rmsRel, info.nTarget);
    fprintf('  va %.3f  vr %.3f  um/s\n', p.va, p.vr);
    fprintf('  kap %.4f  kpa %.4f  krp %.4f  kpr %.4f  1/s\n', ...
            p.kap, p.kpa, p.krp, p.kpr);
    fprintf('  kpd %.4f  kdp %.4f  kds %.5f  ksd %.5f  1/s\n', ...
            p.kpd, p.kdp, p.kds, p.ksd);
    fprintf('  occupancy A %.3f  R %.3f  P %.3f  D %.3f  S %.3f\n', info.occAn);
    fprintf('  net drift v_eff = %+.4f um/s\n', m.veff);
end
end

% ------------------------------------------------------------------------
function f = p211_resid(x, d, so)
p = p211_unpack(x);
if ~isfinite(p.va) || ~isfinite(p.vr) || p.va <= 0 || p.vr <= 0
    f = 1e6; return
end
m = p211_sim(p, so);
t = p211_compare(d, m);
r = t.rel(~isnan(t.rel));
if isempty(r), f = 1e6; return, end
f = sum(r.^2);
% a model that classifies nothing as moving tells us nothing
if isnan(m.thetaA), f = 1e6; end
if ~isfinite(f), f = 1e6; end
end

% ------------------------------------------------------------------------
function p = p211_unpack(x)
p.va = exp(x(1));
p.vr = exp(x(2));

w  = exp(x(3:7) - max(x(3:7)));
pi0 = w/sum(w);
pi0 = max(pi0, 1e-6);  pi0 = pi0/sum(pi0);

p.kap = exp(x(8));
p.krp = exp(x(9));
p.kpd = exp(x(10));
p.kds = exp(x(11));

% detailed balance on the tree closes the remaining four
p.kpa = p.kap*pi0(1)/pi0(3);
p.kpr = p.krp*pi0(2)/pi0(3);
p.kdp = p.kpd*pi0(3)/pi0(4);
p.ksd = p.kds*pi0(4)/pi0(5);
p.pi  = pi0;
end
