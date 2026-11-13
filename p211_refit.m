function R = p211_refit(d, p0, opt)
%P211_REFIT  Warm-started refit of one condition under a specified observation model.
%
%   R = p211_refit(d, p0, opt)
%     d    one element of p211_data_pooled (or p211_data)
%     p0   starting parameter struct (e.g. the reference vector)
%     opt.lam        if set, movie initial states are drawn from the open-
%                    population composition p211_spatcomp(p, lam) instead of
%                    the stationary distribution of the chain
%     opt.exhaustive renormalise whole-track classes (default true)
%     opt.nstart     perturbed restarts in addition to p0 (default 6)
%     opt.maxfun     simplex evaluations per pass (default 1500)
%     opt.N          trajectories per objective evaluation (default 1500)
%     opt.starts     extra starting structs (cell), e.g. rescaled slow rates
%
%   Same ten-constant parametrisation and common-random-number objective as
%   p211_fit.  Every candidate is re-scored on an independent sample of
%   20,000 trajectories; R.cand lists them all, R.p is the best.

if nargin < 3, opt = struct; end
if ~isfield(opt,'exhaustive'), opt.exhaustive = true; end
if ~isfield(opt,'nstart'), opt.nstart = 6;    end
if ~isfield(opt,'maxfun'), opt.maxfun = 1500; end
if ~isfield(opt,'N'),      opt.N      = 1500; end
if ~isfield(opt,'starts'), opt.starts = {};   end
if ~isfield(opt,'lam'),    opt.lam    = [];   end

so = struct('N',opt.N,'vthr',d.vthr,'seed',20261104,'exhaustive',opt.exhaustive);
obj = @(x) resid(x, d, so, opt.lam);

X0 = pack(p0);
for j = 1:numel(opt.starts), X0(end+1,:) = pack(opt.starts{j}); end %#ok<AGROW>
rs = RandStream('mt19937ar','Seed',777);
base = X0;
for j = 1:opt.nstart
    X0(end+1,:) = base(1 + mod(j-1, size(base,1)),:) + 0.3*randn(rs, 1, size(base,2)); %#ok<AGROW>
end

o = optimset('Display','off','MaxFunEvals',opt.maxfun,'MaxIter',opt.maxfun,'TolX',1e-4,'TolFun',1e-5);
nc = size(X0,1);  R.cand = struct('p',{},'f',{},'fchk',{},'m',{});
for j = 1:nc
    [xh, fh] = fminsearch(obj, X0(j,:), o);
    [xh, fh] = fminsearch(obj, xh, o);
    p = unpack(xh);
    [fchk, m] = score(p, d, opt);
    R.cand(j) = struct('p',p,'f',fh,'fchk',fchk,'m',m);
    fprintf('  %-12s start %2d  objective %.4f (check %.4f)\n', d.name, j, fh, fchk);
end
[~, ib] = min([R.cand.fchk]);
R.p = R.cand(ib).p;  R.fchk = R.cand(ib).fchk;  R.m = R.cand(ib).m;
[R.f0, R.m0] = score(p0, d, opt);            % the starting vector, same scoring
R.table = p211_compare(d, R.m);
r = R.table.rel(~isnan(R.table.rel));  R.rms = sqrt(mean(r.^2));
R.opt = opt;  R.name = d.name;
end

% ------------------------------------------------------------------------
function [f, m] = score(p, d, opt)
o = struct('N',20000,'vthr',d.vthr,'seed',99001,'exhaustive',opt.exhaustive);
if ~isempty(opt.lam), o.init = p211_spatcomp(p, opt.lam); end
m = p211_sim(p, o);
t = p211_compare(d, m);
f = sum(t.rel(~isnan(t.rel)).^2);
m = rmfield(m, 'X');        % the simulated trajectories are not stored with the fit
end

function f = resid(x, d, so, lam)
p = unpack(x);
if ~isfinite(p.va) || ~isfinite(p.vr) || p.va <= 0 || p.vr <= 0, f = 1e6; return, end
if ~isempty(lam), so.init = p211_spatcomp(p, lam); end
if isfield(so,'init') && any(~isfinite(so.init)), f = 1e6; return, end
m = p211_sim(p, so);
t = p211_compare(d, m);
r = t.rel(~isnan(t.rel));
f = sum(r.^2);
if isnan(m.thetaA) || ~isfinite(f), f = 1e6; end
end

function x = pack(p)
pi0 = p211_occ(p);
x = [log(p.va) log(p.vr) log(max(pi0,1e-9)) log(p.kap) log(p.krp) log(p.kpd) log(p.kds)];
end

function p = unpack(x)
p.va = exp(x(1));  p.vr = exp(x(2));
w  = exp(x(3:7) - max(x(3:7)));  pi0 = w/sum(w);
pi0 = max(pi0, 1e-6);  pi0 = pi0/sum(pi0);
p.kap = exp(x(8));  p.krp = exp(x(9));  p.kpd = exp(x(10));  p.kds = exp(x(11));
p.kpa = p.kap*pi0(1)/pi0(3);  p.kpr = p.krp*pi0(2)/pi0(3);
p.kdp = p.kpd*pi0(3)/pi0(4);  p.ksd = p.kds*pi0(4)/pi0(5);
p.pi  = pi0;
end
