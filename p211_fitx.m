function R = p211_fitx(d, p0, opt)
%P211_FITX  Fit one condition under the corrected observation rules, with options.
%
%   R = p211_fitx(d, p0, opt)
%     d    one element of p211_data_pooled (targets may be perturbed)
%     p0   starting parameter struct
%     opt.mi       structure: 1 = A,R,P (6 constants); 2 = A,R,P,S (8);
%                  3 = five states (10; default)
%     opt.fitcv    true: also fit a lognormal run-speed CV (one more constant)
%     opt.cv0      starting CV (default 0.3)
%     opt.lam      if set, movie initial states from the open population
%     opt.nstart   perturbed restarts besides p0 (default 4)
%     opt.maxfun   simplex evaluations per pass (default 1500)
%     opt.N        trajectories per evaluation (default 1500)
%   Every candidate is re-scored on 20,000 trajectories (seed 99001);
%   R.p is the best, with R.p.speedcv set when opt.fitcv.  MATLAB R2026a.

if nargin < 3, opt = struct; end
def = struct('mi',3,'fitcv',false,'cv0',0.3,'lam',[],'nstart',4,'maxfun',1500,'N',1500);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opt,fn{i}), opt.(fn{i}) = def.(fn{i}); end, end

so = struct('N',opt.N,'vthr',d.vthr,'seed',20261104,'exhaustive',true);
x0 = pack(p0, opt);
rs = RandStream('mt19937ar','Seed',555);
X0 = x0;
if opt.fitcv
    % the constant-speed model is the limit CV -> 0; start next to it as well
    xs = x0;  xs(end) = log(0.05);  X0(end+1,:) = xs;
end
for j = 1:opt.nstart, X0(end+1,:) = x0 + 0.3*randn(rs,1,numel(x0)); end %#ok<AGROW>
o = optimset('Display','off','MaxFunEvals',opt.maxfun,'MaxIter',opt.maxfun,'TolX',1e-4,'TolFun',1e-5);
C = struct('p',{},'f',{});
for j = 1:size(X0,1)
    xh = fminsearch(@(x) resid(x, d, so, opt), X0(j,:), o);
    xh = fminsearch(@(x) resid(x, d, so, opt), xh, o);
    p = unpack(xh, opt);
    C(end+1) = struct('p',p,'f',score(p, d, opt)); %#ok<AGROW>
end
f0 = score(unpack(x0, opt), d, opt);
[fb, ib] = min([C.f]);
if f0 < fb, R.p = unpack(x0, opt); R.f = f0; else, R.p = C(ib).p; R.f = fb; end
R.cand = C;  R.opt = opt;  R.name = d.name;
[~, R.m, R.table] = score(R.p, d, opt);
r = R.table.rel(~isnan(R.table.rel));  R.rms = sqrt(mean(r.^2));  R.nt = numel(r);
end

% ------------------------------------------------------------------------
function [f, m, t] = score(p, d, opt)
o = struct('N',20000,'vthr',d.vthr,'seed',99001,'exhaustive',true);
if ~isempty(opt.lam), o.init = p211_spatcomp(p, opt.lam); end
if isfield(p,'speedcv'), o.speedcv = p.speedcv; end
m = p211_sim(p, o);  t = p211_compare(d, m);
f = sum(t.rel(~isnan(t.rel)).^2);
m = rmfield(m, 'X');        % the simulated trajectories are not stored with the fit
end

function f = resid(x, d, so, opt)
p = unpack(x, opt);
if ~isfinite(p.va) || ~isfinite(p.vr) || p.va <= 0 || p.vr <= 0, f = 1e6; return, end
if ~isempty(opt.lam)
    so.init = p211_spatcomp(p, opt.lam);
    if any(~isfinite(so.init)), f = 1e6; return, end
end
if isfield(p,'speedcv'), so.speedcv = p.speedcv; end
m = p211_sim(p, so);  t = p211_compare(d, m);
r = t.rel(~isnan(t.rel));  f = sum(r.^2);
if isnan(m.thetaA) || ~isfinite(f), f = 1e6; end
end

function x = pack(p, opt)
pi0 = max(p211_occ(p), 1e-6);
x = [log(p.va) log(p.vr) log(pi0) log(p.kap) log(p.krp) log(max(p.kpd,1e-6)) log(max(p.kds,1e-6))];
if opt.fitcv
    cv = opt.cv0;  if isfield(p,'speedcv'), cv = p.speedcv; end
    x(end+1) = log(cv);
end
end

function p = unpack(x, opt)
p.va = exp(x(1));  p.vr = exp(x(2));
w = exp(x(3:7) - max(x(3:7)));  pi0 = w/sum(w);
switch opt.mi
    case 1, pi0(4) = 0;  pi0(5) = 0;
    case 2, pi0(4) = 1e-4*sum(pi0);
end
pi0 = max(pi0, 1e-9);  pi0 = pi0/sum(pi0);
p.kap = exp(x(8));  p.krp = exp(x(9));  p.kpd = exp(x(10));  p.kds = exp(x(11));
if opt.mi == 1, p.kpd = 1e-9;  p.kds = 1e-9; end
if opt.mi == 2, p.kds = 1e3*p.kpd; end          % D never retains cargo
p.kpa = p.kap*pi0(1)/pi0(3);  p.kpr = p.krp*pi0(2)/pi0(3);
p.kdp = p.kpd*pi0(3)/max(pi0(4),1e-12);  p.ksd = p.kds*pi0(4)/max(pi0(5),1e-12);
p.pi = pi0;
if opt.fitcv, p.speedcv = min(exp(x(12)), 2); end
end
