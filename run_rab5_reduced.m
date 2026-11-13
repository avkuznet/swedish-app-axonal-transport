function run_rab5_reduced(k, mi)
%RUN_RAB5_REDUCED  Reduced chains for Rab5, where the immobile states are least identified.
%
%   run_rab5_reduced(k, mi)   k = 3 (Rab5) or 4 (Rab5+APPswe)
%                             mi = 1: A, R, P            (6 constants)
%                             mi = 2: A, R, P, S         (8 constants)
%   The Rab5 summaries give eight distinct targets and at most seven
%   independent constraints (pooled track length and pause frequency, no
%   whole-track classes), so the ten-constant chain is under-determined by
%   three.  Each reduced chain is fitted with the corrected targets
%   (p211_data_pooled, exhaustive classes) from the 5-state refit
%   (run_refit12) projected onto the reduced structure and from random
%   starts, and scored on the common 20,000-trajectory sample.
%   Same structures as p211_fitx (mi = 1, 2).  Output p211_rab5red_<k>_<mi>.mat.
%   MATLAB R2026a.

Dp = p211_data_pooled();  d = Dp(k);
L = load(sprintf('p211_refit12_%d.mat', k));  p5 = L.R.p;
sp = struct('N',1500,'vthr',d.vthr,'seed',20261104,'exhaustive',true);
ss = struct('N',20000,'vthr',d.vthr,'seed',99001,'exhaustive',true);
rs = RandStream('mt19937ar','Seed',4040 + 10*k + mi);
X0 = pack(p5);
nr = 5;
for j = 1:nr, X0(end+1,:) = X0(1,:) + 0.5*randn(rs,1,11); end %#ok<AGROW>
o = optimset('Display','off','MaxFunEvals',1500,'MaxIter',1500,'TolX',1e-4,'TolFun',1e-5);
best = inf;  C = struct('p',{},'f',{});
for j = 1:size(X0,1)
    xh = fminsearch(@(x) score(unpack(x,mi), d, sp), X0(j,:), o);
    xh = fminsearch(@(x) score(unpack(x,mi), d, sp), xh, o);
    ph = unpack(xh, mi);  fh = score(ph, d, ss);
    C(end+1) = struct('p',ph,'f',fh); %#ok<AGROW>
    fprintf('  %s mi=%d start %d: objective %.4f\n', d.name, mi, j, fh);
    if fh < best, best = fh; bp = ph; end
end
m = p211_sim(bp, ss);  t = p211_compare(d, m);  r = t.rel(~isnan(t.rel));
R = struct('name',d.name,'mi',mi,'p',bp,'f',best,'rms',sqrt(mean(r.^2)),'table',t,'cand',C);
save(sprintf('p211_rab5red_%d_%d.mat', k, mi), 'R');
ns = [3 4];
fprintf('%s %d-state: objective %.4f, rms %.1f%%\n', d.name, ns(mi), best, 100*R.rms);
end

function p = unpack(x, mi)
p.va = exp(x(1)); p.vr = exp(x(2));
w = exp(x(3:7) - max(x(3:7)));  pi0 = w/sum(w);
switch mi
    case 1, pi0(4) = 0; pi0(5) = 0;
    case 2, pi0(4) = 1e-4*sum(pi0);
end
pi0 = max(pi0,1e-9); pi0 = pi0/sum(pi0);
p.kap = exp(x(8)); p.krp = exp(x(9)); p.kpd = exp(x(10)); p.kds = exp(x(11));
if mi == 1
    p.kpd = 1e-9; p.kds = 1e-9;
else
    p.kds = 1e3*p.kpd;                       % D never retains cargo
end
p.kpa = p.kap*pi0(1)/pi0(3);  p.kpr = p.krp*pi0(2)/pi0(3);
p.kdp = p.kpd*pi0(3)/max(pi0(4),1e-12);  p.ksd = p.kds*pi0(4)/max(pi0(5),1e-12);
p.pi = pi0;
end

function x = pack(p)
pi0 = max(p211_occ(p), 1e-6);
x = [log(p.va) log(p.vr) log(pi0) log(p.kap) log(p.krp) log(max(p.kpd,1e-6)) log(max(p.kds,1e-6))];
end

function f = score(p, d, so)
if ~isfinite(p.va) || ~isfinite(p.vr) || p.va <= 0 || p.vr <= 0, f = 1e6; return, end
m = p211_sim(p, so);  t = p211_compare(d, m);  r = t.rel(~isnan(t.rel));
if isempty(r) || isnan(m.thetaA), f = 1e6; return, end
f = sum(r.^2);  if ~isfinite(f), f = 1e6; end
end
