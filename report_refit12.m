function Q = report_refit12(outfile)
%REPORT_REFIT12  Summarise the corrected refits of run_refit12.
%
%   For every condition: objective and rms of the reference vector and of
%   the refit under the same corrected rules (pooled targets, exhaustive
%   classes; APP also sampled from the open population at t1/2 = 3.63 h),
%   the refitted constants, analytic drift, and the residual table.
%   Then: retrograde-entry and event-flux ratios for Rab5 and lysosomes
%   over all candidates within 1.5x of the best objective; and the spatial
%   predictions of the population-consistent APP chains.  MATLAB R2026a.

if nargin < 1, outfile = 'report_refit12_output.txt'; end
S = load('p211_fitted.mat');  P0 = S.R.P;  Dp = p211_data_pooled();
lam = log(2)/(3.63*3600);
fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('P211  corrected refits (pooled targets, exhaustive classes; APP: open-population sampling)\n');
pr('%s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
R = cell(1,6);
for k = 1:6
    L = load(sprintf('p211_refit12_%d.mat', k));  R{k} = L.R;
    % the starting (reference) vector is itself a candidate under the same rules
    R{k}.cand(end+1) = struct('p',P0{k},'f',NaN,'fchk',R{k}.f0,'m',R{k}.m0);
    if R{k}.f0 < R{k}.fchk
        R{k}.p = P0{k};  R{k}.fchk = R{k}.f0;  R{k}.m = R{k}.m0;
        R{k}.table = p211_compare(Dp(k), R{k}.m);
        rr = R{k}.table.rel(~isnan(R{k}.table.rel));  R{k}.rms = sqrt(mean(rr.^2));
        R{k}.kept = 'reference vector kept';
    else
        R{k}.kept = 'refit';
    end
    r = R{k};  p = r.p;
    pi0 = p211_occ(p);
    pr('\n%-12s  reference vector %.4f -> best %.4f (%s)   rms %.1f%%\n', r.name, r.f0, r.fchk, r.kept, 100*r.rms);
    pr('   va %.3f vr %.3f | kap %.3g kpa %.3g krp %.3g kpr %.3g kpd %.3g kdp %.3g kds %.3g ksd %.3g\n', ...
       p.va, p.vr, p.kap, p.kpa, p.krp, p.kpr, p.kpd, p.kdp, p.kds, p.ksd);
    pr('   occupancy %s  v_eff (analytic) %+.4f um/s\n', sprintf('%.3f ', pi0), p.va*pi0(1) - p.vr*pi0(2));
    if k <= 2
        pr('   open-population composition %s\n', sprintf('%.3f ', p211_spatcomp(p, lam)));
    end
    t = r.table;
    for j = 1:numel(t.name)
        if ~isnan(t.obs(j)), pr('     %-22s meas %7.3f  model %7.3f  (%+.0f%%)\n', t.name{j}, t.obs(j), t.mod(j), 100*t.rel(j)); end
    end
end

% ---------------- motor-entry contrasts
pr('\nRetrograde entry, APPswe condition / control (best refits; range over candidates within 1.5x of each best)\n');
pairs = {[3 4], [5 6]};  lab = {'Rab5','lysosomes'};
for i = 1:2
    a = R{pairs{i}(1)};  b = R{pairs{i}(2)};
    ga = a.cand([a.cand.fchk] <= 1.5*a.fchk);  gb = b.cand([b.cand.fchk] <= 1.5*b.fchk);
    rr = [];  ff = [];
    for u = 1:numel(ga)
        for v = 1:numel(gb)
            pa = ga(u).p;  pb = gb(v).p;  oa = p211_occ(pa);  ob = p211_occ(pb);
            rr(end+1) = pb.kpr/pa.kpr;                    %#ok<AGROW>
            ff(end+1) = (pb.kpr*ob(3))/(pa.kpr*oa(3));    %#ok<AGROW>
        end
    end
    oa = p211_occ(a.p);  ob = p211_occ(b.p);
    pr('   %-9s kpr ratio %.2f (range %.2f-%.2f over %d pairs) | event flux kpr*piP ratio %.2f (range %.2f-%.2f)\n', ...
       lab{i}, b.p.kpr/a.p.kpr, min(rr), max(rr), numel(rr), (b.p.kpr*ob(3))/(a.p.kpr*oa(3)), min(ff), max(ff));
    pr('             reference vectors: kpr ratio %.2f, flux ratio %.2f\n', P0{pairs{i}(2)}.kpr/P0{pairs{i}(1)}.kpr, ...
       (P0{pairs{i}(2)}.kpr*p211_occ(P0{pairs{i}(2)})*[0;0;1;0;0])/(P0{pairs{i}(1)}.kpr*p211_occ(P0{pairs{i}(1)})*[0;0;1;0;0]));
end
pr('   APP kpr ratio: refit %.2f (reference %.2f)\n', R{2}.p.kpr/R{1}.p.kpr, P0{2}.kpr/P0{1}.kpr);

% ---------------- spatial predictions of the population-consistent APP chains
pr('\nSpatial predictions at t1/2 = 3.63 h\n');
pr('   measured: areas APPwt 56.0/43.0, APPswe 54.0/34.0; lengths 379, 216 um; ratios 0.964, 0.791\n');
lab2 = {'reference vectors','population-consistent refit'};
for j = 1:2
    if j == 1, pw = P0{1}; ps = P0{2}; else, pw = R{1}.p; ps = R{2}.p; end
    w = p211_profile(pw);  s = p211_profile(ps);
    fw = p211_profile_finite(pw, struct('Lax',Inf));  fs = p211_profile_finite(ps, struct('Lax',Inf));
    pr('   %-28s ell %6.0f / %5.0f um | areas %.1f/%.1f, %.1f/%.1f | ratios %.3f %.3f | J/C %.4f / %.4f um/s (ratio %.3f)\n', ...
       lab2{j}, w.ell, s.ell, w.aucProx, w.aucDist, s.aucProx, s.aucDist, s.aucProx/w.aucProx, s.aucDist/w.aucDist, ...
       lam*w.ell, lam*s.ell, s.ell/w.ell);
    Q.ell(j,:) = [w.ell s.ell];  Q.ratio(j,:) = [s.aucProx/w.aucProx s.aucDist/w.aucDist];
end
% half-life sensitivity of the consistent chains (the consistency itself depends on lam)
pr('   population-consistent chains refitted at 3.63 h, re-evaluated at other half-lives (composition re-derived):\n');
for th = [0.5 1 2 3.63 8 24]
    l2 = log(2)/(th*3600);
    w = p211_profile(R{1}.p, struct('lam',l2));  s = p211_profile(R{2}.p, struct('lam',l2));
    pr('     t1/2 %5.2f h  ell %6.0f / %5.0f um  distal ratio %.3f\n', th, w.ell, s.ell, s.aucDist/w.aucDist);
end
Q.R = R;
fclose(fid);
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
