function Q = report_rev1(outfile)
%REPORT_REV1  Model selection, run-to-run speed variability and perturbation-refit ranges.
%
%   Reads p211_red_*, p211_rab5red_*, p211_refit12_*, p211_cv_* and
%   p211_boot_* and writes report_rev1_output.txt.  MATLAB R2026a.
%
%   The p211_boot_* files hold the perturbation refits of Sec. 2.3 (targets
%   multiplied by exp(0.1 z), model refitted).  With ten refits per
%   condition the 2.5th and 97.5th percentiles printed below equal the
%   smallest and largest values, which is what the paper reports.
%
%   Model selection.  The objective F = sum(rel^2) is turned into a chi-square
%   with an assumed relative error sigma per target, chi2 = F/sigma^2, so that
%   AIC = chi2 + 2k and BIC = chi2 + k ln n, with k constants and n the number
%   of independent targets (APP 11, lysosomes 10, Rab5 7).  sigma = 0.10 is the
%   error assumed throughout; 0.15 matches the experiment-to-experiment spread
%   of the APPwt track classes in Feole et al. (2024).

if nargin < 1, outfile = 'report_rev1_output.txt'; end
fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('P211 short paper: model selection, run-speed variability, perturbation refits   %s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
names = {'APPwt','APPswe','Rab5','Rab5+APPswe','Lyso+APPwt','Lyso+APPswe'};
nind = [11 11 7 7 10 10];  kk = [6 8 10];
lam = log(2)/(3.63*3600);
S0 = load('p211_fitted.mat');
P5 = cell(1,6);  F5 = zeros(1,6);
for k = 1:6
    L = load(sprintf('p211_refit12_%d.mat',k));
    if k == 6, P5{k} = S0.R.P{6}; F5(k) = L.R.f0; else, P5{k} = L.R.p; F5(k) = L.R.fchk; end
end

% ------------------------------------------------------------ model selection
F = nan(6,3);  PR = cell(6,3);
for k = 1:6
    for mi = 1:2
        if k == 3 || k == 4, fn = sprintf('p211_rab5red_%d_%d.mat', k, mi); else, fn = sprintf('p211_red_%d_%d.mat', k, mi); end
        if exist(fn,'file'), L = load(fn); F(k,mi) = L.R.f; PR{k,mi} = L.R.p; end
    end
    F(k,3) = F5(k);  PR{k,3} = P5{k};
end
Q.F = F;
for sg = [0.10 0.15]
    pr('\nModel selection, sigma = %.2f   (objective | dAIC | dBIC relative to the best structure; 3-, 4-, 5-state)\n', sg);
    for k = 1:6
        chi = F(k,:)/sg^2;  aic = chi + 2*kk;  bic = chi + kk*log(nind(k));
        pr('  %-12s n=%2d  F %s | dAIC %s | dBIC %s\n', names{k}, nind(k), sprintf('%7.4f ', F(k,:)), ...
           sprintf('%+7.1f ', aic - min(aic)), sprintf('%+7.1f ', bic - min(bic)));
    end
end
pr('\nkpr ratio (APPswe / control) and event-flux ratio by structure\n');
for mi = 1:3
    for pair = {[3 4],[5 6]}
        a = pair{1}(1);  b = pair{1}(2);
        if isempty(PR{a,mi}) || isempty(PR{b,mi}), continue, end
        oa = p211_occ(PR{a,mi});  ob = p211_occ(PR{b,mi});
        pr('  %d-state  %-9s kpr ratio %.2f  event flux %.2f\n', 2+mi, names{a}, PR{b,mi}.kpr/PR{a,mi}.kpr, ...
           (PR{b,mi}.kpr*ob(3))/(PR{a,mi}.kpr*oa(3)));
    end
end
for mi = 1:3
    if isempty(PR{1,mi}) || isempty(PR{2,mi}), continue, end
    o1 = p211_occ(PR{1,mi});  o2 = p211_occ(PR{2,mi});
    pr('  %d-state  APP drift %+.4f / %+.4f\n', 2+mi, PR{1,mi}.va*o1(1)-PR{1,mi}.vr*o1(2), PR{2,mi}.va*o2(1)-PR{2,mi}.vr*o2(2));
end
% attenuation lengths with the preferred structure of each APP condition
% (five states for APPwt, four for APPswe; Appendix)
if ~isempty(PR{1,3}) && ~isempty(PR{2,2})
    w = p211_profile(PR{1,3});  s = p211_profile(PR{2,2});
    pr('  preferred structures (APPwt 5-state, APPswe 4-state): attenuation lengths %.0f / %.1f um\n', w.ell, s.ell);
end

% ------------------------------------------------------------ stochastic run speeds
Dp = p211_data_pooled();
pr('\nFive states with a fitted lognormal run-speed CV (one more constant)\n');
PC = cell(1,6);  nab = 0;
for k = 1:6
    fn = sprintf('p211_cv_%d.mat', k);  if ~exist(fn,'file'), continue, end
    L = load(fn);  p = L.R.p;  PC{k} = p;  o = p211_occ(p);  t = L.R.table;
    jA = find(strcmp(t.name,'seg speed antero'));  jR = find(strcmp(t.name,'seg speed retro'));
    chi = [F5(k) L.R.f]/0.1^2;  aic = chi + 2*[10 11];
    pr('  %-12s CV %.2f  objective %.4f (constant speed %.4f, dAIC %+.1f)  rms %.1f%%  va %.3f vr %.3f  segmental speed meas/model %.2f/%.2f, %.2f/%.2f  drift %+.4f\n', ...
       names{k}, p.speedcv, L.R.f, F5(k), aic(2)-aic(1), 100*L.R.rms, p.va, p.vr, t.obs(jA), t.mod(jA), t.obs(jR), t.mod(jR), p.va*o(1)-p.vr*o(2));
    nab = nab + (abs(t.rel(jA)) < 0.05) + (abs(t.rel(jR)) < 0.05);
end
pr('  segmental-speed medians reproduced within 5%%: %d of 12\n', nab);
if ~isempty(PC{1}) && ~isempty(PC{2})
    w = p211_profile(PC{1});  s = p211_profile(PC{2});
    pr('  APP with CV: drift ratio %.0f; attenuation lengths %.0f / %.0f um; distal ratio %.3f (mean-speed spatial model)\n', ...
       (PC{1}.va*p211_occ(PC{1})*[1;0;0;0;0] - PC{1}.vr*p211_occ(PC{1})*[0;1;0;0;0]) / ...
       (PC{2}.va*p211_occ(PC{2})*[1;0;0;0;0] - PC{2}.vr*p211_occ(PC{2})*[0;1;0;0;0]), w.ell, s.ell, s.aucDist/w.aucDist);
end
for pair = {[3 4],[5 6]}
    a = pair{1}(1);  b = pair{1}(2);
    if isempty(PC{a}) || isempty(PC{b}), continue, end
    oa = p211_occ(PC{a});  ob = p211_occ(PC{b});
    pr('  with CV: %-9s kpr ratio %.2f  event flux %.2f\n', names{a}, PC{b}.kpr/PC{a}.kpr, (PC{b}.kpr*ob(3))/(PC{a}.kpr*oa(3)));
end

% ------------------------------------------------------------ perturbation refits
pr('\nPerturbation refits: targets x exp(0.1 z), compositions renormalized, five-state refit\n');
pr('  (median [2.5th, 97.5th percentile]; with ten refits the bracket is the smallest and largest value)\n');
B = cell(1,6);
for k = [1 2 5 6]
    fs = dir(sprintf('p211_boot_%d_*.mat', k));  B{k} = {};
    for i = 1:numel(fs), L = load(fs(i).name); B{k}{end+1} = L.R.p; end
    pr('  %-12s %d refits\n', names{k}, numel(B{k}));
end
% run speeds of APPwt over its perturbation refits (Sec. 3.2)
if ~isempty(B{1})
    vaw = cellfun(@(p) p.va, B{1});  vrw = cellfun(@(p) p.vr, B{1});
    pr('  APPwt run speeds: vr %.2f-%.2f, va %.2f-%.2f um/s; vr > va in %d of %d refits\n', ...
       min(vrw), max(vrw), min(vaw), max(vaw), nnz(vrw > vaw), numel(vaw));
end
qs = @(v) sprintf('%.3g [%.3g, %.3g]', median(v), prctile(v,2.5), prctile(v,97.5));
der = @(p) deal(p.va*p211_occ(p)*[1;0;0;0;0] - p.vr*p211_occ(p)*[0;1;0;0;0], p211_occ(p)*[1;1;0;0;0], p.va/p.kap, p.vr/p.krp);
nb = min(numel(B{1}), numel(B{2}));
if nb > 0
    V = zeros(nb, 12);
    for b = 1:nb
        [v1, m1, ra1, rr1] = der(B{1}{b});  [v2, m2, ra2, rr2] = der(B{2}{b});
        o1 = p211_occ(B{1}{b});  o2 = p211_occ(B{2}{b});
        w = p211_profile(B{1}{b});  s = p211_profile(B{2}{b});
        V(b,:) = [v1 v2 v1/v2 m1/m2 ra1/ra2 rr2/rr1 B{2}{b}.va/B{1}{b}.va (o1(1)/o1(2)) (o2(1)/o2(2)) w.ell s.ell s.aucDist/w.aucDist];
    end
    lab = {'drift APPwt (um/s)','drift APPswe (um/s)','drift ratio wt/swe','motile occupancy ratio wt/swe', ...
           'anterograde run-scale ratio wt/swe','retrograde run-scale ratio swe/wt','va ratio swe/wt','piA/piR APPwt','piA/piR APPswe', ...
           'attenuation length APPwt (um)','attenuation length APPswe (um)','transport-only distal ratio'};
    for j = 1:numel(lab), pr('  %-36s %s\n', lab{j}, qs(V(:,j))); end
    pr('  drift ratio > 10 in %d of %d refits; APPswe drift between -0.01 and +0.01 um/s in %d\n', nnz(V(:,3) > 10 | V(:,2) <= 0), nb, nnz(abs(V(:,2)) < 0.01));
    Q.Vapp = V;
end
nl = min(numel(B{5}), numel(B{6}));
if nl > 0
    W = zeros(nl, 2);  RL = zeros(nl, 1);
    for b = 1:nl
        pa = B{5}{b};  pb = B{6}{b};  oa = p211_occ(pa);  ob = p211_occ(pb);
        W(b,:) = [pb.kpr/pa.kpr, (pb.kpr*ob(3))/(pa.kpr*oa(3))];
        RL(b)  = (pb.vr/pb.krp)/(pa.vr/pa.krp);       % mean retrograde run length, APPswe/APPwt
    end
    pr('  lysosome kpr ratio %s; event-flux ratio %s; both > 1 in %d of %d\n', qs(W(:,1)), qs(W(:,2)), nnz(all(W > 1, 2)), nl);
    pr('  lysosome retrograde run-length ratio (with APPswe / with APPwt) %s\n', qs(RL));
    Q.Wlys = W;  Q.RLlys = RL;
end
fclose(fid);
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
