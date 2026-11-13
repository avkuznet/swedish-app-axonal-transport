function R = run_product_model(outfile)
%RUN_PRODUCT_MODEL  A literature-constrained fluorescence model for Fig. 2C of Feole et al. (2024).
%
%   Transport of intact APP: the population-consistent chains of
%   run_refit12 (movies and space sample one open population), unchanged.
%   Fluorescence: intact APP plus ONE additional tag-bearing pool, chosen
%   and constrained by measurements outside Fig. 2 of Feole et al. (2024):
%
%   (a) Candidates and their reach sqrt(D/lambda) from the literature
%       retained surface APP     Moghekar 2011: 5.3% of cellular APP,
%                                t1/2 = 19 h; membrane D (Belapurkar 2023)
%       membrane C83/C99-tag     Belapurkar 2023 D (membrane APP), and
%                                Wang 2024: CTFs barely detectable without
%                                gamma-inhibition, i.e. short-lived
%       soluble AICD-tag         Wang 2024 (AICD-mScarlet formed from a
%                                C-terminal tag), Tam 2014 (diffuse
%                                cytosolic tag after Swedish APP cleavage)
%   (b) Source: soma (Niederst 2015, ~80% of axonally released Abeta
%       processed in the soma), so the pool decays as exp(-x/ell_F) with
%       ell_F = sqrt(D_F/lambda_F).
%   (c) Genotype: a C-terminal tag is released by both alpha (C83) and
%       beta (C99) cleavage.  In neurons C83 is similar for the variants
%       and C99 is 2.8-fold higher for APPswe; with C83:C99 = 3:1 for APPwt
%       the tag-releasing cleavage is r = (3 + 2.8)/(3 + 1) = 1.45 times
%       higher for APPswe (Wang 2024).  The product-to-intact odds at the
%       hillock are therefore tied, odds_swe = r * odds_wt, instead of fitted.
%
%   Models fitted to Fig. 2C with the tau operator (kappa always fitted):
%     M0   intact APP only                                   k = 1
%     M2   two carrier populations                           k = 3
%     M3   common product share phi, length ell_F            k = 3
%     M4   product share tied by r (r = 1.45; 1.16-1.9; 1, 2.8)   k = 3
%          (Eqs. (10)-(11) of the paper; r is theta there)
%     M3g  free shares phi_w, phi_s                          k = 4
%   Each is then used, unchanged, to predict the APPwt/APPwt control
%   (Fig. S4C), and its ell_F is translated into the product diffusivity
%   that the measured product half-lives would require.  MATLAB R2026a.

if nargin < 1, outfile = 'run_product_model_output.txt'; end
S0 = load('p211_spatial2.mat');  R2 = S0.R;
F = R2.F;  xg = R2.xg;  At = R2.At;  neff = R2.neff;
use = F.x >= 5;  xd = F.x(use);  [~, id] = ismember(xd, xg);
yd = [F.wt(use); F.swe(use)];  n = numel(yd);
L1 = load('p211_refit12_1.mat');  L2 = load('p211_refit12_2.mat');
Pc = {L1.R.p, L2.R.p};
G = p211_figS4c();
A4 = interp1(G.x, G.tau, xg, 'pchip');  A4 = A4/A4(1);
k4 = G.x >= 5;  [~, i4] = ismember(G.x(k4), xg);  obs4 = [G.gfp(k4); G.rfp(k4)];

fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('===========================================================\n');
pr(' P211  literature-constrained fluorescence model (population-consistent chains)\n');
pr(' %s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
pr('===========================================================\n');

% ---------------------------------------------------- (a) reach of each candidate pool
pr('\n(a) reach sqrt(D/lambda) of candidate pools (um)\n');
DmW = 0.041;  DmS = 0.007;                 % Belapurkar 2023, Table 2, medians (um^2/s)
l19 = log(2)/(19*3600);  l1 = log(2)/3600;  l25 = log(2)/(2.5*3600);
pr('   retained surface APP, t1/2 19 h:  WT %.0f, Swe %.0f  (abundance 5.3%% of cellular APP; Moghekar 2011)\n', ...
   sqrt(DmW/l19), sqrt(DmS/l19));
pr('   membrane CTF-tag, t1/2 1 h:       WT %.0f, Swe %.0f   (upper IQR D = 0.150: %.0f)\n', ...
   sqrt(DmW/l1), sqrt(DmS/l1), sqrt(0.150/l1));

% ---------------------------------------------------- (b) fits to Fig. 2C
rs = [1.45 1.16 1.9 1 2.8];
for lamh = [3.63 2.5]
    lam = log(2)/(lamh*3600);
    B = p211_spatial2('chains', Pc, lam, xg);  cW = B.cW;  cS = B.cS;
    pr('\n(b) Fig. 2C, tau operator, intact-APP half-life %.2f h: chains ell %.0f / %.0f um\n', lamh, B.ellW, B.ellS);
    mods = {'M0',[],'M2',[],'M3',[],'M4',1.45,'M4',1.16,'M4',1.9,'M4',1,'M4',2.8,'M3g',[]};
    fits = {};
    for j = 1:2:numel(mods)
        f = fitm(mods{j}, mods{j+1}, At, cW, cS, xg, id, yd);
        f.aic = neff*log(f.rss/n) + 2*f.k;  fits{end+1} = f; %#ok<AGROW>
    end
    amin = min(cellfun(@(f) f.aic, fits));
    for j = 1:numel(fits)
        f = fits{j};
        % prediction of the APPwt/APPwt control, every constant fixed
        sh = A4.*f.wtShape;  s2 = [sh(i4); sh(i4)];  p1 = f.kap*s2;
        kr = (s2'*obs4)/(s2'*s2);
        pr('   %-12s k=%d rms %.4f  areas %.1f/%.1f %.1f/%.1f  distal %.3f prox %.3f  dAIC %+5.1f | S4C rms %.4f (kappa refit %.4f) | %s\n', ...
           f.name, f.k, sqrt(f.rss/n), f.auc, f.auc(4)/f.auc(2), f.auc(3)/f.auc(1), f.aic - amin, ...
           sqrt(mean((p1 - obs4).^2)), sqrt(mean((kr*s2 - obs4).^2)), f.txt);
        if lamh == 3.63, R.fits{j} = f; end
    end
end
pr('   measured: areas 56.0/43.0 54.0/34.0, distal 0.791, prox 0.964; S4C: rms between its two APPwt channels %.3f\n', ...
   sqrt(mean((G.gfp(k4) - G.rfp(k4)).^2)));

% ---------------------------------------------------- (c) implied product lifetime
f = R.fits{4};                              % M4, r = 1.45
cW = p211_spatial2('chains', Pc, log(2)/(3.63*3600), xg);  cS = cW.cS;  cW = cW.cW;
pr('\n(c) M4 (r = 1.45): ell_F = %.0f um; product half-life implied for D_F = 1 / 3 / 10 / 30 um^2/s: %s h\n', ...
   f.lF, sprintf('%.2f ', log(2)*f.lF^2./[1 3 10 30]/3600));
pr('    D_F required for ell_F = %.0f um at product t1/2 = 0.5 / 1 / 2.5 / 6 h: %s um^2/s\n', ...
   f.lF, sprintf('%.1f ', f.lF^2*log(2)./([0.5 1 2.5 6]*3600)));
pr('    hillock product shares phi_w %.3f, phi_s %.3f; distal share of APPwt signal from the product %.2f\n', ...
   f.pw, f.ps, f.distShareW);
% product-to-intact ratio at the hillock and the product lifetime it implies
% if nearly all APP is eventually cleaved (Moghekar 2011: ~90% released as
% ectodomain within 9 h), so that product/intact = tau_F/tau_M
for tM = [2.5 3.63]
    tF = tM*f.pw/(1-f.pw);
    pr('    intact t1/2 %.2f h: product/intact %.2f -> product t1/2 %.1f h -> D_F = %.1f um^2/s\n', ...
       tM, f.pw/(1-f.pw), tF, log(2)*f.lF^2/(tF*3600));
end

% ---------------------------------------------------- (d) held-out test of APPswe
% Fit kappa, phi_w and ell_F to the APPwt/APPwt control alone (Fig. S4C,
% no Swedish protein, its own tau profile), then predict BOTH Fig. 2C
% profiles with phi_s from r; only kappa is refitted for Fig. 2C, which
% does not affect the APPswe/APPwt ratios.
pr('\n(d) product parameters from the APPwt/APPwt control only, APPswe predicted\n');
sig = @(z) 1./(1 + exp(-z));
shp = @(u) A4.*((1-sig(u(2)))*cW + sig(u(2))*exp(-xg/exp(u(3))));
obj4 = @(u) sum(([exp(u(1))*sel1(shp(u),i4); exp(u(1))*sel1(shp(u),i4)] - obs4).^2);
opt = optimset('TolX',1e-9,'TolFun',1e-12,'MaxFunEvals',4e4,'MaxIter',4e4,'Display','off');
best = inf;
for st = {[log(0.8) 0 log(300)], [log(0.8) 1 log(600)], [log(0.8) -1 log(100)], [log(0.8) 2 log(1500)], [log(0.8) -2 log(50)]}
    [u, v] = fminsearch(obj4, st{1}, opt);  if v < best, best = v; ub = u; end
end
pw4 = sig(ub(2));  lF4 = exp(ub(3));
pr('    S4C fit: rms %.4f  phi_w %.3f  ell_F %.0f um\n', sqrt(best/numel(obs4)), pw4, lF4);
% profile of the S4C objective over phi_w, and the APPswe prediction along it
pr('    phi_w  ell_F  S4C rms | predicted Fig. 2C: distal ratio  prox ratio  rms (kappa refit)   [r = 1.45]\n');
for pwf = [0.3 0.5 0.6 0.7 pw4 0.8 0.85 0.9]
    o2 = @(v) sum(([exp(v(1))*sel1(A4.*((1-pwf)*cW + pwf*exp(-xg/exp(v(2)))),i4); exp(v(1))*sel1(A4.*((1-pwf)*cW + pwf*exp(-xg/exp(v(2)))),i4)] - obs4).^2);
    bb = inf;
    for st = {[log(0.8) log(100)], [log(0.8) log(500)], [log(0.8) log(2000)]}
        [v, val] = fminsearch(o2, st{1}, opt);  if val < bb, bb = val; vb = v; end
    end
    lFi = exp(vb(2));  e = exp(-xg/lFi);  ow = pwf/(1-pwf);  psi = 1.45*ow/(1+1.45*ow);
    w = At.*((1-pwf)*cW + pwf*e);  s = At.*((1-psi)*cS + psi*e);
    yy = [w(id); s(id)];  kk = (yy'*yd)/(yy'*yy);
    iP = xg <= 100;  iD = xg >= 100;
    pr('    %.3f %6.0f  %.4f  | %.3f  %.3f  %.4f\n', pwf, lFi, sqrt(bb/numel(obs4)), ...
       trapz(xg(iD),s(iD))/trapz(xg(iD),w(iD)), trapz(xg(iP),s(iP))/trapz(xg(iP),w(iP)), sqrt(mean((kk*yy - yd).^2)));
end
pr('    measured Fig. 2C: distal 0.791 (replicates 0.76 +/- 0.14), prox 0.964\n');
fclose(fid);
save('p211_product_model.mat','R');
end

% ------------------------------------------------------------------------
function f = fitm(name, r, A, cW, cS, xg, id, yd)
sig = @(z) 1./(1 + exp(-z));
switch name
    case 'M0',  starts = {log(0.8)};
    case 'M2',  starts = {[log(0.8) 0 0], [log(0.8) -2 -2], [log(0.8) -2 1], [log(0.8) 1 -2]};
    case {'M3','M4'}, starts = {[log(0.8) 0 log(300)], [log(0.8) 1 log(600)], [log(0.8) -1 log(100)], [log(0.8) 2 log(1000)]};
    case 'M3g', starts = {[log(0.8) 0 0 log(300)], [log(0.8) 1 1 log(600)], [log(0.8) 1 2 log(500)], [log(0.8) -1 -1 log(100)]};
end
opt = optimset('TolX',1e-9,'TolFun',1e-12,'MaxFunEvals',4e4,'MaxIter',4e4,'Display','off');
pred = @(u) model(name, r, u, A, cW, cS, xg, sig);
best = inf;
for s = 1:numel(starts)
    [u, v] = fminsearch(@(u) sum((sel(pred(u), id) - yd).^2), starts{s}, opt);
    if v < best, best = v; ub = u; end
end
[out, fw, fs, info] = pred(ub); %#ok<ASGLU>
f = info;  f.k = numel(ub);  f.rss = best;  f.u = ub;
iP = xg <= 100;  iD = xg >= 100;
f.auc = [trapz(xg(iP),fw(iP)) trapz(xg(iD),fw(iD)) trapz(xg(iP),fs(iP)) trapz(xg(iD),fs(iD))];
if strcmp(name,'M4'), f.name = sprintf('M4 r=%.2f', r); else, f.name = name; end
end

function y = sel1(v, id)
y = v(id);
end

function y = sel(c, id)
y = [c{1}(id); c{2}(id)];
end

function [c, fw, fs, info] = model(name, r, u, A, cW, cS, xg, sig)
kap = exp(u(1));  info.kap = kap;  info.lF = NaN;  info.pw = 0;  info.ps = 0;
switch name
    case 'M0'
        w = cW;  s = cS;  info.txt = sprintf('kappa %.3f', kap);
    case 'M2'
        mw = sig(u(2));  ms = sig(u(3));
        w = (1-mw)*cW + mw*cS;  s = ms*cW + (1-ms)*cS;
        info.txt = sprintf('kappa %.3f m_w %.3f m_s %.3f', kap, mw, ms);
    case {'M3','M4','M3g'}
        if strcmp(name,'M3g')
            pw = sig(u(2));  ps = sig(u(3));  lF = exp(u(4));
        else
            pw = sig(u(2));  lF = exp(u(3));
            if strcmp(name,'M4'), ow = pw/(1-pw);  ps = r*ow/(1 + r*ow); else, ps = pw; end
        end
        e = exp(-xg/lF);
        w = (1-pw)*cW + pw*e;  s = (1-ps)*cS + ps*e;
        info.lF = lF;  info.pw = pw;  info.ps = ps;
        info.txt = sprintf('kappa %.3f phi_w %.3f phi_s %.3f ell_F %.0f um', kap, pw, ps, lF);
end
info.wtShape = w;
iD = xg >= 100;
info.distShareW = trapz(xg(iD), info.pw*exp(-xg(iD)/max(info.lF,1)))/max(trapz(xg(iD), w(iD)), eps);
fw = kap*A.*w;  fs = kap*A.*s;  c = {fw, fs};
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
