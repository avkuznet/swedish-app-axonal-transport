function R = run_spatial2(outfile)
%RUN_SPATIAL2  Builds p211_spatial2.mat: digitised profiles, grid, tau profile.
%
%   ROLE IN THE SHORT PAPER.  run_short_all does not call this script.  It
%   is kept because it produced p211_spatial2.mat, from which the scripts
%   of the paper (run_product_model, run_coupled_product,
%   make_short_figures) read four fields only:
%       R.F     digitised Fig. 2C profiles (p211_fig2c)
%       R.xg    model grid, 0:0.5:200 um
%       R.At    tau profile on that grid, tau(x)/tau(0)
%       R.neff  effective sample size of the residuals
%   The candidate fits that it also stores and prints (M0 ... M3 below)
%   use the REFERENCE vectors of p211_fitted.mat, not the adopted fits,
%   and the paper quotes none of them.  The fits of the paper are in
%   run_product_model.
%
%   WHAT IT DOES.  Fig. 2C of Feole et al. (2024) gives the whole profile
%   of APPwt, of APPswe and of tau along the first 200 um, in the same
%   axons.  p211_fig2c holds it, digitised and checked against the
%   published areas to within 3%.  This routine asks what the profiles
%   demand of the model, using two facts about the experiment.
%
%   (1) Geometry.  Fixed-cell intensity is concentration times the amount
%       of axon in the pixel, and an axon narrows away from the hillock.
%       Tau was stained in the same axons as an internal control; a
%       microtubule-bound protein at near-uniform linear density reports
%       that taper directly.  The observation operator for the spatial
%       assay is therefore
%
%           f_g(x) = kappa * A(x) * c_g(x)/c_g(0),   A(x) = tau(x)/tau(0),
%
%       with kappa the hillock-to-soma ratio, a property of the
%       normalisation and the same for both variants.  A = 1 omits the
%       tau correction.
%
%   (2) Co-expression.  The movies that fixed the chain came from neurons
%       expressing one variant only; Fig. 2 co-expressed APPwt-GFP and
%       APPswe-tRFP in the same neurons.  APP dimerises and the two
%       variants share the secretory pathway, so in a co-expressing axon
%       neither protein need travel only on carriers with its own
%       single-variant kinetics.  The simplest description is two carrier
%       populations, one moving with the fitted APPwt chain (W) and one
%       with the fitted APPswe chain (S), each protein appearing in both:
%
%           c_wt  = (1 - m_w) cW + m_w cS
%           c_swe =  m_s cW + (1 - m_s) cS
%
%       cW and cS are the five-state profiles of p211_profile, normalised
%       to one at the hillock, with NO new kinetic constant; m_w and m_s
%       are the fractions of each protein's hillock signal found on the
%       other population.  m_w = m_s = 0 is transport alone.
%
%   Candidates, each with kappa fitted (one nuisance constant):
%     M0     transport alone                                    0 extra
%     M0lam  M0 with the turnover rate free                     1 extra
%     M1     shared carriers, m_w = m_s                         1 extra
%     M2     shared carriers, m_w and m_s free                  2 extra
%     M3     genotype-independent second route, weight phi and
%            its own decay length ell_B                         2 extra
%   each with A = 1 and with the tau operator.  A further variant, M2L,
%   replaces cS by exp(-x/L) with L free, which asks whether the data
%   choose the decay length the APPswe chain predicts.
%
%   Residuals of a digitised smooth curve are serially correlated, so the
%   information criteria are also reported with an effective sample size
%   n_eff = n(1 - rho)/(1 + rho), rho the lag-1 autocorrelation.

if nargin < 1, outfile = 'run_spatial2_output.txt'; end

S = load('p211_fitted.mat');  R0 = S.R;
P = R0.P;
LAM = log(2)/(3.63*3600);
F = p211_fig2c();

xg = (0:0.5:200)';                         % model grid
use = F.x >= 5;                             % fitted data points
xd = F.x(use);
[~, id] = ismember(xd, xg);
yd = [F.wt(use); F.swe(use)];
n = numel(yd);

A1 = ones(size(xg));
At = interp1(F.x, F.tau, xg, 'pchip');  At = At/At(1);

base = p211_spatial2('chains', P, LAM, xg);
cW = base.cW;  cS = base.cS;

fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('===========================================================\n');
pr(' P211  spatial test against the measured profiles (Fig. 2C)\n');
pr(' %s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
pr('===========================================================\n');
pr('\n data: %s\n', F.src);
pr(' %d points (x = 5..200 um, both genotypes); lam = %.3e 1/s\n', n, LAM);
pr(' chain decay lengths: APPwt %.0f um, APPswe %.1f um\n', base.ellW, base.ellS);

% ------------------------------------------------------------- models
M = {'M0','M0lam','M1','M2','M3'};
lab = {'transport alone','free turnover','shared carriers, m_w=m_s', ...
       'shared carriers, m_w, m_s','second route (phi, ell_B)'};
ops = {A1, At};  opn = {'A = 1','tau operator'};
fit = cell(numel(M), 2);
for io = 1:2
    for im = 1:numel(M)
        fit{im,io} = p211_spatial2('fit', M{im}, ops{io}, cW, cS, P, xg, id, yd);
    end
end

% effective sample size from the residuals of the best two-parameter fit
rho = lag1(fit{4,2}.res, numel(xd));
neff = n*(1 - rho)/(1 + rho);
pr(' lag-1 residual autocorrelation %.2f  ->  n_eff = %.0f of %d\n', rho, neff, n);

for io = 1:2
    pr('\n--- %s -------------------------------------------------------\n', opn{io});
    pr('  %-6s %-28s %2s %8s %13s %13s %8s %9s %9s\n', 'model','', 'k', ...
       'rms', 'AUC wt', 'AUC swe', 'ratio', 'dAIC', 'dAIC_eff');
    aic = zeros(1,numel(M));  aice = aic;
    for im = 1:numel(M)
        f = fit{im,io};
        aic(im)  = n*log(f.rss/n) + 2*f.k;
        aice(im) = neff*log(f.rss/n) + 2*f.k;
    end
    for im = 1:numel(M)
        f = fit{im,io};
        pr('  %-6s %-28s %2d %8.4f %6.1f /%5.1f %6.1f /%5.1f %8.3f %9.1f %9.1f\n', ...
           M{im}, lab{im}, f.k, sqrt(f.rss/n), f.auc(1), f.auc(2), f.auc(3), ...
           f.auc(4), f.ratio, aic(im)-min(aic), aice(im)-min(aice));
    end
    pr('  %-6s %-28s %2s %8s %6.1f /%5.1f %6.1f /%5.1f %8.3f\n', 'data', ...
       'published Fig. 2D', '', '', 56.0, 43.0, 54.0, 34.0, 34/43);
    for im = 1:numel(M)
        pr('     %-6s %s\n', M{im}, fit{im,io}.ptxt);
    end
    fitAIC{io} = aic; fitAICe{io} = aice;                         %#ok<AGROW>
end

% ----------------------------------------- uncertainty in the weights
pr('\n--- moving-block bootstrap of the shared-carrier weights ----------\n');
pr('  (residual blocks of 20 um resampled within each genotype; 400 refits)\n');
rng(21101);
NB = 400;  blk = 4;
boot = struct();
for io = 1:2
    f2 = fit{4,io};  f1 = fit{3,io};
    th2 = zeros(NB,3);  th1 = zeros(NB,2);  rb = zeros(NB,1);
    for b = 1:NB
        ystar = f2.yhat + blockres(f2.res, numel(xd), blk);
        g2 = p211_spatial2('fit', 'M2', ops{io}, cW, cS, P, xg, id, ystar, f2.u);
        th2(b,:) = g2.theta;  rb(b) = g2.ratio;
        ystar = f1.yhat + blockres(f1.res, numel(xd), blk);
        g1 = p211_spatial2('fit', 'M1', ops{io}, cW, cS, P, xg, id, ystar, f1.u);
        th1(b,:) = g1.theta;
    end
    q = @(v) pctl(v,[2.5 97.5]);
    pr('  %s\n', opn{io});
    pr('    M1  m        %.3f   [%.3f, %.3f]\n', f1.theta(2), q(th1(:,2)));
    pr('    M2  m_w      %.3f   [%.3f, %.3f]\n', f2.theta(2), q(th2(:,2)));
    pr('    M2  m_s      %.3f   [%.3f, %.3f]\n', f2.theta(3), q(th2(:,3)));
    pr('    M2  distal ratio %.3f [%.3f, %.3f]   measured %.3f\n', ...
       f2.ratio, q(rb), 34/43);
    boot(io).th1 = th1;  boot(io).th2 = th2;  boot(io).ratio = rb;
end

% --------------------- does the data choose the APPswe chain's length?
pr('\n--- the proximal length scale: chosen by the data, or by the chain? ---\n');
Lg = (10:2:200)';
for io = 1:2
    rssL = zeros(size(Lg));
    for j = 1:numel(Lg)
        eS = exp(-xg/Lg(j));
        g = p211_spatial2('fit', 'M2', ops{io}, cW, eS, P, xg, id, yd);
        rssL(j) = g.rss;
    end
    [rmin, jb] = min(rssL);
    % 95% profile interval, chi-square(1) on the effective sample
    crit = rmin*exp(3.84/neff);
    inL = Lg(rssL <= crit);
    pr('  %-13s best L = %3.0f um, 95%% profile interval [%3.0f, %3.0f] um\n', ...
       opn{io}, Lg(jb), min(inL), max(inL));
    LL(io) = struct('L',Lg,'rss',rssL,'best',Lg(jb),'lo',min(inL),'hi',max(inL)); %#ok<AGROW>
end
pr('  predicted by the APPswe chain: %.1f um\n', base.ellS);

% --------------------------------------- weights as fractions of flux
pr('\n--- what the weights mean as fractions of exported flux -----------\n');
f2 = fit{4,2};
mw = f2.theta(2);  ms = f2.theta(3);
qw = (mw/base.kS)/((1-mw)/base.kW + mw/base.kS);
qs = (ms/base.kW)/(ms/base.kW + (1-ms)/base.kS);
pr('  signal per unit flux at the hillock: W %.0f s/um,  S %.0f s/um  (ratio %.0f)\n', ...
   base.kW, base.kS, base.kS/base.kW);
pr('  APPwt : %.0f%% of hillock signal on S carriers = %.2f%% of its exported flux\n', ...
   100*mw, 100*qw);
pr('  APPswe: %.0f%% of hillock signal on W carriers = %.1f%% of its exported flux\n', ...
   100*ms, 100*qs);

% ------------------------------------------------------------- save
R = struct('F',F,'xg',xg,'At',At,'cW',cW,'cS',cS,'base',base, ...
           'models',{M},'labels',{lab},'fit',{fit},'opn',{opn}, ...
           'rho',rho,'neff',neff,'boot',boot,'Lprof',LL, ...
           'AIC',{fitAIC},'AICeff',{fitAICe},'flux',[qw qs]);
save('p211_spatial2.mat','R');
pr('\n saved p211_spatial2.mat\n');
fclose(fid);
end

% ========================================================================
function q = pctl(v, p)
%PCTL  Percentiles by linear interpolation (base MATLAB, no toolbox).
v = sort(v(:));  nv = numel(v);
q = interp1((0:nv-1)'/(nv-1)*100, v, p);
end

function rho = lag1(res, m)
r1 = res(1:m);  r2 = res(m+1:end);
c = @(r) sum(r(1:end-1).*r(2:end))/sum(r.^2);
rho = max(0, min(0.95, 0.5*(c(r1 - mean(r1)) + c(r2 - mean(r2)))));
end

function e = blockres(res, m, L)
%BLOCKRES  Moving-block resample of residuals, separately per genotype.
e = zeros(size(res));
for g = 0:1
    r = res(g*m + (1:m));
    out = zeros(m,1);  pos = 0;
    while pos < m
        st = randi(m - L + 1);
        take = min(L, m - pos);
        out(pos + (1:take)) = r(st + (0:take-1));
        pos = pos + take;
    end
    e(g*m + (1:m)) = out;
end
end

function fprintf2(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
