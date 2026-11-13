function R = run_coupled_product(outfile)
%RUN_COUPLED_PRODUCT  Transport of intact APP explicitly coupled to a diffusing tag fragment.
%
%   Intact APP of genotype g: the population-consistent five-state chains
%   (run_refit12), total density C_g(x) = C_g0 exp(-x/ell_g), lost at the
%   rate lam of Eqs. (3)-(7).  Tag-bearing soluble fragment F_g:
%
%       D_F F'' - lam_F F + y lam C_g(x) = 0,   F_g(0) = r_g rho C_g0,  F bounded,
%
%   rho = lam/lam_F (fragment lifetime / intact lifetime), ell_F =
%   sqrt(D_F/lam_F).  y is the fraction of the intact-APP loss in the axon
%   that leaves a tag fragment in the axon: y = 0 is the soma-only product
%   of Eqs. (10)-(11) (M4 of run_product_model); Niederst et al. 2015 put ~80% of the processing
%   that supplies axonal products in the soma, i.e. y ~ 0.2; y = 1 is
%   complete local conversion.  r_wt = 1, r_swe = 1.45 (Wang et al. 2024).
%   Solution:
%       F = a exp(-x/ell_g) + (r_g rho C_g0 - a) exp(-x/ell_F),
%       a = y rho C_g0 / (1 - ell_F^2/ell_g^2).
%   Signal f_g = kappa A(x) [C_g + F_g]/[C_g0 + F_g0]; fitted: kappa, rho,
%   ell_F.  Fits to Fig. 2C with and without the tau operator, and the
%   unrefitted prediction of the APPwt/APPwt control (Fig. S4C).
%   MATLAB R2026a.

if nargin < 1, outfile = 'run_coupled_product_output.txt'; end
S0 = load('p211_spatial2.mat');  R2 = S0.R;
F = R2.F;  xg = R2.xg;  At = R2.At;  neff = R2.neff;
use = F.x >= 5;  xd = F.x(use);  [~, id] = ismember(xd, xg);
yd = [F.wt(use); F.swe(use)];  n = numel(yd);
L1 = load('p211_refit12_1.mat');  L2 = load('p211_refit12_2.mat');
lam = log(2)/(3.63*3600);
B = p211_spatial2('chains', {L1.R.p, L2.R.p}, lam, xg);
G = p211_figS4c();
A4 = interp1(G.x, G.tau, xg, 'pchip');  A4 = A4/A4(1);
k4 = G.x >= 5;  [~, i4] = ismember(G.x(k4), xg);  obs4 = [G.gfp(k4); G.rfp(k4)];
r = [1 1.45];  ell = [B.ellW B.ellS];

fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('P211  intact APP coupled to a diffusing tag fragment   %s   MATLAB %s\n', ...
   char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
pr('intact-APP lengths %.0f / %.0f um (population-consistent chains, t1/2 3.63 h); r = 1, 1.45\n', ell);
pr('measured Fig. 2C areas 56.0/43.0 54.0/34.0, distal 0.791, prox 0.964; S4C channel-to-channel rms %.3f\n', ...
   sqrt(mean((G.gfp(k4) - G.rfp(k4)).^2)));
opt = optimset('TolX',1e-9,'TolFun',1e-12,'MaxFunEvals',4e4,'MaxIter',4e4,'Display','off');
ops = {At, ones(size(xg))};  opn = {'tau operator', 'A = 1'};
R = struct([]);
for io = 1:2
    A = ops{io};
    pr('\n--- %s ---\n', opn{io});
    pr('   y     rms     areas wt      areas swe     distal  prox   dAIC(vs y=0) | rho  t_F/t_M  ell_F  D_F(t_M 3.63 h) | phi_w phi_s | S4C rms\n');
    a0 = NaN;
    for y = [0 0.2 0.5 1]
        best = inf;
        for st = {[log(0.8) log(3) log(500)], [log(0.8) log(1) log(200)], [log(0.8) log(10) log(1000)], [log(0.8) log(0.3) log(100)]}
            [u, v] = fminsearch(@(u) sum((sig2(u, y, A, B, r, ell, xg, id) - yd).^2), st{1}, opt);
            if v < best, best = v; ub = u; end
        end
        [~, fw, fs, sh] = sig2(ub, y, A, B, r, ell, xg, id);
        kap = exp(ub(1));  rho = exp(ub(2));  lF = exp(ub(3));
        iP = xg <= 100;  iD = xg >= 100;
        auc = [trapz(xg(iP),fw(iP)) trapz(xg(iD),fw(iD)) trapz(xg(iP),fs(iP)) trapz(xg(iD),fs(iD))];
        aic = neff*log(best/n) + 2*3;  if y == 0, a0 = aic; end
        % S4C: APPwt alone in both channels, its own tau profile, kappa fixed
        % the control is predicted with the same geometry rule as the fit:
        % its own tau profile when the fit used tau, none when it did not
        if all(A == 1), s4 = sh{1}(:); else, s4 = A4(:).*sh{1}(:); end
        p4 = kap*[s4(i4); s4(i4)];
        F0 = rho;  phw = F0/(1 + F0);  phs = 1.45*rho/(1 + 1.45*rho);
        tM = 3.63;  tF = rho*tM;  DF = lF^2*log(2)/(tF*3600);
        pr('  %4.1f  %.4f  %5.1f/%5.1f  %5.1f/%5.1f   %.3f  %.3f  %+6.1f      | %5.2f %5.2f  %5.0f  %6.1f | %.3f %.3f | %.4f\n', ...
           y, sqrt(best/n), auc, auc(4)/auc(2), auc(3)/auc(1), aic - a0, rho, rho, lF, DF, phw, phs, sqrt(mean((p4 - obs4).^2)));
        R(end+1).op = opn{io};  R(end).y = y;  R(end).rms = sqrt(best/n);  R(end).auc = auc; %#ok<AGROW>
        R(end).rho = rho;  R(end).lF = lF;  R(end).kap = kap;  R(end).s4c = sqrt(mean((p4 - obs4).^2));
    end
end
% APPwt/tau decline over the first 100 um implied by each y (tau operator removed)
pr('\nAPPwt intensity relative to geometry, fall over 0-100 um (best fit, tau operator): measured APPwt/tau 22%%\n');
for j = find(strcmp({R.op}, 'tau operator'))
    [~, ~, ~, sh] = sig2([log(R(j).kap) log(R(j).rho) log(R(j).lF)], R(j).y, At, B, r, ell, xg, id);
    pr('   y = %.1f: %.1f%%\n', R(j).y, 100*(1 - interp1(xg, sh{1}, 100)/sh{1}(1)));
end
fclose(fid);
save('p211_coupled_product.mat','R');
end

function [yv, fw, fs, sh] = sig2(u, y, A, B, r, ell, xg, id)
kap = exp(u(1));  rho = exp(u(2));  lF = exp(u(3));
c = {B.cW, B.cS};  sh = cell(1,2);
for g = 1:2
    C = c{g};                                  % C_g / C_g0
    den = 1 - lF^2/ell(g)^2;
    if abs(den) < 1e-6, den = 1e-6; end
    a = y*rho/den;
    Fg = a*C + (r(g)*rho - a)*exp(-xg/lF);
    sh{g} = (C + Fg)/(1 + r(g)*rho);
end
fw = kap*A.*sh{1};  fs = kap*A.*sh{2};
yv = [fw(id); fs(id)];
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
