function N = make_short_figures(outdir)
%MAKE_SHORT_FIGURES  Figures 2-5 and the numbers quoted with them in the short paper.
%
%   One parameter set throughout: the corrected fits of run_refit12
%   (pooled targets, exhaustive whole-track classes; APP movies sampled
%   from the open, turning-over population at t1/2 = 3.63 h).  For
%   lysosomes with APPswe the reference vector remained the best found
%   under the corrected rules and is used.  Writes figN.png (300 dpi) and
%   short_numbers.txt to outdir.  MATLAB R2026a, base MATLAB only.

if nargin < 1, outdir = fullfile('..','figures'); end
if ~exist(outdir,'dir'), mkdir(outdir); end
S0 = load('p211_fitted.mat');  P0 = S0.R.P;
Dp = p211_data_pooled();  D0 = p211_data();
lam = log(2)/(3.63*3600);
P = cell(1,6);
for k = 1:5, L = load(sprintf('p211_refit12_%d.mat',k)); P{k} = L.R.p; end
P{6} = P0{6};
names = {'APPwt','APPswe','Rab5','Rab5+APPswe','Lyso+APPwt','Lyso+APPswe'};
cg = [0.27 0.45 0.22];  cr = [0.85 0.22 0.15];  gry = [0.5 0.5 0.5];
mn = char(8722);   % ISO 80000-2 labels: italic quantities (\it), upright subscripts (\rm), true minus
um = [char(181) 'm'];   % upright micro sign (TeX \mu renders slanted)
fid = fopen(fullfile(outdir,'short_numbers.txt'),'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('Numbers for the short paper   %s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);

% ======================================================= fits and Table 1
M = cell(1,6);  T = cell(1,6);  rmsr = zeros(1,6);  occ = zeros(6,5);
for k = 1:6
    o = struct('N',20000,'vthr',Dp(k).vthr,'seed',99001,'exhaustive',true);
    if k <= 2, o.init = p211_spatcomp(P{k}, lam); end
    M{k} = p211_sim(P{k}, o);  T{k} = p211_compare(Dp(k), M{k});
    r = T{k}.rel(~isnan(T{k}.rel));  rmsr(k) = sqrt(mean(r.^2));
    occ(k,:) = p211_occ(P{k});
end
pr('\nTable 1 (constants; s^-1 unless stated)\n');
f = {'va','vr','kap','kpa','krp','kpr','kpd','kdp','kds','ksd'};
for k = 1:6
    p = P{k};  v = cellfun(@(q) p.(q), f);
    pr('%-12s %s | occ %s | veff %+.4f | rms %.1f%%\n', names{k}, sprintf('%.3g ', v), sprintf('%.3f ', occ(k,:)), ...
       p.va*occ(k,1) - p.vr*occ(k,2), 100*rmsr(k));
end
pr('\nTable 2 (APP statistics: measured / model)\n');
for j = 1:numel(T{1}.name)
    if ~isnan(T{1}.obs(j)), pr('%-22s %7.3f %7.3f | %7.3f %7.3f\n', T{1}.name{j}, T{1}.obs(j), T{1}.mod(j), T{2}.obs(j), T{2}.mod(j)); end
end
pr('\nSupplementary statistics, Rab5 and lysosomes (measured / model)\n');
for k = 3:6
    pr('%s\n', names{k});
    for j = 1:numel(T{k}.name)
        if ~isnan(T{k}.obs(j)), pr('   %-22s %7.3f %7.3f (%+.0f%%)\n', T{k}.name{j}, T{k}.obs(j), T{k}.mod(j), 100*T{k}.rel(j)); end
    end
end
nceil = 0;
for k = 1:6, nceil = nceil + (Dp(k).vsA > P{k}.va) + (Dp(k).vsR > P{k}.vr); end
pr('\nmeasured median segmental speed above the fitted run speed: %d of 12 directional comparisons\n', nceil);

% ---------------------------------------------------------- Fig. 2
fh = figure('Visible','off','Units','centimeters','Position',[0 0 17 7.5],'Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on
mk = {'o','o','s','s','^','^'};  cc = {cg,cr,cg,cr,cg,cr};
for k = 1:6
    t = T{k};  g = ~isnan(t.obs) & ~isnan(t.mod);
    plot(t.obs(g), t.mod(g), mk{k}, 'MarkerFaceColor', cc{k}, 'MarkerEdgeColor','k', 'MarkerSize',5);
end
z = [0.08 40];  plot(z, z, 'k--');  plot(z, 1.25*z, 'k:');  plot(z, z/1.25, 'k:');
set(gca,'XScale','log','YScale','log','FontSize',8,'Box','on','XTick',[0.1 1 10],'YTick',[0.1 1 10]);  xlim(z); ylim(z); axis square
xlabel('measured'); ylabel('model');  legend(names,'Location','southeast','FontSize',6.5,'Box','off');
text(0.11, 22, 'dotted: \pm25%', 'FontSize',7, 'Color', gry);
title('(a)','FontWeight','normal');
nexttile;
b = bar(categorical(names, names), 100*rmsr, 'FaceColor','flat');
for k = 1:6, b.CData(k,:) = cc{k}; end
ylabel('rms relative error (%)'); set(gca,'FontSize',8); title('(b)','FontWeight','normal'); ylim([0 15]);
isominus(fh); exportgraphics(fh, fullfile(outdir,'fig2.png'), 'Resolution', 300); close(fh);

% ======================================================= Fig. 3: motion without progress
pw = P{1};  ps = P{2};  ow = occ(1,:);  os = occ(2,:);
vw = pw.va*ow(1) - pw.vr*ow(2);  vs = ps.va*os(1) - ps.vr*os(2);
dsum = @(d) (1 - d.fracS)*(d.thetaA*d.vsA - d.thetaR*d.vsR);
U = p211_suppdata();  V = p211_figS1c();
s1 = @(g) [U.S1.(g).frac(1)*V.ant.(g)(2) - U.S1.(g).frac(3)*V.ret.(g)(2), ...
           U.S1.(g).frac(1)*mean(V.ant.(g)) - U.S1.(g).frac(3)*mean(V.ret.(g))];
sw = s1('wt');  ss = s1('swe');
nrm = @(f) f/sum(f);  ex = @(f) [1 0 -1]*nrm(f(:));
exW = [ex(U.main.wt) ex(U.S2.gfp) ex(U.S2.trfp) ex(U.S1.wt.frac) ex(U.S3.wt)];
exS = [ex(U.main.swe) ex(U.S1.swe.frac) ex(U.S3.swe)];
mot = [ow(1)+ow(2), os(1)+os(2)];  runA = [pw.va/pw.kap, ps.va/ps.kap];  runR = [pw.vr/pw.krp, ps.vr/ps.krp];
pr('\nFig. 3\n   drift (analytic) %+.4f / %+.4f um/s, ratio %.0f\n', vw, vs, vw/vs);
pr('   summary-statistic drift %+.4f / %+.4f\n', dsum(D0(1)), dsum(D0(2)));
pr('   Fig. S1 velocity-weighted drift: APPwt %+.3f to %+.3f, APPswe %+.3f to %+.3f\n', min(sw), max(sw), min(ss), max(ss));
pr('   track excess APPwt %s; APPswe %s\n', sprintf('%+.2f ', exW), sprintf('%+.2f ', exS));
pr('   motile occupancy %.3f -> %.3f (%.1f-fold); va %.3f -> %.3f (%.0f%%); vr %.3f -> %.3f\n', mot, mot(1)/mot(2), pw.va, ps.va, 100*ps.va/pw.va, pw.vr, ps.vr);
pr('   anterograde run scale %.2f -> %.2f um (%.1f-fold); retrograde %.2f -> %.2f um\n', runA, runA(1)/runA(2), runR);
pr('   piA/piR %.2f -> %.2f; stationary occupancy %.3f -> %.3f\n', ow(1)/ow(2), os(1)/os(2), ow(5), os(5));
fh = figure('Visible','off','Units','centimeters','Position',[0 0 17 6.2],'Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on
est = [vw vs; dsum(D0(1)) dsum(D0(2)); mean(sw) mean(ss)];
bb = bar(est, 'grouped');  bb(1).FaceColor = cg;  bb(2).FaceColor = cr;
xe = [bb(1).XEndPoints; bb(2).XEndPoints];
errorbar(xe(:,3)', est(3,:), est(3,:) - [min(sw) min(ss)], [max(sw) max(ss)] - est(3,:), 'k.', 'LineStyle','none');
yline(0,'k-'); set(gca,'XTick',1:3,'XTickLabel',{'model','summary','Fig. S1'},'FontSize',8,'Box','on');
ylabel(['net drift (' um ' s^{' mn '1})']); legend({'APPwt','APPswe'},'Location','southwest','FontSize',7,'Box','off');
title('(a)','FontWeight','normal');
nexttile; hold on
plot(1, exW(1), 'o', 'MarkerFaceColor', cg, 'MarkerEdgeColor','k');  plot(1+0.08*(1:4), exW(2:5), 'o', 'MarkerEdgeColor', cg, 'LineWidth',1);
plot(2, exS(1), 'o', 'MarkerFaceColor', cr, 'MarkerEdgeColor','k');  plot(2+0.08*(1:2), exS(2:3), 'o', 'MarkerEdgeColor', cr, 'LineWidth',1);
yline(0,'k-'); xlim([0.6 2.6]); set(gca,'XTick',[1 2],'XTickLabel',{'APPwt','APPswe'},'FontSize',8,'Box','on');
ylabel(['\itf\rm_A ' mn ' \itf\rm_R (track classes)']); text(0.7, -0.12, 'filled: fitted; open: unfitted', 'FontSize',6.5, 'Color', gry);
title('(b)','FontWeight','normal');
nexttile; hold on
fc = [mot(2)/mot(1), ps.va/pw.va, ps.vr/pw.vr, runA(2)/runA(1), runR(2)/runR(1), (os(1)/os(2))/(ow(1)/ow(2)), vs/vw];
lab = {'\it\pi\rm_A + \it\pi\rm_R','\itv\rm_a','\itv\rm_r','\itv\rm_a/\itk\rm_{ap}','\itv\rm_r/\itk\rm_{rp}','\it\pi\rm_A/\it\pi\rm_R','\itv\rm_{eff}'};
barh(1:7, log10(fc), 'FaceColor', [0.6 0.6 0.6]);  xline(0,'k-');
set(gca,'YTick',1:7,'YTickLabel',lab,'YDir','reverse','FontSize',8,'Box','on');
xlabel('log_{10}(APPswe / APPwt)'); title('(c)','FontWeight','normal');
isominus(fh); exportgraphics(fh, fullfile(outdir,'fig3.png'), 'Resolution', 300); close(fh);

% ======================================================= Fig. 4: what the images see
S2 = load('p211_spatial2.mat');  R2 = S2.R;  F = R2.F;  xg = R2.xg;  At = R2.At;
use = F.x >= 5;  [~, id] = ismember(F.x(use), xg);  yd = [F.wt(use); F.swe(use)];
B = p211_spatial2('chains', {P{1}, P{2}}, lam, xg);  cW = B.cW;  cS = B.cS;
G = p211_figS4c();  A4 = interp1(G.x, G.tau, xg, 'pchip');  A4 = A4/A4(1);
k4 = G.x >= 5;  [~, i4] = ismember(G.x(k4), xg);  obs4 = [G.gfp(k4); G.rfp(k4)];
r = 1.45;  sg = @(z) 1./(1 + exp(-z));
opt = optimset('TolX',1e-9,'TolFun',1e-12,'MaxFunEvals',4e4,'MaxIter',4e4,'Display','off');
% transport only (kappa fitted)
y0 = [At(id).*cW(id); At(id).*cS(id)];  k0 = (y0'*yd)/(y0'*y0);
% soma-derived fragment, odds_swe = r*odds_wt
mdl = @(u) deal((1-sg(u(2)))*cW + sg(u(2))*exp(-xg/exp(u(3))), ...
                (1 - r*sg(u(2))/(1-sg(u(2)))/(1 + r*sg(u(2))/(1-sg(u(2)))))*cS + ...
                (r*sg(u(2))/(1-sg(u(2)))/(1 + r*sg(u(2))/(1-sg(u(2)))))*exp(-xg/exp(u(3))));
best = inf;
for st = {[log(0.8) 0 log(300)], [log(0.8) 1 log(600)], [log(0.8) -1 log(100)], [log(0.8) 2 log(1000)]}
    [u, v] = fminsearch(@(u) obj4(u, mdl, At, id, yd), st{1}, opt);  if v < best, best = v; ub = u; end
end
[w4, s4] = mdl(ub);  kap = exp(ub(1));  phw = sg(ub(2));  lF = exp(ub(3));  ow4 = phw/(1-phw);  phs = r*ow4/(1+r*ow4);
iP = xg <= 100;  iD = xg >= 100;  ar = @(v) [trapz(xg(iP),v(iP)) trapz(xg(iD),v(iD))];
a0w = ar(k0*At.*cW);  a0s = ar(k0*At.*cS);  a4w = ar(kap*At.*w4);  a4s = ar(kap*At.*s4);
pc = kap*A4.*w4;  p4 = [pc(i4); pc(i4)];
lw = p211_profile(P{1});  ls = p211_profile(P{2});
pr('\nFig. 4\n   transport only: ell %.0f / %.0f um; areas %.1f/%.1f %.1f/%.1f; ratios prox %.3f dist %.3f; rms %.3f\n', ...
   lw.ell, ls.ell, a0w, a0s, a0s(1)/a0w(1), a0s(2)/a0w(2), sqrt(mean((k0*y0 - yd).^2)));
pr('   fragment model: kappa %.3f phi_w %.3f phi_s %.3f ell_F %.0f; areas %.1f/%.1f %.1f/%.1f; ratios prox %.3f dist %.3f; rms %.4f\n', ...
   kap, phw, phs, lF, a4w, a4s, a4s(1)/a4w(1), a4s(2)/a4w(2), sqrt(best/numel(yd)));
pr('   WT/WT control, no refit: rms %.4f (channel-to-channel %.3f)\n', sqrt(mean((p4 - obs4).^2)), sqrt(mean((G.gfp(k4)-G.rfp(k4)).^2)));
pr('   J/C = lam*ell: %.4f / %.4f um/s, ratio %.3f; images as transport profiles 216/379 = %.2f\n', lam*lw.ell, lam*ls.ell, ls.ell/lw.ell, 216/379);
fh = figure('Visible','off','Units','centimeters','Position',[0 0 17 13],'Color','w');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on
plot(F.x, F.wt, 'o', 'Color', cg, 'MarkerSize',3);  plot(F.x, F.swe, 's', 'Color', cr, 'MarkerSize',3);  plot(F.x, F.tau, '.', 'Color', gry);
plot(xg, k0*At.*cW, '--', 'Color', cg, 'LineWidth',1.3);  plot(xg, k0*At.*cS, '--', 'Color', cr, 'LineWidth',1.3);
xlabel(['\itx\rm (' um ')']); ylabel('\itf\rm/\itf\rm_0'); ylim([0 1.5]); set(gca,'FontSize',8,'Box','on');
legend({'APPwt','APPswe','tau','transport only'},'FontSize',6.5,'Box','off','Location','northeast'); title('(a)','FontWeight','normal');
nexttile; hold on
plot(F.x, F.wt, 'o', 'Color', cg, 'MarkerSize',3);  plot(F.x, F.swe, 's', 'Color', cr, 'MarkerSize',3);
plot(xg, kap*At.*w4, '-', 'Color', cg, 'LineWidth',1.3);  plot(xg, kap*At.*s4, '-', 'Color', cr, 'LineWidth',1.3);
plot(xg, kap*At.*phw.*exp(-xg/lF), ':', 'Color', cg, 'LineWidth',1.2);  plot(xg, kap*At.*phs.*exp(-xg/lF), ':', 'Color', cr, 'LineWidth',1.2);
xlabel(['\itx\rm (' um ')']); ylabel('\itf\rm/\itf\rm_0'); ylim([0 1.1]); set(gca,'FontSize',8,'Box','on');
legend({'APPwt','APPswe','model','','fragment'},'FontSize',6.5,'Box','off','Location','northeast'); title('(b)','FontWeight','normal');
nexttile; hold on
plot(G.x, G.gfp, 'o', 'Color', cg, 'MarkerSize',3);  plot(G.x, G.rfp, 'd', 'Color', [0.7 0.2 0.7], 'MarkerSize',3);  plot(G.x, G.tau, '.', 'Color', gry);
plot(xg, pc, '-', 'Color', 'k', 'LineWidth',1.3);
xlabel(['\itx\rm (' um ')']); ylabel('\itf\rm/\itf\rm_0'); ylim([0 1.5]); set(gca,'FontSize',8,'Box','on');
legend({'APPwt-GFP','APPwt-tRFP','tau','prediction'},'FontSize',6.5,'Box','off','Location','northeast'); title('(c)','FontWeight','normal');
nexttile;
bv = [ls.ell/lw.ell, 216/379];
bc = bar(1:2, bv, 0.55, 'FaceColor','flat');  bc.CData = [0.4 0.4 0.4; 0.8 0.8 0.8];
set(gca,'XTick',1:2,'XTickLabel',{'movies','images as transport'},'FontSize',8,'Box','on','YScale','log'); ylim([0.01 1]);
ylabel('(\itJ\rm/\itC\rm)_{swe}/(\itJ\rm/\itC\rm)_{wt}');
text(1:2, bv*1.35, compose('%.2f', bv), 'HorizontalAlignment','center','FontSize',8); title('(d)','FontWeight','normal');
isominus(fh); exportgraphics(fh, fullfile(outdir,'fig4.png'), 'Resolution', 300); close(fh);

% ======================================================= Fig. 5: motor-pool test
Rr = cell(1,6);
for k = 3:6
    L = load(sprintf('p211_refit12_%d.mat',k)); Rr{k} = L.R;
end
o6 = struct('N',20000,'vthr',Dp(6).vthr,'seed',99001,'exhaustive',true);
pairs = {[3 4],[5 6]};  rng_k = zeros(2,2);  rng_f = zeros(2,2);  pt_k = zeros(1,2);  pt_f = zeros(1,2);
for i = 1:2
    a = pairs{i}(1);  b = pairs{i}(2);
    ca = cands(Rr{a}, P0{a}, P{a});  cb = cands(Rr{b}, P0{b}, P{b});
    rr = [];  ff = [];
    for u = 1:numel(ca), for v = 1:numel(cb)
        oa = p211_occ(ca{u});  ob = p211_occ(cb{v});
        rr(end+1) = cb{v}.kpr/ca{u}.kpr;  ff(end+1) = (cb{v}.kpr*ob(3))/(ca{u}.kpr*oa(3)); %#ok<AGROW>
    end, end
    pt_k(i) = P{b}.kpr/P{a}.kpr;  pt_f(i) = (P{b}.kpr*occ(b,3))/(P{a}.kpr*occ(a,3));
    rng_k(i,:) = [min(rr) max(rr)];  rng_f(i,:) = [min(ff) max(ff)];
end
red = arrayfun(@(mi) load(sprintf('p211_rab5red_%d_%d.mat',3,mi)), 1:2);  red2 = arrayfun(@(mi) load(sprintf('p211_rab5red_%d_%d.mat',4,mi)), 1:2);
rk3 = [red2(1).R.p.kpr/red(1).R.p.kpr, red2(2).R.p.kpr/red(2).R.p.kpr, pt_k(1)];
pr('\nFig. 5\n   kpr ratio Rab5 %.2f (%.2f-%.2f), lysosomes %.2f (%.2f-%.2f); APP %.2f\n', pt_k(1), rng_k(1,:), pt_k(2), rng_k(2,:), P{2}.kpr/P{1}.kpr);
pr('   event-flux ratio Rab5 %.2f (%.2f-%.2f), lysosomes %.2f (%.2f-%.2f)\n', pt_f(1), rng_f(1,:), pt_f(2), rng_f(2,:));
pr('   Rab5 kpr ratio with 3-, 4-, 5-state chains: %s\n', sprintf('%.2f ', rk3));
fh = figure('Visible','off','Units','centimeters','Position',[0 0 12 7],'Color','w'); hold on
patch([0.4 2.6 2.6 0.4], [0.1 0.1 1 1], [0.93 0.93 0.93], 'EdgeColor','none');
text(0.5, 0.15, 'fixed shared pool predicts < 1', 'FontSize',7, 'Color', gry);
for i = 1:2
    h1 = errorbar(i-0.12, pt_k(i), pt_k(i)-rng_k(i,1), rng_k(i,2)-pt_k(i), 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
    h2 = errorbar(i+0.12, pt_f(i), pt_f(i)-rng_f(i,1), rng_f(i,2)-pt_f(i), 'ks', 'MarkerFaceColor','w', 'MarkerSize',5);
end
h3 = plot(1-0.12+[-0.1 -0.05], rk3(1:2), 'k^', 'MarkerSize',4);
yline(1,'k-'); set(gca,'YScale','log','XTick',[1 2],'XTickLabel',{'Rab5','lysosomes'},'FontSize',8,'Box','on'); xlim([0.4 2.6]); ylim([0.1 20]);   % headroom keeps the legend clear of the Rab5 bar (top 4.97)
ylabel('APPswe / control');
legend([h1 h2 h3], {'\itk\rm_{pr} (entry rate)','\itk\rm_{pr}\it\pi\rm_P (event flux)','\itk\rm_{pr}, 3- and 4-state chains'},'FontSize',6.5,'Box','off','Location','northwest');
isominus(fh); exportgraphics(fh, fullfile(outdir,'fig5.png'), 'Resolution', 300); close(fh);
fclose(fid);
N.rms = rmsr;  N.P = P;
end

% ------------------------------------------------------------------------
function v = obj4(u, mdl, At, id, yd)
[w, s] = mdl(u);  kap = exp(u(1));
y = [kap*At(id).*w(id); kap*At(id).*s(id)];  v = sum((y - yd).^2);
end

function C = cands(R, pref, pbest)
% candidate vectors within 1.5x of the best objective, including the reference
f = [R.cand.fchk];  fb = min([f R.f0]);
C = {R.cand(f <= 1.5*fb).p};
if R.f0 <= 1.5*fb, C{end+1} = pref; end
if isempty(C), C = {pbest}; end
end

function isominus(fh)
% ISO 80000-2: freeze tick labels and replace hyphen-minus before a digit
% (including TeX exponents such as 10^{-1}) by the minus sign U+2212.
drawnow;
for ax = findall(fh, 'Type', 'axes')'
    for r = {'XTickLabel', 'YTickLabel'}
        lbl = get(ax, r{1});
        if isempty(lbl), continue; end
        set(ax, r{1}, regexprep(cellstr(lbl), '-(?=\d)', char(8722)));
    end
end
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
