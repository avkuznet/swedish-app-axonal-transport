function R = run_drift_replicates(outfile)
%RUN_DRIFT_REPLICATES  Is the collapse of net APPswe progress replicated?
%
%   The drift contrast of Sec. 3.2 rests on one set of recordings, the
%   lentiviral 2 frames/s movies of Fig. 1D and Fig. 3 of Feole et al.
%   (2024); the model fit and the summary-statistic estimate both use
%   them.  That study contains three further APP recordings that no fit used:
%     Fig. S1B-D  plasmid expression, 4 frames/s, 60 s (APPwt, APPswe)
%     Fig. S2C    APPwt under GFP and under tRFP, 30 s
%     Fig. S3B    axons in the axonal compartment of microfluidic chambers,
%                 far from the cell bodies, 30 s (APPwt, APPswe)
%   This script asks what they say about net progress.
%
%   (1) Directional track excess fA - fR, the fraction of tracks called
%       anterograde minus the fraction called retrograde.  Model-free and
%       available for every recording; its sign is the sign of the net
%       direction of the moving tracks, though it is not a velocity.
%   (2) The summary-statistic drift of Sec. 3.2 for Fig. 1D,
%           v = (1 - fS) (thetaA vsA - thetaR vsR),
%       reproduced as a check.
%   (3) A drift for Fig. S1, which reports per-particle average
%       velocities (Fig. S1C) instead of segmental speeds:
%           v = fA vbarA - fR vbarR,
%       stationary tracks contributing nothing.  vbar is estimated from
%       the box plots (median, and (Q1 + median + Q3)/3, the quartile
%       estimator of the mean); the range uses the quartiles.
%   (4) Calibration of estimator (3) on simulated tracks of the reference
%       chains (p211_fitted.mat, the starting points of run_refit12, not
%       the adopted fits) passed through the operator at the Fig. S1 settings, taking
%       a particle's average velocity to be its net displacement over the
%       movie divided by the duration.  With all 20,000 tracks it is
%       compared with v_eff; drawing 10 particles per class, as in Fig.
%       S1C, gives the sampling spread of the estimate.
%
%   The paper quotes parts (1)-(3), which use the published data only.
%   The model lines of parts (1) and (4) use the reference chains.
%
%   MATLAB R2026a, base MATLAB only.

if nargin < 1, outfile = 'run_drift_replicates_output.txt'; end
S = load('p211_fitted.mat');  D = S.R.D;  P = S.R.P;
U = p211_suppdata();
V = p211_figS1c();          % digitised Fig. S1C quartiles (from digitize_figS1c)
fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('===========================================================\n');
pr(' P211  replication of the APP drift contrast in unfitted recordings of Feole et al. (2024)\n');
pr(' %s   MATLAB %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
pr('===========================================================\n');

% ------------------------------------------------------------ (1)
nrm = @(f) f/sum(f);
rec = { 'APPwt ', 'Fig. 1D, lentiviral, 2 fps, 30 s (fitted)', U.main.wt
        'APPwt ', 'Fig. S2C GFP, 2 fps, 30 s',                  U.S2.gfp
        'APPwt ', 'Fig. S2C tRFP, 2 fps, 30 s',                 U.S2.trfp
        'APPwt ', 'Fig. S1B, plasmid, 4 fps, 60 s',             U.S1.wt.frac
        'APPwt ', 'Fig. S3B, microfluidic, distal axons',       U.S3.wt
        'APPswe', 'Fig. 1D, lentiviral, 2 fps, 30 s (fitted)', U.main.swe
        'APPswe', 'Fig. S1B, plasmid, 4 fps, 60 s',             U.S1.swe.frac
        'APPswe', 'Fig. S3B, microfluidic, distal axons',       U.S3.swe };
pr('\n(1) directional track excess fA - fR  (classes renormalised to sum to one)\n');
R.excess = zeros(size(rec,1),1);
for i = 1:size(rec,1)
    f = nrm(rec{i,3});  R.excess(i) = f(1) - f(3);
    pr('   %s  %-44s  fA %.3f fS %.3f fR %.3f   excess %+.3f\n', rec{i,1}, rec{i,2}, f, R.excess(i));
end
pr('   range, APPwt  %+.3f to %+.3f;  APPswe %+.3f to %+.3f\n', ...
   min(R.excess(1:5)), max(R.excess(1:5)), min(R.excess(6:8)), max(R.excess(6:8)));
o2 = struct('N',20000,'seed',99001);
o4 = struct('N',20000,'fps',U.S1.fps,'T',U.S1.T,'seed',99001);
for k = 1:2
    o2.vthr = D(k).vthr;  o4.vthr = D(k).vthr;
    m2 = p211_sim(P{k}, o2);  m4 = p211_sim(P{k}, o4);
    R.model(k).m2 = m2;  R.model(k).m4 = m4;
    pr('   reference chain %-6s excess %+.3f at 2 fps, 30 s; %+.3f at 4 fps, 60 s   (v_eff %+.4f um/s)\n', ...
       D(k).name, m2.fracA - m2.fracR, m4.fracA - m4.fracR, m2.veff);
end

% ------------------------------------------------------------ (2)
pr('\n(2) summary-statistic drift for the fitted recordings (Sec. 3.2)\n');
for k = 1:2
    t = D(k);
    v = (1 - t.fracS)*(t.thetaA*t.vsA - t.thetaR*t.vsR);
    pr('   %-6s  (1 - %.2f)(%.3f x %.2f - %.3f x %.2f) = %+.4f um/s\n', D(k).name, ...
       t.fracS, t.thetaA, t.vsA, t.thetaR, t.vsR, v);
    R.vsummary(k) = v;
end

% ------------------------------------------------------------ (3)
pr('\n(3) Fig. S1: track classes (S1B) weighted by per-particle average velocities (S1C, n = 10 per class)\n');
nm = {'wt','swe'};
for k = 1:2
    f = U.S1.(nm{k}).frac;  a = V.ant.(nm{k});  r = V.ret.(nm{k});
    pr('   APP%-3s  anterograde [Q1 median Q3] = [%.2f %.2f %.2f], retrograde [%.2f %.2f %.2f] um/s\n', nm{k}, a, r);
    vMed = f(1)*a(2) - f(3)*r(2);
    vQm  = f(1)*mean(a) - f(3)*mean(r);
    vLo  = f(1)*a(1) - f(3)*r(3);
    vHi  = f(1)*a(3) - f(3)*r(1);
    pr('           drift: median-based %+.3f, quartile-mean %+.3f, quartile range [%+.3f, %+.3f] um/s\n', vMed, vQm, vLo, vHi);
    R.S1(k) = struct('vMed',vMed,'vQm',vQm,'vLo',vLo,'vHi',vHi);
end

% ------------------------------------------------------------ (4)
pr('\n(4) estimator (3) on simulated tracks of the reference chains (p211_fitted.mat) at 4 fps, 60 s\n');
nb = 4000;  rng(20261106);
for k = 1:2
    m = R.model(k).m4;  X = m.X;  T = U.S1.T;  dtf = 1/U.S1.fps;
    spd = mean(abs(diff(X,1,1)/dtf),1);  vnet = (X(end,:) - X(1,:))/T;
    isA = spd >= 0.1 & vnet > 0.01;  isR = spd >= 0.1 & vnet < -0.01;
    fA = mean(isA);  fR = mean(isR);  va = vnet(isA);  vr = -vnet(isR);
    qa = quantile(va,[.25 .5 .75]);  qr = quantile(vr,[.25 .5 .75]);
    vAll = mean(vnet);
    pr('   %-6s  v_eff %+.4f | mean net velocity of all tracks %+.4f | fA mean - fR mean %+.4f | median-based %+.4f | quartile-mean %+.4f um/s\n', ...
       D(k).name, m.veff, vAll, fA*mean(va) - fR*mean(vr), fA*qa(2) - fR*qr(2), fA*mean(qa) - fR*mean(qr));
    % ten particles per class, as in Fig. S1C; class fractions as reported
    est = zeros(nb,2);
    for b = 1:nb
        sa = va(randi(numel(va),10,1));  sr = vr(randi(numel(vr),10,1));
        est(b,1) = fA*median(sa) - fR*median(sr);
        est(b,2) = fA*mean(quantile(sa,[.25 .5 .75])) - fR*mean(quantile(sr,[.25 .5 .75]));
    end
    ci = quantile(est,[.025 .975]);
    pr('           10 particles per class: median-based 95%% range [%+.3f, %+.3f], quartile-mean [%+.3f, %+.3f]; P(estimate < 0) = %.2f / %.2f\n', ...
       ci(1,1), ci(2,1), ci(1,2), ci(2,2), mean(est(:,1) < 0), mean(est(:,2) < 0));
    R.boot(k).ci = ci;  R.boot(k).pneg = mean(est < 0);
    R.boot(k).pObs = [mean(est(:,1) <= R.S1(k).vMed) mean(est(:,2) <= R.S1(k).vQm)];
    pr('           fraction of draws at or below the Fig. S1 value: %.3f (median-based), %.3f (quartile-mean)\n', R.boot(k).pObs);
end
fclose(fid);
end

function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
