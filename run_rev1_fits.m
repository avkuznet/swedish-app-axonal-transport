function run_rev1_fits(job, k, arg)
%RUN_REV1_FITS  Reduced chains, run-speed variability and perturbation refits (Sec. 2.3, Appendix).
%
%   run_rev1_fits('reduced', k, mi)   k = 1,2,5,6; mi = 1 (A,R,P) or 2 (A,R,P,S)
%                                     -> p211_red_<k>_<mi>.mat
%   run_rev1_fits('cv', k)            five states + fitted lognormal run-speed CV
%                                     -> p211_cv_<k>.mat
%   run_rev1_fits('boot', k, b)       perturbation refits b (vector): targets
%                                     multiplied by exp(0.1 z), compositions
%                                     renormalized, five-state refit from the
%                                     corrected fit -> p211_boot_<k>_<b>.mat
%   Rab5 reduced chains are in p211_rab5red_* (run_rab5_reduced).
%   All use the corrected rules; APP (k = 1, 2) is sampled from the open
%   population at t1/2 = 3.63 h.  MATLAB R2026a.

lam = log(2)/(3.63*3600);
Dp = p211_data_pooled();  d = Dp(k);
L = load(sprintf('p211_refit12_%d.mat', k));  p0 = L.R.p;
if k == 6, S0 = load('p211_fitted.mat'); p0 = S0.R.P{6}; end
opt = struct();  if k <= 2, opt.lam = lam; end
switch job
    case 'reduced'
        opt.mi = arg;  opt.nstart = 5;
        R = p211_fitx(d, p0, opt);
        save(sprintf('p211_red_%d_%d.mat', k, arg), 'R');
        fprintf('%s mi=%d: objective %.4f rms %.1f%%\n', d.name, arg, R.f, 100*R.rms);
    case 'cv'
        opt.fitcv = true;  opt.nstart = 4;
        R = p211_fitx(d, p0, opt);
        save(sprintf('p211_cv_%d.mat', k), 'R');
        fprintf('%s cv: objective %.4f rms %.1f%% cv %.3f\n', d.name, R.f, 100*R.rms, R.p.speedcv);
    case 'boot'
        opt.nstart = 0;  opt.maxfun = 800;
        f = {'thetaA','thetaR','thetaP','fracA','fracS','fracR','LA','LR','fpA','fpR','frev','vsA','vsR','Lpool','fpPool'};
        for b = arg
            rs = RandStream('mt19937ar','Seed',7000 + 100*k + b);
            db = d;
            for i = 1:numel(f)
                if isfield(db, f{i}) && ~isnan(db.(f{i})), db.(f{i}) = db.(f{i})*exp(0.1*randn(rs)); end
            end
            s = db.thetaA + db.thetaR + db.thetaP;
            db.thetaA = db.thetaA/s;  db.thetaR = db.thetaR/s;  db.thetaP = db.thetaP/s;
            if ~isnan(db.fracA)
                s = db.fracA + db.fracS + db.fracR;
                db.fracA = db.fracA/s;  db.fracS = db.fracS/s;  db.fracR = db.fracR/s;
            end
            R = p211_fitx(db, p0, opt);  R.target = db;
            save(sprintf('p211_boot_%d_%02d.mat', k, b), 'R');
            fprintf('%s boot %d: objective %.4f\n', d.name, b, R.f);
        end
end
end
