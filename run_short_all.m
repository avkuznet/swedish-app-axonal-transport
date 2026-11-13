function run_short_all(stage)
%RUN_SHORT_ALL  Reproduce the figures and numbers of the short paper from this folder.
%
%   run_short_all              figures and analyses from the saved fits (minutes)
%   run_short_all('figures')   Figs. 1-5 and short_numbers.txt only
%   run_short_all('analyses')  the analysis scripts whose numbers the text quotes
%   run_short_all('refit')     ALSO refits all six conditions, the reduced chains,
%                              the run-speed CV fits and the perturbation refits (many hours;
%                              overwrites the saved .mat fits).  Each MATLAB session
%                              needs about 1-1.5 GB: run at most 3 sessions at once.
%
%   Run from this folder.  Figures go to ../figures.  MATLAB R2026a, base MATLAB.

if nargin < 1, stage = 'all'; end
t0 = tic;
if strcmp(stage, 'refit')
    fprintf('refitting the six conditions (run_refit12) ...\n');
    for k = 1:6, run_refit12(k); end
    fprintf('reduced Rab5 chains (run_rab5_reduced) ...\n');
    for k = [3 4], for mi = [1 2], run_rab5_reduced(k, mi); end, end
    fprintf('reduced chains, run-speed CV and perturbation refits (run_rev1_fits) ...\n');
    for k = [1 2 5 6], for mi = [1 2], run_rev1_fits('reduced', k, mi); end, end
    for k = 1:6, run_rev1_fits('cv', k); end
    % Ten perturbation refits per condition.  The index sets the random
    % seed of the perturbed targets.  The saved set for APPwt (k = 1) has
    % indices 1-4 and 7-12, not 1-10, so the same indices are used here.
    bset = {[1:4 7:12], 1:10, [], [], 1:10, 1:10};
    for k = [1 2 5 6], run_rev1_fits('boot', k, bset{k}); end
    stage = 'all';
end
if any(strcmp(stage, {'all','analyses'}))
    report_refit12();               % corrected fits, motor-entry ratios, consistent spatial lengths
    report_rab5_reduced();          % 3-, 4-, 5-state Rab5 chains
    run_rev12_diag();               % chain vs open-population composition
    run_drift_replicates();         % drift in recordings no fit used
    run_finite_axon_consistent();   % terminal boundary conditions
    run_product_model();            % soma-derived product, held-out control
    run_coupled_product();          % product coupled to transported APP; no-tau check
    report_rev1();                  % model selection, run-speed CV, perturbation-refit ranges
end
if any(strcmp(stage, {'all','figures'}))
    make_fig1();
    make_short_figures();
end
fprintf('run_short_all(''%s'') finished in %.0f s\n', stage, toc(t0));
end
