function report_rab5_reduced(outfile)
%REPORT_RAB5_REDUCED  3-, 4- and 5-state chains for Rab5 under the corrected targets.
%   Reads p211_rab5red_<k>_<mi>.mat (run_rab5_reduced) and the 5-state
%   refits p211_refit12_<k>.mat; reports objective, rms, retrograde-entry
%   ratio kpr and event-flux ratio kpr*piP (APPswe condition / control),
%   and drift.  MATLAB R2026a.
if nargin < 1, outfile = 'report_rab5_reduced_output.txt'; end
fid = fopen(outfile,'w');
pr = @(varargin) fprintf2(fid, varargin{:});
pr('P211  reduced chains for Rab5 (corrected targets: 8 distinct, 7 independent)   %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')));
lab = {'3-state A,R,P (6 constants)','4-state A,R,P,S (8)','5-state (10)'};
for mi = 1:3
    for k = [3 4]
        if mi < 3, L = load(sprintf('p211_rab5red_%d_%d.mat', k, mi)); P{k} = L.R.p; f(k) = L.R.f; rr(k) = L.R.rms; %#ok<AGROW>
        else, L = load(sprintf('p211_refit12_%d.mat', k)); P{k} = L.R.p; f(k) = L.R.fchk; rr(k) = L.R.rms; end
        o{k} = p211_occ(P{k}); %#ok<AGROW>
    end
    pr('%-28s objective %.4f / %.4f  rms %.1f%% / %.1f%% | kpr ratio %.2f | event-flux ratio %.2f | v_eff %+.4f / %+.4f\n', ...
       lab{mi}, f(3), f(4), 100*rr(3), 100*rr(4), P{4}.kpr/P{3}.kpr, (P{4}.kpr*o{4}(3))/(P{3}.kpr*o{3}(3)), ...
       P{3}.va*o{3}(1)-P{3}.vr*o{3}(2), P{4}.va*o{4}(1)-P{4}.vr*o{4}(2));
end
fclose(fid);
end
function fprintf2(fid, varargin)
fprintf(varargin{:});  fprintf(fid, varargin{:});
end
