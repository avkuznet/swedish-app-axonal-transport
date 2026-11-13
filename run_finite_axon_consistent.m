function run_finite_axon_consistent(outfile)
%RUN_FINITE_AXON_CONSISTENT  Effect of the terminal boundary conditions on the imaged segment (Sec. 2.4).
%
%   Uses the adopted APP chains (p211_refit12_1, _2).  Reports the decay
%   lengths of the terminal and soma modes, then, for axon lengths of
%   250-1000 um, the largest relative change of the 0-100 and 100-200 um
%   window areas over terminal return fractions rho = 0, 0.5, 1 and distal
%   retrograde supplies jin = 0-2 (p211_profile_finite), relative to the
%   semi-infinite axon.  MATLAB R2026a.
if nargin < 1, outfile = 'run_finite_axon_consistent_output.txt'; end
L1 = load('p211_refit12_1.mat');  L2 = load('p211_refit12_2.mat');  P = {L1.R.p, L2.R.p};
fid = fopen(outfile,'w');
q1 = p211_profile_finite(P{1}, struct('Lax',Inf));  q2 = p211_profile_finite(P{2}, struct('Lax',Inf));
fprintf(fid,'terminal mode lengths: APPwt %.1f um, APPswe %.1f um; soma mode %.0f / %.0f um; balance %.1e %.1e\n', q1.ell2, q2.ell2, q1.ell1, q2.ell1, q1.balance, q2.balance);
base = [q1.aucProx q1.aucDist q2.aucProx q2.aucDist];
for La = [250 300 400 500 1000]
    mx = 0;
    for rho = [0 0.5 1], for ji = [0 0.1 0.25 0.5 1 2]
        w = p211_profile_finite(P{1}, struct('Lax',La,'rho',rho,'jin',ji));
        s = p211_profile_finite(P{2}, struct('Lax',La,'rho',rho,'jin',ji));
        mx = max(mx, max(abs([w.aucProx w.aucDist s.aucProx s.aucDist]./base - 1)));
    end, end
    fprintf(fid,'Lax %4d um: largest relative change of a window area %.4f\n', La, mx);
end
fclose(fid);  type(outfile);
end
