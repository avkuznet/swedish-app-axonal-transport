function pi0 = p211_occ(p)
%P211_OCC  Closed-form stationary occupancies of the five-state chain.
%
%   The allowed transitions form a tree (A-P, R-P, P-D, D-S), so detailed
%   balance holds at steady state and the occupancies follow without
%   solving an eigenproblem:
%
%     pi_A/pi_P = kpa/kap,  pi_R/pi_P = kpr/krp,
%     pi_D/pi_P = kpd/kdp,  pi_S/pi_D = kds/ksd.
%
%   If the parameter struct carries kar or kra -- optional direct A <-> R
%   transitions, which the paper does not use -- the
%   graph contains a cycle, detailed balance no longer holds, and the
%   occupancies are taken from the null space of the generator instead.
%   The two routes agree to machine precision when kar = kra = 0.

if (isfield(p,'kar') && p.kar > 0) || (isfield(p,'kra') && p.kra > 0)
    pi0 = occ_null(p);
    return
end

wP = 1;
wA = p.kpa/p.kap;
wR = p.kpr/p.krp;
wD = p.kpd/p.kdp;
wS = wD*p.kds/p.ksd;
w  = [wA wR wP wD wS];
pi0 = w/sum(w);
end

% ------------------------------------------------------------------------
function pi0 = occ_null(p)
%OCC_NULL  Stationary distribution of the generator, cycle or no cycle.
Q = p211_gen(p);
A = [Q'; ones(1,5)];
b = [zeros(5,1); 1];
pi0 = (A\b)';
pi0 = max(pi0, 0);
pi0 = pi0/sum(pi0);
end
