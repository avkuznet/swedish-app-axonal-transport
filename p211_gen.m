function Q = p211_gen(p)
%P211_GEN  Generator matrix of the five-state cargo chain.
%
%   States  1 A, 2 R, 3 P, 4 D, 5 S.
%
%   The baseline structure is the tree A <-> P <-> D <-> S with R <-> P.
%   Two optional fields, kar (A -> R) and kra (R -> A), add the direct
%   switch that a tug of war resolved inside one frame would produce.
%   They close a cycle A -> R -> P -> A, so detailed balance is lost and
%   the occupancies must come from the null space rather than from the
%   product of ratios along the tree.  Omitting them, or setting them to
%   zero, recovers the baseline generator exactly.

Q = zeros(5,5);
Q(1,3) = p.kap;  Q(3,1) = p.kpa;
Q(2,3) = p.krp;  Q(3,2) = p.kpr;
Q(3,4) = p.kpd;  Q(4,3) = p.kdp;
Q(4,5) = p.kds;  Q(5,4) = p.ksd;

if isfield(p,'kar') && p.kar > 0, Q(1,2) = p.kar; end
if isfield(p,'kra') && p.kra > 0, Q(2,1) = p.kra; end

Q = Q - diag(sum(Q,2));
end
