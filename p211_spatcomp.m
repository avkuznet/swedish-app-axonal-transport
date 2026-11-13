function c = p211_spatcomp(p, lam)
%P211_SPATCOMP  Local state composition of the open, turning-over population.
%
%   c = p211_spatcomp(p, lam) returns the fractions [A R P D S] of the
%   steady spatial solution of Eqs. (3)-(7) with first-order loss lam (1/s)
%   in every state.  On a uniform axon that solution is a single decaying
%   mode, so its composition does not depend on position.  It differs from
%   the stationary distribution of the switching chain (p211_occ) wherever
%   a transition is not fast compared with lam; for instance
%       s/d = kds/(ksd + lam)   instead of   piS/piD = kds/ksd.
%   A movie recorded at any position samples this composition, which is
%   what makes the movie and the spatial calculation one population.

M = (p.kdp + p.kds + lam) - p.ksd*p.kds/(p.ksd + lam);
G = (p.kpa + p.kpr + p.kpd + lam) - p.kdp*p.kpd/M;
B = [ (-(p.kap+lam) + p.kpa*p.kap/G)/p.va ,  p.kpa*p.krp/G/p.va ; ...
      (-p.kpr*p.kap/G)/p.vr               , ((p.krp+lam) - p.kpr*p.krp/G)/p.vr ];
[V, Dg] = eig(B);  [~, i] = min(real(diag(Dg)));
w = real(V(:,i));  w = w/w(1);
a = 1;  r = w(2);
pp = (p.kap*a + p.krp*r)/G;
dd = p.kpd*pp/M;
ss = p.kds*dd/(p.ksd + lam);
c = [a r pp dd ss];  c = c/sum(c);
end
