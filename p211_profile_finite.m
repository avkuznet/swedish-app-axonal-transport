function out = p211_profile_finite(p, opt)
%P211_PROFILE_FINITE  Steady five-state profile on a finite axon 0 <= x <= Lax.
%
%   out = p211_profile_finite(p, opt)
%
%   Same steady equations as p211_profile (Eqs. 3-7 of the manuscript),
%   with P, D and S eliminated algebraically, so that
%
%       d/dx [a; r] = B [a; r],
%
%   but on a finite axon instead of a semi-infinite one.  Both spatial
%   modes are then retained,
%
%       [a; r](x) = c1 w1 exp(mu1 x) + c2 w2 exp(mu2 (x - Lax)),
%
%   mu1 < 0 decaying from the soma, mu2 > 0 decaying from the terminal
%   (written relative to x = Lax so that nothing overflows), and the two
%   constants are fixed by one condition at each end:
%
%     soma      x = 0    a(0) = 1      (anterograde supply; the retrograde
%                                       flux arriving at the soma leaves)
%     terminal  x = Lax  vr r(Lax) = rho va a(Lax) + jin va
%
%   rho is the fraction of the anterograde flux arriving at the terminal
%   that is sent back retrogradely (0: absorbing terminal, 1: reflecting,
%   zero net flux), and jin is an additional retrograde supply from the
%   terminal, in units of the anterograde flux entering at the soma.
%   Lax -> Inf with rho = jin = 0 recovers p211_profile exactly.
%
%   opt.lam   loss rate, 1/s (default: in-axon half-life 3.63 h)
%   opt.Lax   axon length, um (default 1000)
%   opt.rho   terminal return fraction (default 0)
%   opt.jin   distal retrograde supply (default 0)
%   opt.L     imaged window, um (default 200);  opt.f0 (0.78); opt.nx (2001)
%
%   out.J is the net number flux va a - vr r; out.balance checks the
%   conservation identity dJ/dx = -lam C (relative error).

if nargin < 2, opt = struct; end
if ~isfield(opt,'lam'), opt.lam = log(2)/(3.63*3600); end
if ~isfield(opt,'Lax'), opt.Lax = 1000; end
if ~isfield(opt,'rho'), opt.rho = 0;    end
if ~isfield(opt,'jin'), opt.jin = 0;    end
if ~isfield(opt,'L'),   opt.L   = 200;  end
if ~isfield(opt,'f0'),  opt.f0  = 0.78; end
if ~isfield(opt,'nx'),  opt.nx  = 2001; end
lam = opt.lam;  Lax = opt.Lax;

M = (p.kdp + p.kds + lam) - p.ksd*p.kds/(p.ksd + lam);
G = (p.kpa + p.kpr + p.kpd + lam) - p.kdp*p.kpd/M;
B = [ (-(p.kap+lam) + p.kpa*p.kap/G)/p.va ,  p.kpa*p.krp/G/p.va ; ...
      (-p.kpr*p.kap/G)/p.vr               , ((p.krp+lam) - p.kpr*p.krp/G)/p.vr ];
[V, Dg] = eig(B);  mu = real(diag(Dg));  V = real(V);
[mu, k] = sort(mu);  V = V(:,k);           % mu(1) < 0 < mu(2)
w1 = V(:,1)/V(1,1);  w2 = V(:,2);

if isinf(Lax)
    c = [1; 0];
else
    % a(0) = 1:  c1 + c2 w2(1) e^{-mu2 Lax} = 1
    % terminal:  vr r(L) - rho va a(L) = jin va
    E1 = exp(mu(1)*Lax);  E2 = exp(-mu(2)*Lax);
    tL = @(w) p.vr*w(2) - opt.rho*p.va*w(1);
    Mx = [ 1,            w2(1)*E2 ;
           tL(w1)*E1,    tL(w2)   ];
    c = Mx \ [1; opt.jin*p.va];
end

x = linspace(0, min(opt.L, Lax), opt.nx)';
m1 = exp(mu(1)*x);
if isinf(Lax), m2 = zeros(size(x)); else, m2 = exp(mu(2)*(x - Lax)); end
a = c(1)*w1(1)*m1 + c(2)*w2(1)*m2;
r = c(1)*w1(2)*m1 + c(2)*w2(2)*m2;
pp = (p.kap*a + p.krp*r)/G;
dd = p.kpd*pp/M;
ss = p.kds*dd/(p.ksd + lam);
C  = a + r + pp + dd + ss;
J  = p.va*a - p.vr*r;
f  = opt.f0 * C/C(1);

iP = x <= 100;  iD = x >= 100;
out.x = x;  out.a = a;  out.r = r;  out.C = C;  out.J = J;  out.f = f;
out.aucProx = trapz(x(iP), f(iP));
out.aucDist = trapz(x(iD), f(iD));
out.ellWin  = 100/log(out.aucProx/out.aucDist);   % length implied by the two windows
out.mu = mu;  out.ell1 = -1/mu(1);  out.ell2 = 1/mu(2);  out.c = c;
dJ = gradient(J, x);
out.balance = max(abs(dJ + lam*C))/max(abs(lam*C));
end
