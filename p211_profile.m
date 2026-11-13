function out = p211_profile(p, opt)
%P211_PROFILE  Steady-state cargo distribution along the axon.
%
%   out = p211_profile(p, opt)
%
%   The mean-field limit of the same five-state chain, with a first-order
%   turnover lam acting on every state, on x >= 0:
%
%     d/dx( va a) = -(kap+lam) a + kpa p
%     d/dx(-vr r) = -(krp+lam) r + kpr p
%        0        = kap a + krp r - (kpa+kpr+kpd+lam) p + kdp d
%        0        = kpd p - (kdp+kds+lam) d + ksd s
%        0        = kds d - (ksd+lam) s
%
%   P, D and S carry no flux, so the last three equations are algebraic
%   and eliminate, leaving a linear 2x2 system in (a, r).  Its matrix has
%   exactly one eigenvalue with negative real part whenever lam > 0; the
%   other diverges.  The 0-200 um window imaged in Fig. 2 is a window of a
%   much longer axon, not its terminal, so the bounded solution is the
%   decaying mode alone.  That choice carries no free parameter: the whole
%   profile is then
%
%       f(x) = f0 * exp(mu1 * x),      ell = -1/mu1
%
%   with ell the transport decay length.  For small lam, ell -> v_eff/lam,
%   where v_eff = va*pi_A - vr*pi_R is the net drift of the five-state
%   chain.  A cargo whose net drift collapses cannot hold a distal
%   population against turnover, whatever its instantaneous speeds are.
%
%   opt.lam  turnover rate, 1/s.  Declared, not fitted.  The default is
%            the in-axon half-life of 3.63 h from a log-linear fit to Table 1
%            of Morin et al. (1993), which gives lam = 5.30e-5 1/s.  It is
%            the only published figure measured on APP inside an axon.
%   opt.L    length of the imaged window, um (default 200)
%   opt.f0   hillock-to-soma intensity ratio, read off Fig. 2C (0.78), the
%            same for both variants, so AUC ratios do not depend on it.

if nargin < 2, opt = struct; end
if ~isfield(opt,'lam'), opt.lam = log(2)/(3.63*3600); end
if ~isfield(opt,'L'),   opt.L   = 200;  end
if ~isfield(opt,'f0'),  opt.f0  = 0.78; end
if ~isfield(opt,'nx'),  opt.nx  = 2001; end

% lam may be a scalar, applied to every state, or a 1x5 vector giving the
% loss rate of each state separately in the order A, R, P, D, S.  The
% paper uses the scalar form only (one loss rate for all five states).
if isscalar(opt.lam)
    lamv = opt.lam*ones(1,5);
else
    lamv = opt.lam(:)';
end
lam = lamv(1);                      % kept for the messages below
lA = lamv(1); lR = lamv(2); lP = lamv(3); lD = lamv(4); lS = lamv(5);
L = opt.L;

% optional direct A <-> R transitions (p.kar, p.kra); the paper does not use them.
% They enter the two flux equations only; P still receives kap*a + krp*r
% and loses kpa+kpr+kpd+lam, so G and M are untouched.
kar = 0; if isfield(p,'kar'), kar = p.kar; end
kra = 0; if isfield(p,'kra'), kra = p.kra; end

M = (p.kdp + p.kds + lD) - p.ksd*p.kds/(p.ksd + lS);
G = (p.kpa + p.kpr + p.kpd + lP) - p.kdp*p.kpd/M;

B = [ (-(p.kap+kar+lA) + p.kpa*p.kap/G)/p.va , ( p.kpa*p.krp/G + kra)/p.va ; ...
      (-p.kpr*p.kap/G - kar)/p.vr            , ((p.krp+kra+lR) - p.kpr*p.krp/G)/p.vr ];

[V, Dg] = eig(B);
mu = diag(Dg);

% keep the bounded (decaying) mode
[~, ib] = min(real(mu));
mu1 = real(mu(ib));
w   = real(V(:,ib));
w   = w/w(1);                 % normalise so a(0) = 1

if mu1 >= 0
    warning('p211_profile:unbounded','no decaying mode; check lam > 0');
end

x = linspace(0, L, opt.nx)';
e = exp(mu1*x);
a = w(1)*e;
r = w(2)*e;

pp = (p.kap*a + p.krp*r)/G;
dd = p.kpd*pp/M;
ss = p.kds*dd/(p.ksd + lS);

I = a + r + pp + dd + ss;
f = opt.f0 * I/I(1);

iP = x <= 100;  iD = x >= 100;
out.x = x; out.a = a; out.r = r; out.p = pp; out.d = dd; out.s = ss;
out.I = I; out.f = f;
out.aucProx = trapz(x(iP), f(iP));
out.aucDist = trapz(x(iD), f(iD));
out.mu   = mu1;
out.ell  = -1/mu1;                                   % decay length, um
out.lam  = lam;
out.lamv = lamv;
pi0 = p211_occ(p);
out.veff = p.va*pi0(1) - p.vr*pi0(2);
out.ellAsym = out.veff/lam;                          % leading behaviour
end
