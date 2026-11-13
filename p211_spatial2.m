function varargout = p211_spatial2(cmd, varargin)
%P211_SPATIAL2  Axonal profile of APP in co-expressing axons.
%
%   Spatial models for the fixed-cell profiles (Sec. 2.5).  The paper uses
%   the 'chains' command only, for the normalised transport profiles of
%   Eq. (10); 'fit' and 'predict' serve run_spatial2.  Two carrier populations are
%   present in an axon that expresses both variants: W, moving with the
%   five-state chain fitted to APPwt movies, and S, moving with the chain
%   fitted to APPswe movies.  Each protein is distributed over both.  With
%   cW and cS the steady profiles of the two chains (p211_profile),
%   normalised to one at the hillock,
%
%       c_wt  = (1 - m_w) cW + m_w cS
%       c_swe =  m_s cW + (1 - m_s) cS
%
%   and the fixed-cell intensity is
%
%       f_g(x) = kappa * A(x) * c_g(x),     A(x) = tau(x)/tau(0)
%
%   where A is the axon's own geometry, reported by tau in the same axons,
%   and kappa the hillock-to-soma intensity ratio.
%
%   B = p211_spatial2('chains', P, lam, xg)
%       normalised profiles B.cW, B.cS on grid xg for the APP pair
%       P{1} (wt), P{2} (swe); decay lengths B.ellW, B.ellS; and the
%       hillock signal per unit exported flux, B.kW, B.kS.
%
%   f = p211_spatial2('fit', name, A, cW, cS, P, xg, id, yd [, u0])
%       least-squares fit of one candidate (see run_spatial2) to data yd
%       sampled at grid indices id; kappa is always fitted.
%
%   [y, fw, fs, theta, txt] = p211_spatial2('predict', name, u, A, cW, cS, P, xg, id)

switch cmd
    case 'chains',  varargout{1} = chains(varargin{:});
    case 'fit',     varargout{1} = fitmodel(varargin{:});
    case 'predict', [varargout{1:5}] = predict(varargin{:});
    otherwise, error('p211_spatial2: unknown command %s', cmd);
end
end

% ------------------------------------------------------------------------
function B = chains(P, lam, xg)
o = struct('lam',lam,'L',xg(end),'nx',numel(xg),'f0',1);
pw = p211_profile(P{1}, o);  ps = p211_profile(P{2}, o);
B.cW = pw.I/pw.I(1);  B.cS = ps.I/ps.I(1);
B.ellW = pw.ell;  B.ellS = ps.ell;
% hillock signal per unit exported flux, I(0)/J(0), J = va a - vr r
B.kW = pw.I(1)/(P{1}.va*pw.a(1) - P{1}.vr*pw.r(1));
B.kS = ps.I(1)/(P{2}.va*ps.a(1) - P{2}.vr*ps.r(1));
end

% ------------------------------------------------------------------------
function f = fitmodel(name, A, cW, cS, P, xg, id, yd, u0)
switch name
    case 'M0',    k = 1;  ustart = {log(0.8)};
    case 'M0lam', k = 2;  ustart = {[log(0.8) log(log(2)/3600)], [log(0.8) log(1e-3)], ...
                                    [log(0.8) log(1e-5)]};
    case 'M1',    k = 2;  ustart = {[log(0.8) 0], [log(0.8) -2], [log(0.8) 2]};
    case 'M2',    k = 3;  ustart = {[log(0.8) 0 0], [log(0.8) -2 -2], [log(0.8) -2 1], ...
                                    [log(0.8) 1 -2]};
    case 'M3',    k = 3;  ustart = {[log(0.8) 0 log(50)], [log(0.8) 1 log(200)], ...
                                    [log(0.8) -1 log(20)]};
    case 'M3g',   k = 4;  ustart = {[log(0.8) 0 0 log(200)], [log(0.8) 1 1 log(500)], ...
                                    [log(0.8) 1 2 log(300)], [log(0.8) 2 1 log(300)], ...
                                    [log(0.8) -1 -1 log(50)]};
    otherwise, error('p211_spatial2: unknown model %s', name);
end
if nargin >= 9 && ~isempty(u0), ustart = {u0}; end
opt = optimset('TolX',1e-9,'TolFun',1e-12,'MaxFunEvals',4e4,'MaxIter',4e4,'Display','off');
best = inf;
for s = 1:numel(ustart)
    obj = @(u) sum((predict(name,u,A,cW,cS,P,xg,id) - yd).^2);
    [u, val] = fminsearch(obj, ustart{s}, opt);
    if val < best, best = val; ub = u; end
end
[yhat, fw, fs, theta, ptxt] = predict(name, ub, A, cW, cS, P, xg, id);
f.k = k;  f.u = ub;  f.rss = best;  f.yhat = yhat;  f.res = yd - yhat;
f.theta = theta;  f.ptxt = ptxt;  f.fw = fw;  f.fs = fs;
iP = xg <= 100;  iD = xg >= 100;
f.auc = [trapz(xg(iP),fw(iP)) trapz(xg(iD),fw(iD)) trapz(xg(iP),fs(iP)) trapz(xg(iD),fs(iD))];
f.ratio = f.auc(4)/f.auc(2);
end

% ------------------------------------------------------------------------
function [y, fw, fs, theta, ptxt] = predict(name, u, A, cW, cS, P, xg, id)
sig = @(z) 1./(1 + exp(-z));
kap = exp(u(1));
switch name
    case 'M0'
        w = cW;  s = cS;  theta = kap;
        ptxt = sprintf('kappa %.3f', kap);
    case 'M0lam'
        lam = exp(u(2));
        B = chains(P, lam, xg);  w = B.cW;  s = B.cS;  theta = [kap lam];
        ptxt = sprintf('kappa %.3f  lam %.3e 1/s (half-life %.1f h)  ell %.0f / %.1f um', ...
                       kap, lam, log(2)/lam/3600, B.ellW, B.ellS);
    case 'M1'
        m = sig(u(2));
        w = (1-m)*cW + m*cS;  s = m*cW + (1-m)*cS;  theta = [kap m];
        ptxt = sprintf('kappa %.3f  m %.3f', kap, m);
    case 'M2'
        mw = sig(u(2));  ms = sig(u(3));
        w = (1-mw)*cW + mw*cS;  s = ms*cW + (1-ms)*cS;  theta = [kap mw ms];
        ptxt = sprintf('kappa %.3f  m_w %.3f  m_s %.3f', kap, mw, ms);
    case 'M3'
        phi = sig(u(2));  lB = exp(u(3));  e = exp(-xg/lB);
        w = (1-phi)*cW + phi*e;  s = (1-phi)*cS + phi*e;  theta = [kap phi lB];
        ptxt = sprintf('kappa %.3f  phi %.3f  ell_B %.1f um', kap, phi, lB);
    case 'M3g'
        % second route with a genotype-specific share: a fluorescent product
        % whose abundance relative to intact APP may differ between variants
        pw = sig(u(2));  ps = sig(u(3));  lB = exp(u(4));  e = exp(-xg/lB);
        w = (1-pw)*cW + pw*e;  s = (1-ps)*cS + ps*e;  theta = [kap pw ps lB];
        ptxt = sprintf('kappa %.3f  phi_w %.3f  phi_s %.3f  ell_B %.1f um', kap, pw, ps, lB);
end
fw = kap*A.*w;  fs = kap*A.*s;
y = [fw(id); fs(id)];
end
