function out = p211_sim(p, opt)
%P211_SIM  Stochastic five-state cargo trajectories plus observation operator.
%
%   out = p211_sim(p, opt)
%
%   States
%     1  A  anterograde, moving at +v_a       (kinesin engaged)
%     2  R  retrograde,  moving at -v_r       (dynein-dynactin engaged)
%     3  P  paused, still attached to the microtubule
%     4  D  detached from the microtubule
%     5  S  stationary, docked, long-lived
%
%   The allowed transitions form a tree, A <-> P <-> D <-> S with R <-> P,
%   so the chain obeys detailed balance and its occupancies are available
%   in closed form (p211_occ).  A reversal can happen only by way of P: a
%   cargo cannot switch motors without first disengaging the one it has.
%
%   Fields of p (um, s)
%     va, vr       instantaneous speeds while engaged
%     kap, kpa     A -> P, P -> A
%     krp, kpr     R -> P, P -> R
%     kpd, kdp     P -> D, D -> P
%     kds, ksd     D -> S, S -> D
%
%   Fields of opt
%     N      number of tracks                      (default 4000)
%     T      movie duration, s                     (default 30)
%     fps    frame rate, 1/s                       (default 2)
%     dt     integration step, s                   (default 0.02)
%     vthr   velocity cut for the published histograms, um/s
%     vpause frame speed below which a frame counts as paused (0.1 um/s)
%     gap    tracker linking gap, in frames (default 0, i.e. off)
%     tmin   minimum time in motion, s, for a moving track to contribute
%            segmental statistics (default 0, i.e. off)
%     seed   RNG seed, fixed so the optimiser sees a smooth objective
%     init   initial-state distribution [A R P D S] (default: stationary
%            distribution of the chain, p211_occ)
%     exhaustive  true: renormalise the whole-track classes over the
%            classified tracks (default false)
%
%   Optional fields of p
%     kar, kra   direct A -> R and R -> A.  Absent or zero gives the tree,
%                which is the only structure the paper uses.
%
%   OBSERVATION OPERATOR.  Nothing here reads the hidden state.
%   Trajectories are sampled at the movie frame rate and every statistic
%   is then computed from the sampled positions exactly as the source
%   paper computes it.  Points (i)-(iii) describe the operator; (iv) and
%   (v) are optional extensions.
%
%   (i)  A pause frequency measured from 0.5 s frames with a displacement
%        threshold is not the microscopic rate kap.  It is what kap looks
%        like through that filter, and the two differ by a factor that
%        depends on every other rate in the chain.
%
%   (ii) The published analysis is hierarchical.  Tracks are first split
%        into stationary and moving; the segmental statistics are then
%        reported for the moving subset only.  Time budgets, track
%        lengths, pause and reversal frequencies and segmental speeds are
%        therefore all CONDITIONAL on a track having been called moving.
%        Applying them to the whole population instead makes the
%        anterograde, retrograde and pause fractions of APPswe (which sum
%        to one over 75% of the time in motion) irreconcilable with its
%        78% stationary tracks.  Conditioned correctly, the two are
%        statements about different subpopulations and sit together
%        without strain.
%
%   (iii) Particle-tracking software links detections across frames in
%        which the particle was not found, up to a maximum gap.  Every
%        excursion the cargo made inside a bridged gap is then replaced by
%        a straight line between its endpoints, so the reported path is
%        shorter than the path actually walked.  opt.gap switches this on:
%        any run of at most opt.gap consecutive sub-threshold frames,
%        flanked by detected motion, is collapsed to its net displacement
%        before the path length is summed.  out.LA0 and out.LR0 always
%        carry the ungapped values, so the bias is reported rather than
%        merely assumed.
%
%   (iv) opt.speedcv (optional): every run draws its own speed multiplier
%        from a lognormal distribution with mean 1 and coefficient of
%        variation speedcv; v_a and v_r are then mean run speeds and
%        Eq. (2) is unchanged.
%
%   (v)  opt.speedmix (optional; off by default; not used in the paper)
%        draws the multiplier from a mixture of normal modes with relative
%        centres C, weights W and relative width S, normalised so that the
%        mean multiplier is one.
%        Fields: aC, aW, aS (anterograde), rC, rW, rS (retrograde).

if nargin < 2, opt = struct; end
if ~isfield(opt,'N'),      opt.N      = 4000;  end
if ~isfield(opt,'T'),      opt.T      = 30;    end
if ~isfield(opt,'fps'),    opt.fps    = 2;     end
if ~isfield(opt,'dt'),     opt.dt     = 0.02;  end
if ~isfield(opt,'vthr'),   opt.vthr   = 0.5;   end
if ~isfield(opt,'vpause'), opt.vpause = 0.1;   end
if ~isfield(opt,'gap'),    opt.gap    = 0;     end
if ~isfield(opt,'tmin'),   opt.tmin   = 0;     end
if ~isfield(opt,'seed'),   opt.seed   = 20261104; end

N = opt.N;  T = opt.T;  dtf = 1/opt.fps;
spf = round(dtf/opt.dt);              % integration steps per movie frame
dt  = dtf/spf;
nf  = round(T/dtf) + 1;
ns  = (nf-1)*spf;

vel = [p.va, -p.vr, 0, 0, 0];

% generator (p211_gen adds the direct A <-> R switch if kar/kra are set)
Q = p211_gen(p);

% exact one-step transition matrix: no time-step error in the state
% sequence, only in where inside a step the jump lands
Pm = expm(Q*dt);
Pm = max(Pm, 0);
Pm = Pm ./ sum(Pm, 2);
cp = cumsum(Pm, 2);
cp(:,end) = 1;

pi0 = p211_occ(p);
% opt.init (optional) replaces the stationary distribution of the chain as
% the distribution of initial states, e.g. the local composition of the
% open, turning-over population of the spatial model (p211_spatcomp).
if isfield(opt,'init') && ~isempty(opt.init)
    pi0 = opt.init(:)'/sum(opt.init);
end
rs  = RandStream('mt19937ar','Seed',opt.seed);

st = 1 + sum(rand(rs,N,1) > cumsum(pi0(:))', 2)';   % 1 x N
st = min(max(st,1),5);
x  = zeros(1,N);

X  = zeros(nf,N);
fi = 1;
occ = zeros(1,5);

% per-run speed multipliers (opt.speedmix); all ones otherwise
useSM = isfield(opt,'speedmix') && ~isempty(opt.speedmix);
% opt.speedcv (optional): every run draws its speed multiplier from a
% lognormal distribution with mean 1 and coefficient of variation speedcv,
% the same in both directions; va and vr are then mean run speeds.
useCV = isfield(opt,'speedcv') && opt.speedcv > 0;
if useCV
    sg = sqrt(log(1 + opt.speedcv^2));
    opt.speedmix = struct('aC',[],'aW',[],'aS',sg,'rC',[],'rW',[],'rS',sg);
    useSM = true;
end
mlt = ones(1,N);
if useSM
    SM = opt.speedmix;
    mlt(st == 1) = drawmult(rs, nnz(st == 1), SM.aC, SM.aW, SM.aS);
    mlt(st == 2) = drawmult(rs, nnz(st == 2), SM.rC, SM.rW, SM.rS);
end

for s = 1:ns
    x = x + vel(st).*mlt*dt;
    occ = occ + accumarray(st(:), 1, [5 1])';
    u = rand(rs,N,1);
    prev = st;
    st = 1 + sum(u > cp(st,:), 2)';
    if useSM
        nA = st == 1 & prev ~= 1;  nR = st == 2 & prev ~= 2;
        if any(nA), mlt(nA) = drawmult(rs, nnz(nA), SM.aC, SM.aW, SM.aS); end
        if any(nR), mlt(nR) = drawmult(rs, nnz(nR), SM.rC, SM.rW, SM.rS); end
    end
    if mod(s,spf) == 0
        fi = fi + 1;
        X(fi,:) = x;
    end
end
occ = occ/sum(occ);

% ------------------------------------------------------- observation ----
dX = diff(X,1,1);          % (nf-1) x N  per-frame displacement
U  = dX/dtf;               % per-frame velocity

% --- level one: split tracks into stationary and moving, as in the SI
spd  = mean(abs(U),1);
vnet = (X(end,:) - X(1,:))/T;
isSta = spd < 0.1;
isMov = ~isSta;
isAnt = isMov & vnet >  0.01;
isRet = isMov & vnet < -0.01;

out.fracA = mean(isAnt);
out.fracS = mean(isSta);
out.fracR = mean(isRet);
% opt.exhaustive: moving tracks whose net velocity lies within +/-0.01 um/s
% fit neither published class; the published proportions sum to one, so
% they are renormalised over the classified tracks.
if isfield(opt,'exhaustive') && opt.exhaustive
    ncl = out.fracA + out.fracS + out.fracR;
    out.fracA = out.fracA/ncl;  out.fracS = out.fracS/ncl;  out.fracR = out.fracR/ncl;
end
out.fracUncl = 1 - mean(isAnt | isSta | isRet);

% --- level two: everything below is conditional on a moving track
if ~any(isMov)
    out.thetaA = NaN; out.thetaR = NaN; out.thetaP = NaN;
    out.LA = NaN; out.LR = NaN; out.LA0 = NaN; out.LR0 = NaN;
    out.fpA = NaN; out.fpR = NaN;
    out.frev = NaN; out.vsA = NaN; out.vsR = NaN; out.fpAll = NaN; out.Lall = NaN;
    out.occ = occ; out.veff = p.va*occ(1) - p.vr*occ(2); out.X = X;
    out.fracKept = NaN;
    return
end

% --- the source's duration rule.  Its Supporting Information states, of
%     the transport-dynamics analysis, that tracks were first divided into
%     stationary or moving and that "all the tracks moving at < 10 s were
%     excluded".  A moving track must therefore have spent at least tmin
%     seconds actually in motion to contribute a segmental statistic.  The
%     cut is applied to the moving subset only, which is where the source
%     applies it; the whole-track proportions reported above come from the
%     net analysis, which the same paragraph treats separately.
if opt.tmin > 0
    tmove = sum(abs(U) > opt.vpause, 1) * dtf;      % seconds in motion
    keep  = isMov & (tmove >= opt.tmin);
    if ~any(keep), keep = isMov; end                % never empty the subset
else
    keep = isMov;
end
out.fracKept = mean(keep)/max(mean(isMov), eps);

Um   = U(:,keep);
movA = Um >  opt.vpause;
movR = Um < -opt.vpause;
pau  = ~movA & ~movR;

out.thetaA = mean(movA(:));
out.thetaR = mean(movR(:));
out.thetaP = mean(pau(:));

% --- track length: total path, which is what the tracking software
%     reports ("the distance traveled by a particle along each vector").
%     plen0 is the raw sum of frame displacements; plen is the same after
%     the tracker's linking gap has been bridged.  With opt.gap = 0 the
%     two are identical and nothing below changes.
plen0 = sum(abs(dX),1);
if opt.gap > 0
    plen = gapclose(dX, opt.gap, opt.vpause*dtf);
else
    plen = plen0;
end
out.LA0 = pickmed(plen0, isAnt & keep);
out.LR0 = pickmed(plen0, isRet & keep);
out.LA = pickmed(plen, isAnt & keep);
out.LR = pickmed(plen, isRet & keep);
out.Lall = pickmed(plen, keep);          % pooled over directions (Fig. 6C)

% --- pause frequency per 10 s: entries into a pause run
pauAll = ~(U > opt.vpause) & ~(U < -opt.vpause);
nP = sum(pauAll(2:end,:) & ~pauAll(1:end-1,:), 1) + pauAll(1,:);
out.fpA   = cmed(nP(isAnt & keep)) * (10/T);
out.fpR   = cmed(nP(isRet & keep)) * (10/T);
out.fpAll = cmed(nP(keep)) * (10/T);

% --- reversal frequency per 10 s: sign changes between successive motile
%     frames, ignoring the pauses that sit between them
sg = zeros(size(U));
sg(U >  opt.vpause) =  1;
sg(U < -opt.vpause) = -1;
last = zeros(1,N);
nrev = zeros(1,N);
for k = 1:size(sg,1)
    sk  = sg(k,:);
    nz  = sk ~= 0;
    chg = nz & (last ~= 0) & (sk ~= last);
    nrev = nrev + chg;
    last(nz) = sk(nz);
end
out.frev = cmed(nrev(keep)) * (10/T);

% --- segmental speed distribution, conditioned exactly as published
sa = Um(movA);   sa = sa(sa  > opt.vthr);
sr = -Um(movR);  sr = sr(sr  > opt.vthr);
out.vsA = medsafe(sa);
out.vsR = medsafe(sr);

% --- hidden-state quantities, for the spatial model
out.occ  = occ;
out.veff = p.va*occ(1) - p.vr*occ(2);     % net drift, um/s
out.X    = X;
end

% ------------------------------------------------------------------------
function plen = gapclose(dX, gap, thr)
%GAPCLOSE  Path length after the tracker has bridged its linking gaps.
%
%   dX  (F x N) per-frame displacements
%   gap  longest run of undetected frames the software will link across
%   thr  displacement below which a frame is treated as a non-detection
%
%   A maximal run of at most `gap` consecutive sub-threshold frames that
%   is flanked on BOTH sides by detected motion is replaced by a single
%   straight step from the frame before the run to the frame after it.
%   The path contributed by that run therefore drops from sum|dx| to
%   |sum dx|, which is the entire effect of gap closing on a total-path
%   statistic.  A run touching either end of the track has no flank and
%   is left alone, exactly as a tracker would leave it.
%
%   The loop runs over frames, not over tracks, so its cost is set by the
%   movie length (60 frames here) and not by the number of trajectories.

[F, N] = size(dX);
sub = abs(dX) < thr;

plen = sum(abs(dX), 1);
corr = zeros(1, N);

runLen = zeros(1, N);       % frames in the run currently open
runAbs = zeros(1, N);       % sum |dx| over it
runNet = zeros(1, N);       % sum  dx  over it
runOK  = false(1, N);       % was it preceded by a detected frame?

for k = 1:F
    s = sub(k,:);
    d = dX(k,:);

    % a detected frame closes whatever run was open and flanks it right
    close = ~s & (runLen > 0);
    ok = close & runOK & (runLen <= gap);
    corr(ok) = corr(ok) + (runAbs(ok) - abs(runNet(ok)));
    runLen(close) = 0; runAbs(close) = 0; runNet(close) = 0;
    runOK(close) = false;

    % a detected frame also licenses the NEXT run to be bridged
    runOK(~s) = true;

    % extend or open a run
    runLen(s) = runLen(s) + 1;
    runAbs(s) = runAbs(s) + abs(d(s));
    runNet(s) = runNet(s) + d(s);
end
% runs still open at the last frame have no right flank: not bridged

plen = plen - corr;
end

% ------------------------------------------------------------------------
function m = pickmed(v, idx)
if any(idx), m = median(v(idx)); else, m = NaN; end
end

function m = medsafe(v)
if isempty(v), m = NaN; else, m = median(v); end
end

function m = cmed(v)
%CMED  Median of integer counts, interpolated inside the bin it falls in.
%   The plain median of small integers can only take values 0, 1, 2, ...,
%   so an objective built on it is piecewise constant and traps a simplex
%   search on a flat step.  Spreading each count k over the unit interval
%   [k-1/2, k+1/2] and reading the 50th percentile off the resulting
%   distribution removes the steps without shifting the statistic: it is
%   the ordinary interpolated median for grouped data.
if isempty(v), m = NaN; return, end
v = v(:);
k = unique(v);
cnt = histc(v, k);                        %#ok<HISTC>
pr  = cnt/sum(cnt);
cum = cumsum(pr);
j   = find(cum >= 0.5, 1, 'first');
below = 0; if j > 1, below = cum(j-1); end
m = (k(j) - 0.5) + (0.5 - below)/pr(j);
end

function m = drawmult(rs, n, C, W, S)
%DRAWMULT  n speed multipliers from a mixture of normal modes.
%   C relative centres, W weights, S relative width; the multipliers are
%   normalised so that their population mean is one, and kept positive.
%   With C empty, lognormal multipliers of mean one and log-sd S.
if isempty(C)
    m = exp(S*randn(rs,1,n) - S^2/2);
    return
end
C = C(:)'/sum(W(:)'.*C(:)');           % mean of the mixture = 1
k = 1 + sum(rand(rs,n,1) > cumsum(W(:)')/sum(W), 2)';
m = C(k) + S*randn(rs,1,n);
m = max(m, 0.05);
end
