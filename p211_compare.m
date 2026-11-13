function t = p211_compare(d, m)
%P211_COMPARE  Line up model output against the published statistics.
%
%   t = p211_compare(d, m) returns a struct with one row per observable:
%   name, measured value, model value, and relative residual.  Targets the
%   source does not report come back NaN and are dropped from the fit.

nm  = {'time antero','time retro','time pause', ...
       'tracks antero','tracks stationary','tracks retro', ...
       'track length antero','track length retro', ...
       'pause freq antero','pause freq retro','reversal freq', ...
       'seg speed antero','seg speed retro'};

obs = [d.thetaA d.thetaR d.thetaP d.fracA d.fracS d.fracR ...
       d.LA d.LR d.fpA d.fpR d.frev d.vsA d.vsR];

mod = [m.thetaA m.thetaR m.thetaP m.fracA m.fracS m.fracR ...
       m.LA m.LR m.fpA m.fpR m.frev m.vsA m.vsR];

% pooled targets (p211_data_pooled): Fig. 6C/6D and 7F report one
% distribution over both directions, compared with the pooled model median
if isfield(d,'Lpool')
    nm  = [nm, {'track length pooled','pause freq pooled'}];
    obs = [obs, d.Lpool, d.fpPool];
    mod = [mod, getf(m,'Lall'), getf(m,'fpAll')];
end

sc  = max(abs(obs), 1e-3);
rel = (mod - obs)./sc;

% A model that cannot produce a statistic the experiment reports is not
% excused from it.  Dropping such a target rewards the failure: a chain
% whose anterograde speed falls just below the velocity cut the source
% applied to its own histogram yields no qualifying frames at all, so
% vsA comes back NaN, the target silently leaves the objective, and the
% vector scores better than one that actually fits.  That is how a fit
% with va = 0.497 um/s against a 0.5 um/s cut came to be preferred.
% Such a target is charged a full 100% relative residual instead.
miss = ~isnan(obs) & isnan(mod);
rel(miss) = 1;

rel(isnan(obs)) = NaN;

t.name = nm; t.obs = obs; t.mod = mod; t.rel = rel;
end

function v = getf(m, f)
if isfield(m, f), v = m.(f); else, v = NaN; end
end
