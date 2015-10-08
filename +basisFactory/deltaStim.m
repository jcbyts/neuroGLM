function stim = deltaStim(bt, nT, v)
% Returns a sparse vector with events at binned timings

o = ones(numel(bt), 1);

if nargin < 3
    v = o;
end

assert(numel(o) == numel(v));

if any(bt>nT)
    badix=bt>nT;
%     fprintf('%d events of %d are thrown out because they exceed the trial duration\n',sum(badix), numel(badix))
    bt(badix)=[];
    o(badix)=[];
    v(badix)=[];
end
stim = sparse(bt, o, v, nT, 1);