function stim = boxcarStim(startBinnedTime, endBinnedTime, nT, val)
% Returns a boxcar duration stimulus design
% stim = boxcarStim(startBinnedTime, endBinnedTime, nT, val)

if nargin < 4
    val=1;
end

endBinnedTime(endBinnedTime>nT)=nT;

idx=[];
o=[];
v=[];
for k=1:numel(startBinnedTime)
    id = startBinnedTime(k):endBinnedTime(k);
    idx=[idx; id(:)];
    o = [o; ones(numel(id), 1)];
    v = [v; val(k)*ones(numel(id), 1)];
end

[idx, ii]=unique(idx);
stim = sparse(idx, o(ii), v(ii), nT, 1);