function stim = boxcarStim(startBinnedTime, endBinnedTime, nT)
% Returns a boxcar duration stimulus design
if endBinnedTime>nT
    endBinnedTime=nT;
end
idx = startBinnedTime:endBinnedTime;
o = ones(numel(idx), 1);

stim = sparse(idx, o, o, nT, 1);