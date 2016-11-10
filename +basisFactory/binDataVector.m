function [binnedData, binMat] = binDataVector(dataVector, binSize)
% binnedData = binDataVector(dataVector, binSize)
% bin dataVector into bins of size binSize
% (c) jly 2013
% if binSize==1
%     binnedData = dataVector(:);
%     return
% end
% 
% a = size(dataVector);
% if a(2)>a(1)
%     dataVector = dataVector';
%     a = fliplr(a);
% end

% B = sparse(a(1),a(1));
% 
% for t = 1:a(1)
%     B(t:(t+binSize-1),t) = 1;
% end
% 
% A = B(1:a(1),1:binSize:end);
% 
% binnedData = (dataVector'*A)/binSize;
% binnedData = binnedData(:);
% nd=numel(dataVector);
% nBins=ceil(nd/binSize);
X=filter(ones(binSize,1)/binSize, 1, full(dataVector));
binnedData=X(1:binSize:end);
if nargout==2
    nB=numel(dataVector);
    binMat=filter(ones(binSize,1), 1, full(sparse(1:binSize:nB, 1:ceil(nB/binSize), 1, nB, ceil(nB/binSize))));
end