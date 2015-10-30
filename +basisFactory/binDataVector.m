function [binnedData,A] = binDataVector(dataVector, binSize)
% binnedData = binDataVector(dataVector, binSize)
% bin dataVector into bins of size binSize
% (c) jly 2013
if binSize==1
    binnedData = dataVector(:);
    return
end

a = size(dataVector);
if a(2)>a(1)
    dataVector = dataVector';
    a = fliplr(a);
end

B = sparse(a(1),a(1));

for t = 1:a(1)
    B(t:(t+binSize-1),t) = 1;
end

A = B(1:a(1),1:binSize:end);

binnedData = (dataVector'*A)/binSize;
binnedData = binnedData(:);