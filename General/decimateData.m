function [ newData ] = decimateData( data, reductionFactor )
%convertSamplingTime Converts a dataset with a higher sampling factor to
%one with a reduced sampling factor.
%   Detailed explanation goes here
newData = data(:,1:reductionFactor:length(data));
% for i = 1:size(data,2)/reductionFactor
%     index = reductionFactor * (i - 1) + 1;
%     endIndex = index + reductionFactor - 1;
%     newData(:,i) = mean(data(:, index : endIndex),2);
% 
% end
end

