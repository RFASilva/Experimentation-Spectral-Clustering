function [validdata, index] = remove_nans(data);
%
%

[n d] = size(data);
nvalid = isnan(data);
index = (~(sum(nvalid')>0))';
validdata = data(index,:);

%fprintf('\n%d rows with NaNs removed out of %d rows in total\n\n',sum(~index),n);