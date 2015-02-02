function [data,labels] = pre_process_multivariate()

allData = importdata('wdbc.data');

%% Read Data
data = allData.data;
labels = allData.textdata(:,2);

%% Remove Nans
[n d] = size(data);
nvalid = isnan(data);
index = (~(sum(nvalid')>0))';
data = data(index,:);
labels = labels(index,:);

%% Normalize Data
l = length(data);
min_data = repmat(min(data), l, 1);
max_data = repmat(max(data), l, 1);
data = (data - min_data) ./ (max_data - min_data);

%% Mutltivariate data -- PCA (Principle Component Analysis)
[pc, scores] = princomp(data);
datad = data*pc(:,1:2);

%% Data Visualization
gscatter(datad(:,1),datad(:,2), labels, '', '.');

end