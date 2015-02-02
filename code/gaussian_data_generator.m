function [data,label] = gaussian_data_generator(nclusters, sn, size)

% Data generated from a bivariate spherical Gaussian distribution 

% Input: 
    % nclusters - number of clusters pretended
    % sn - indicate the level o intermix between clusters
    % size - indicate the size of each cluster

% Output: 
%   data  = data generated
%   label = stores the id cluster for each data point

center = [0 0];
p = 1.2;
p = p + (p*sn);

% Adding the factor
if nclusters==3
MU1 = [p p];
SIGMA1 = [0.975 0;
          0 0.975];

MU2 = [-p p];
SIGMA2 = [0.975 0;
          0 0.975];

MU3 = [0 -p];
SIGMA3 = [0.975 0;
          0 0.975];

data = [mvnrnd(MU1,SIGMA1,size(1));
        mvnrnd(MU2,SIGMA2,size(2));
        mvnrnd(MU3,SIGMA3,size(3))];
 
label = ones(sum(size), 1);


label(1:size(1)) = 1;
size_i = size(1)+1;
size_f = size(1)+size(2);
label(size_i:size_f) = 2;

size_i = size_f + 1;
size_f = sum(size);
label(size_i:size_f) = 3;

end

if nclusters==4
MU1 = [p p];
SIGMA1 = [0.975 0;
          0 0.975];

MU2 = [-p p];
SIGMA2 = [0.975 0;
          0 0.975];

MU3 = [-p -p];
SIGMA3 = [0.975 0;
          0 0.975];

MU4 = [p -p];
SIGMA4 = [0.975 0;
          0 0.975];

data = [mvnrnd(MU1,SIGMA1,size(1));
        mvnrnd(MU2,SIGMA2,size(2));
        mvnrnd(MU3,SIGMA3,size(3));
        mvnrnd(MU4,SIGMA4,size(4))];

label = ones(sum(size), 1);


label(1:size(1)) = 1;
size_i = size(1)+1;
size_f = size(1)+size(2);
label(size_i:size_f) = 2;

size_i = size_f + 1;
size_f = size(1)+size(2)+size(3);
label(size_i:size_f) = 3;

size_i = size_f + 1;
size_f = sum(size);
label(size_i:size_f) = 4;
 
end

end