function [C,U] = gcut_fuzzy(A, nGroups)

%% [clusts,distortion]=gcut(A,nClusts)
%%
%% Graph partitioning using fuzzy spectral clustering. 
%% Input: 
%%   A = Affinity matrix
%%   nClusts = number of clusters
%% 
%% Output:
%%   clusts = a cell array with indices for each cluster
%%   distortion = the distortion of the final clustering 
%%
%% Algorithm steps:
%% 1. Obtain Laplacian of the affinity matrix
%% 2. Compute eigenvectors of Laplacian
%% 3. Normalize the rows of the eigenvectors
%% 4. Fuzzy C-Means on the rows of the normalized eigenvectors
%%
%% Original code by Yair Weiss
%% Modified and Updated by Lihi Zelnik-Manor and Ricardo Silva

%%%%%%%% Compute the Laplacian
npix = size(A,1);
useSparse = issparse(A);

dd = 1./(sum(A)+eps);
dd = sqrt(dd);
if(useSparse)
    DD = sparse(1:npix,1:npix,dd);
else
    DD = diag(dd);
end
L = DD*A*DD;

%%%%%%% Compute eigenvectors
if (useSparse)
    opts.issym = 1;
    opts.isreal = 1;
    opts.disp = 0;
    [V,ss] = eigs(L,nClusts,1,opts);
    
%     [VV,ss]=svds(L,nClusts,1,opts);
else
    [V,ss] = svd(L);
    V = V(:,1:nGroups);  
end

%%%%%%% Normalize rows of V
for i=1:size(V,1);
    V(i,:)=V(i,:)/(norm(V(i,:))+1e-10);
end

%%  Fuzzy C-Means

[C, U] = fcm(V,nGroups,2);
