
function [V,ss,L, numberClusters] = evecs(A, flagGroup, nGroups, flagHeuristic)

%% Calculate eigenvectors, eigenvalues of the laplaican of A
%%
%%   [V,ss,L] = evecs(A,nEvecs)
%%  
%%  Input:
%%        A = Affinity matrix
%%        flagGroup = Tell us if the number of the groups is given by input
%%                    or not 
%%        nGroups = number of eigenvectors to compute if we dont find
%%        automatically the number of groups
%%	    flagHeuristic = if flag=1 use eigengap heuristic, if flag=2 use "eigenvalues with value one" heuristic	
%%        
%%  Output:       
%%        V = eigenvectors
%%        ss = eigenvalues
%%        L = Laplacian
%%
%%  Code by Lihi Zelnik-Manor (2005)
%%  Modified by Ricardo Silva

%% Compute the Laplacian
tic;
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
ttt = toc;
 disp(['Laplacian computation took ' num2str(ttt) ' seconds']);

if (useSparse)
    opts.issym = 1;
    opts.isreal = 1;
    opts.disp = 0;
    [V,ss] = eigs(L,nGroups,1,opts);
    
    if flagGroup == 1
        numberClusters = nGroups;
    else
         if flagHeuristic == 1 
            numberClusters = eigengap1(ss, 1, 0, -1);
         end
         if flagHeuristic == 2 
            numberClusters = eigengap2(ss, 1, 0, -1);
         end
    end
%     [VV,ss]=svds(L,nClusts,1,opts);
else
%%%%%%% Compute eigenvectors
tic;
[V,ss] = svd(L);
    if flagGroup == 1
        numberClusters = nGroups;
    else
         if flagHeuristic == 1
            numberClusters = eigengap1(ss, 1, 0, -1);
         end
         if flagHeuristic == 2
            numberClusters = eigengap2(ss, 1, 0, -1);
         end
    end
    V = V(:,1:numberClusters);
ttt = toc;
disp(['eigenvectors computation took ' num2str(ttt) ' seconds']);
end

