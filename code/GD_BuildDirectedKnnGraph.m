function W = GD_BuildDirectedKnnGraph(M,k,which_matrix)
% Usage: W = GD_BuildDirectedKnnGraph(M,k,which_matrix) 
%
% Input: 
% M                = either the distance or the similarity matrix, needs to be square, symmetric, non-negative
% k                = connectivity parameter of the kNN graph
% which_matrix     = either 'sim' or 'dist' (similarity or distance matrix)
% 
% Output: 
% W              = adjacency matrix of the directed knn graph
%
% For a similarity matrix S, returns the directed knn graph, edges are weighted by S. 
% For a distance matrix D, returns the directed (unweighted!) knn graph. If you want to get 
% a weighted graph in this case, you need to take care of transforming D to S yourself and then 
% call the function with a similarity matrix.  
% Self-edges are excluded in both graphs. 

% implemented by brute force sorting

% testing whether matrix is square: 
if (size(M,1) ~= size(M,2))
  error('Matrix not square!')
end
n = size(M,1); 

% to exclude self-edges, set diagonal of sim/dissim matrix to Inf or 0
for it=1:n
  if (strcmp(which_matrix,'sim'))
    M(it,it) = 0; 
  elseif (strcmp(which_matrix, 'dist'))
    M(it,it) = Inf; 
  else
    error('Unknown matrix type')
  end
end


% now do it: 
W = M; 

if (strcmp(which_matrix, 'sim'))
  
  for it = 1:n
    % sort points according to similarities: 
    [sorted_row,order] = sort(M(it,:), 'descend'); 
    
    % for all points which are not among the k nearest neighbors, set W to 0: 
    W(it, order(k+1:end)) = 0;
    %W(it,order(1:k) just stays the same
  end
  
elseif (strcmp(which_matrix, 'dist'))
  
  for it = 1:n
    % sort points according to distances: 
    [sorted_row,order] = sort(M(it,:), 'ascend'); 
    % for all points which are not among the k nearest neighbors, set W to 0: 
    W(it, order(k+1:end)) = 0;
    W(it, order(1:k)) = 1; %unweighted! 
  end
  
else 
  error('build_directed_knn_graph: unknown matrix type')
end


  
