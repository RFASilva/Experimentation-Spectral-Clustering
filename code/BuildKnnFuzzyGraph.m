function W = BuildKnnFuzzyGraph(M,k,n)

% Input: 
    % M - dataset
    % k - the number k of the nearest prototypes
    % n - the number of prototypes returned by the FCM

% Output: 
%   W  = affinity matrix of the knn fuzzy graph

[prototypes, membership] = fcm(M,n,2);
[rows, coluns] = size(M);

[U index]  = sort(membership, 'descend');

W = zeros(rows, rows);

for i=1:rows
    for j=1:i      
        nearest_prototype_i = index(:,i);
        nearest_prototype_j = index(:,j);
        
        p_i = nearest_prototype_i(1:k);
        p_j = nearest_prototype_j(1:k);
            
        if(nearest_prototype_i(1) == nearest_prototype_j(1))
            W(i,j) = 1;  
        elseif( length(find(p_i == p_j ) == 1) > 0)
            temp = [U(p_i,i)  U(p_j,j) ];
            W(i,j) = max(max(temp));
        else
            W(i,j) = 0;
        end
        
        W(j,i) = W(i,j);
    end
end



end