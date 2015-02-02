function ari=eva_ari(cluster_gen,label, fuzzy)

% Compute the adjusted rand index

% Input: 
    % cluster_gen - result of a clustering algorithm
    % label - stores the id cluster for each data point
    % fuzzy - 1 - indicate if the result is a result from the NJW algorithm
    % using fuzzy c-means
    %       - 2 - indicate if the result is a result from the difuzzy
    %       - 3 - otherwise
% Output: 
    % ari - the value of adjusted rand index
    
temp = zeros(length(label),1);

if (fuzzy == 1)
    [Y temp] = max(cluster_gen);
end

if (fuzzy == 2)

    rows = length(cluster_gen);

    %% Difuzzification
    for i=1:rows
        [max_value prot_id] = max(cluster_gen(i,:));
        temp(i) = prot_id;
    end
end
if (fuzzy == 3)
    
    for i=1:length(cluster_gen)
        temp(cluster_gen{i}) = i;
    end
    
end

ari = RandIndex(temp, label);

end

