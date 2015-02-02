function [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data, nGroups, labels, k1, k2, k3, k4, n4)

% Runs the NJW algorithm using different similarities approaches

    % Input: 
        % data - dataset
        % nGroups - the number of clusters
        % labels - stores the id cluster for each data point
        % k1 - k parameter used by the Gaussian Kernel
        % k2 - k parameter used by the k-NN
        % k3 - k parameter used by the mutual k-NN
        % k4 - k parameter used by the fuzzy similarity
        % n4 - t parameter used by the fuzzy similarity (number of
        % prototypes)
        

    % Output: 
        % clusts_STD1 - the clustering obtained using the Gaussian Kernel
        % clusts_STD2 - the clustering obtained using the k-NN
        % clusts_STD3 - the clustering obtained using the mutual k-NN
        % clusts_STD4 - the clustering obtained using the fuzzy similarity


    colors = [1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1; 0,1,1; 0,0,0; 0.8,0.8,0.8; 0.4,0.4,0.4; 0.7,0.7,0.7];

    X = data;
    %% Centralize and scale the data
    X = X - repmat(mean(X),size(X,1),1);
    X = X/max(max(abs(X)));

    %%%%%%%%%%%%%%%%% Build affinity matrices   
    %% Affinity: Fully connected graph by Gaussian Kernel
    D1 = dist2(X,X);              %% Euclidean distance
    %A1 = exp(-D1 /(2*sigma_values(j)^2));
    
    [D_LS,A1,LS] = scale_dist(D1,k1); %% Locally scaled affinity matrix
    clear D_LS; clear LS;
    
    %% Affinity: K-NN graph
    D = sqrt(D1);    
    A2 = GD_BuildSymmetricKnnGraph(D,k2,'dist');
    
    %% Affinity; Mutual K-NN graph
    A3 = GD_BuildMutualKnnGraph(D,k3,'dist');
    
    %% Affinity: K-NN fuzzy similarity measure based on fuzzy c-means (FCM) clustering
    A4 = BuildKnnFuzzyGraph(D,k4,n4);
    
    %% Zero out diagonal
    ZERO_DIAG = ~eye(size(X,1));
    A1 = A1.*ZERO_DIAG;
    A2 = A2.*ZERO_DIAG;
    A3 = A3.*ZERO_DIAG;
    A4 = A4.*ZERO_DIAG;
        
    %% Apply Spectral Clustering: Ng., Jordan and Weiss (NJW) algorithm
    clusts_STD1 = gcut(A1, nGroups);
    clusts_STD2 = gcut(A2, nGroups);
    clusts_STD3 = gcut(A3, nGroups);
    clusts_STD4 = gcut(A4, nGroups);
      
   % figure
    %% Display results
    subplot(1,4,1);
        hold on
        for i=1:length(clusts_STD1),
           plot(X(clusts_STD1{i},1),X(clusts_STD1{i},2),'.','Color',colors(i,:),'MarkerSize',16);
        end
        axis equal;
    subplot(1,4,2);    
    hold on
        for y=1:length(clusts_STD2),
            plot(X(clusts_STD2{y},1),X(clusts_STD2{y},2),'.','Color',colors(y,:),'MarkerSize',16);
        end
        axis equal;
    subplot(1,4,3);    
        hold on
        for w=1:length(clusts_STD3),
           plot(X(clusts_STD3{w},1),X(clusts_STD3{w},2),'.','Color',colors(w,:),'MarkerSize',16);
        end
        axis equal;
    subplot(1,4,4);    
        hold on
        for w=1:length(clusts_STD4),
            plot(X(clusts_STD4{w},1),X(clusts_STD4{w},2),'.','Color',colors(w,:),'MarkerSize',16);
        end
        
        axis equal;
        hold off;
       drawnow;
  end
