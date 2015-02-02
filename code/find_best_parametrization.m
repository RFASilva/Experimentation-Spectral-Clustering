% In this file there are several scripts in order to run the algorithms
% with several parameters in order to find the optimal adjusted rand index
% value. (basically this is just a log file of what was done).


% Find best parametrization

%% Gaussian Equal Size: Three Clusters
%size = [200 200 200];
%best_params = zeros(1000,1);
%[data, label] = gaussian_data_generator(3,1,size);
%exec = 1;

%info = zeros(1000,2);
%for k=1:80
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, 3, 19, k, 1, 1);
%    
%    %% Compute Adjusted Rand Index
%    clusters_gen_STD3 = zeros(600,1);
%
%    for i=1:length(clusts_STD3)
%        clusters_gen_STD3(clusts_STD3{i}) = i;
%    end
%    ari = RandIndex(clusters_gen_STD3, label)
%    best_params(k, 1) = ari;
%    
%    exec = exec + 1;
%end

%for k=2:2:20
%        for p=3:8:80
%            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, 17, 17, 19, k, p);

                %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(600,1);
%                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
%               
%                ari = RandIndex(label, clusters_gen_STD4);
%                
%                best_params(exec,1) = ari;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

%% Gaussian Different Size: Three Clusters

%size = [50 200 100];
%best_params = zeros(1000,1);
%[data, label] = gaussian_data_generator(3,1,size);
%exec = 1;

%info = zeros(1000,2);


%for k=1:2:80
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, 3,7, k, 1, 1);
%    
    %% Compute Adjusted Rand Index
%    clusters_gen_STD3 = zeros(350,1);

%    for i=1:length(clusts_STD3)
%        clusters_gen_STD3(clusts_STD3{i}) = i;
%    end
%    ari = RandIndex(clusters_gen_STD3, label);
%    best_params(k, 1) = ari;
%    
%    exec = exec + 1;
%end

% size = [50 200 100];
% best_params = zeros(1000,1);
% [data, label] = gaussian_data_generator(3,1,size);
% exec = 1;

% info = zeros(1000,2);

% for k=2:2:20
%        for p=3:8:80
            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, 17, 17, 19, k, p);

               %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(350,1);
%                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
%               
%                ari = RandIndex(label, clusters_gen_STD4);
%                
%                best_params(exec,1) = ari;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end


%% Gaussian Equal Size: Four Clusters
%size = [200 200 200 200];
%best_params = zeros(1000,1);
%[data, label] = gaussian_data_generator(4,1,size);
%exec = 1;

%info = zeros(1000,2);
%for k=1:80
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,4,label, 3, 28, k, 1, 1);
    
    %% Compute Adjusted Rand Index
%    clusters_gen_STD3 = zeros(800,1);

%    for i=1:length(clusts_STD3)
%        clusters_gen_STD3(clusts_STD3{i}) = i;
%    end
%    ari = RandIndex(clusters_gen_STD3, label)
%    best_params(k, 1) = ari;
    
%    exec = exec + 1;
%end


%bsize = [200 200 200 200];
%best_params = zeros(1000,1);
%[data, label] = gaussian_data_generator(4,1,size);
%exec = 1;

%info = zeros(1000,2);

%for k=2:2:20
%        for p=3:8:80
%            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,4,label, 17, 17, 19, k, p);
                %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(800,1);
                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
               
%                ari = RandIndex(label, clusters_gen_STD4);
                
%                best_params(exec,1) = ari;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

% size = [50 150 200 100];
% best_params = zeros(1000,1);
% [data, label] = gaussian_data_generator(4,1,size);
% exec = 1;

% info = zeros(1000,2);

% for k=2:2:20
%        for p=3:8:80
%            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,4,label, 17, 17, 19, k, p);
                %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(500,1);
                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
%               
%                ari = RandIndex(label, clusters_gen_STD4);
                
%                best_params(exec,1) = ari;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

% size = [50 150 200 100];
% best_params = zeros(1000,1);
% [data, label] = gaussian_data_generator(4,1,size);
% exec = 1;

% info = zeros(1000,2);
 
%for k=1:80
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,4,label, k, 44, 2, 1, 1);
%    
%    %% Compute Adjusted Rand Index
%    clusters_gen_STD1 = zeros(500,1);
    
%    for i=1:length(clusts_STD1)
%        clusters_gen_STD1(clusts_STD1{i}) = i;
%    end
    
%    ari = RandIndex(clusters_gen_STD1, label)
%    best_params(k, 1) = ari;
    
%    exec = exec + 1;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%% DATA 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find best params for data0, data1, data2 etc...
%best_params1 = zeros(1000,3);
%exec = 1;
 
%for k=1:3:70
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data0,3,group0, k, k, k, 1, 1);
%    
%    %% Compute Adjusted Rand Index
%    clusters_gen_STD1 = zeros(318,1);
%    clusters_gen_STD2 = zeros(318,1);
%    clusters_gen_STD3 = zeros(318,1);
%   
%   
%    for i=1:length(clusts_STD1)
%        clusters_gen_STD1(clusts_STD1{i}) = i;
%    end
%    for i=1:length(clusts_STD2)
%        clusters_gen_STD2(clusts_STD2{i}) = i;
%    end
%    for i=1:length(clusts_STD3)
%        clusters_gen_STD3(clusts_STD3{i}) = i;
%    end
%   
%    ari1 = RandIndex(clusters_gen_STD1, group0);
%    ari2 = RandIndex(clusters_gen_STD2, group0);
%    ari3 = RandIndex(clusters_gen_STD3, group0);
%    
%    best_params1(k, :) = [ari1 ari2 ari3];
%    
%    exec = exec + 1;
%end

% best_params2 = zeros(1000,1);
% exec = 1;
% info2 = zeros(1000,2);

% for k=2:2:40
%        for p=3:4:60
%            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data0,3,group0, 1, 1, 1, k, p);
                %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(303,1);
                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
               
%                ari = RandIndex(group4, clusters_gen_STD4);
               
%                best_params2(exec,1) = ari;
%                info2(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%% DATA FINAL 2 - FUZZY SCHEME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %best_params3 = zeros(1500,1);
 %exec = 1;
 %info3 = zeros(1500,2);

 %for k=2:2:80
 %       for p=3:5:150
 %           
 %           if(p>=k)
 %               [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data_final_2,5,group_final_2, 1, 1, 1, k, p);
 %               %% Compute Adjusted Rand Index
 %               clusters_gen_STD4 = zeros(770,1);
 %               
 %               for i=1:length(clusts_STD4)
 %                   clusters_gen_STD4(clusts_STD4{i}) = i;
 %               end
 %              
 %               ari = RandIndex(group_final_2, clusters_gen_STD4);
 %              
 %               best_params3(exec,1) = ari;
 %               info3(exec,:) = [k p];
 %               exec = exec + 1;
 %           end
 %       end
%end




%% Gaussian Equal Size: Three Clusters
%size = [200 200 200];
%best_params = zeros(1000,1);
%[data, label] = gaussian_data_generator(3,0.1,size);
%exec = 1;

% info = zeros(1500,2);
%for k=2:2:10
%        for p=3:3:40
            
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, 1, 1, 1, k, p);
%                %% Compute Adjusted Rand Index
%                clusters_gen_STD4 = zeros(600,1);
%                
%                for i=1:length(clusts_STD4)
%                    clusters_gen_STD4(clusts_STD4{i}) = i;
%                end
%               
%                ari = RandIndex(label, clusters_gen_STD4);
%                
%                best_params(exec,1) = ari;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

%% Gaussian Equal Size: Three Clusters
%size = [200 200 200];
%best_params3 = zeros(1000,3);
%[data, label] = gaussian_data_generator(3,0.1,size);
%exec = 1;

%for k=1:60
%    [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(data,3,label, k, k, k, 2, 3);
    
   %% Compute Adjusted Rand Index
%    clusters_gen_STD1 = zeros(600,1);
%    clusters_gen_STD2 = zeros(600,1);
%    clusters_gen_STD3 = zeros(600,1);
   
%    for i=1:length(clusts_STD1)
%        clusters_gen_STD1(clusts_STD1{i}) = i;
%    end
%    for i=1:length(clusts_STD2)
%        clusters_gen_STD2(clusts_STD2{i}) = i;
%   end
%    for i=1:length(clusts_STD3)
%        clusters_gen_STD3(clusts_STD3{i}) = i;
%    end
   
%    ari1 = RandIndex(clusters_gen_STD1, label);
%    ari2 = RandIndex(clusters_gen_STD2, label);
%    ari3 = RandIndex(clusters_gen_STD3, label);
    
%    best_params3(k, :) = [ari1 ari2 ari3];
    
%    exec = exec + 1;
%end



%% Breast Cancer Wisconsin (Diagnostic) Data Set 

%best_params = zeros(1000,3);
%exec = 1;

%for k=1:100
%   [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(wdbc,2,wdbc_labels, k, k, k, 1, 1);
    
   %% Compute Adjusted Rand Index
%    ari1 = eva_ari(clusts_STD1, wdbc_labels);
%    ari2 = eva_ari(clusts_STD2, wdbc_labels);
%    ari3 = eva_ari(clusts_STD3, wdbc_labels);
%    
%    best_params(k, :) = [ari1 ari2 ari3];
   
%end

%best_params = zeros(1000,1);
%exec=60;
%info = zeros(1500,2);
 
%for k=2:2:10
%        for p=41:3:60
           
%            if(p>=k)
%                [clusts_STD1, clusts_STD2, clusts_STD3, clusts_STD4] = run_njw(wdbc,2,wdbc_labels, 1, 1, 1, k, p);
                %% Compute Adjusted Rand Index
%                ari4 = eva_ari(clusts_STD4, wdbc_labels);
                
%                best_params(exec,1) = ari4;
%                info(exec,:) = [k p];
%                exec = exec + 1;
%            end
%        end
%end

