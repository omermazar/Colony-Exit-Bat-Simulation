function [cl_data] = cluster_by_dist(data, max_dist, it_max_num)


% function [c_data, c_counts, c_width, c_data_ind] = cluster_by_dist(data, max_dist)
% [c_data, c_counts, c__width] = cluster_by_dist(data, max_dist)
%
% the function cluster the data data by max_dist parameter, with priority
% to first signals
%
% inputs:   data : 1d vector of samples
%            max_dist: the maximum distance allowed in each cluster
%           it_max_num: number of 'forward' iterations
%           i
% Outputs:  c_data = the cluster
%           c_counts: how many samples in each cluster
%           c_width: the distance between the max sample and the min sample in the cluster

% sort the data
% [s_data, isx] = sort(data);

% cont data- unique apearncess
data_counts = nonzeros(accumarray(data',data,[],@numel));
u_data = unique(data); % u_idx - the index to go bak to original data 

% diff_data = diff(u_data);
% [s_counts, s_idx] = sort(data_counts,'descend');

%% True Cluster
remain_data = u_data; % sorted by counts
remain_counts = data_counts';

c_num = numel(u_data);
cl_data(c_num) = struct( ...
    'data' , [], ...
    'counts', [], ...
    'width',[], ...
    'data_ind',[], ...
    'all_ind', []);

c_data = zeros(1,c_num);
c_counts = zeros(1,c_num);
c_width = -1*ones(1,c_num);
% c_data_ind = zeros(c_num, numel(data));

n=1;
while ~isempty(remain_data)
    %     [~, max_ind] = max(remain_counts);
    if numel(remain_data) == 1
        c_first = remain_data(1);
    else
        [~, c_first_ind] = min(diff(remain_data));
        c_first =  remain_data(c_first_ind);
    end %if numel(remain_data) == 1
    c_ind = abs(remain_data - c_first) <= max_dist;
    c_counts = sum(remain_counts(c_ind));
    c_mean = sum(remain_counts(c_ind).*remain_data(c_ind)) / c_counts;
    c_width = max(remain_data(c_ind)) - min(remain_data(c_ind)); 
    
    % forward iteration around the average
    it_flag =1; % condition
    n_it =1;
    
    while it_flag
        it_flag = 0;    % if condition not met - stop iteration
        it_ind = abs(remain_data - c_mean) <= max_dist;
        it_width = max(remain_data(it_ind)) - min(remain_data(it_ind)); 
        
        % iteratrion condition to update data and continue
        if sum(it_ind) > sum(c_ind) ...    
               || ( sum(it_ind)==sum(c_ind) && it_width<c_width )
            % update data
            c_ind = it_ind; 
            c_counts = sum(remain_counts(it_ind ));
            c_mean = sum(remain_counts(c_ind).*remain_data(c_ind)) / c_counts; % 
            c_width = it_width;
            it_flag= 1;
        end % if
        
        % max iteration stopm
        n_it = n_it +1;
        if n_it > it_max_num
             it_flag= 0;
        end %if n_it > = it_max_num
    end %while it_flag
    
%     cl_data(n).counts = sum(remain_counts(c_ind));
    cl_data(n).counts = c_counts;
    cl_data(n).data_ind = dsearchn(data', remain_data(c_ind)');
    cl_data(n).width = c_width; 
    % average in cluster
    cl_data(n).data = c_mean;
    
    % index of full data XXX
 
     cl_data(n).all_ind = find(abs(data - cl_data(n).data) <= max_dist );

    
    % update remain groups
    remain_idx = ~c_ind;
    remain_data = remain_data(remain_idx);
    remain_counts = remain_counts(remain_idx);
    n = n+1;
end % while ~isempty(remain_data)

notempty_idx = find(~cellfun(@isempty,{cl_data.data}) );
cl_data= cl_data(notempty_idx);

%