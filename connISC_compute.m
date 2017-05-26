function [mean_conn,edgelist,conn_mat,mean_conn_internal,conn_mat_internal,isc_inds] = connISC_compute(roits,blacklist)
%CONNICS_COMPUTE Summary of this function goes here
%   Detailed explanation goes here

N_subj = length(roits);
N_roi = size(roits{1},2);
if nargin<2
    blacklist = [];
end

for i=1:N_subj
    if abs(1-var(roits{i}(:,1)))>1e-6 || abs(mean(roits{i}(:,1)))>1e-6
        error('Data is not z-scored!')
    end        
end

N_subj_pairs = N_subj^2 - N_subj;
N_roi_pairs = N_roi*(N_roi-1)/2;

conn_mat = zeros(N_roi_pairs,N_subj_pairs,'single');
conn_mat_internal = zeros(N_roi_pairs,N_subj,'single');
T = size(roits{1},1);
div = T-1;
edgelist = zeros(N_roi_pairs,2);
kk=0;
good_edges = true(N_roi_pairs,1);
for s = 1:N_roi
    for t = (s+1):N_roi       
        kk=kk+1;
        edgelist(kk,1)=s;
        edgelist(kk,2)=t;
        if ~isempty(blacklist)        
            if any(ismember([s,t],blacklist))
                good_edges(kk)=false;
            end
        end
    end
end

fprintf('Computing connectivity data (internal)\n')
for i=1:N_subj
    fprintf('...subj. %i/%i\n',i,N_subj);
    kk=0;
    for s = 1:N_roi
        a = roits{i}(:,s);
        for t = (s+1):N_roi
            b = roits{i}(:,t);
            kk=kk+1;
            conn_mat_internal(kk,i)=sum(sum(a.*b));
        end
    end    
    if max(abs(conn_mat_internal(good_edges,i))/div)>0.99
        warning('Very high correlations found (subject %i), possible bug!',i)
    end    
end
conn_mat_internal = conn_mat_internal/div;
mean_conn_internal = mean(atanh(conn_mat_internal),2);

fprintf('Computing connectivity data (isc)\n')
isc_inds=[];
k=0;
for i=1:N_subj    
    fprintf('...source subj. %i/%i\n',i,N_subj);
    for j=1:N_subj
        if i~=j
            k=k+1;
            isc_inds(k,:)=[i,j];
            kk=0;
            for s = 1:N_roi
                a = roits{i}(:,s);
                for t = (s+1):N_roi
                    b = roits{j}(:,t);
                    kk=kk+1;
                    conn_mat(kk,k)=sum(sum(a.*b));
                end
            end                    
        end
        %nnz(good_edges)
        if k>0 && max(abs(conn_mat(good_edges,k))/div)>0.99
            warning('Very high correlations found (pair %i), possible bug!',k)
        end
    end
end
conn_mat = conn_mat/div;
mean_conn = mean(atanh(conn_mat),2);

end

