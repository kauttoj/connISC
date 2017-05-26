function connISC_group_comparison_run(group1_result_path,group2_result_path,comparison_result_path,link_density,CPUs,DO_MANTEL)

if nargin<5
    CPUs = 6;
end
if nargin<6
    DO_MANTEL=0;
end

if exist([comparison_result_path,filesep,'connISC_comparison_data.mat'],'file')
    fprintf('\n-----Old comparison data found, loading-----\n');
    load([comparison_result_path,filesep,'connISC_comparison_data.mat']);
else
    fprintf('\n-----Computing new comparison data-------\n');
    group1_res = load([group1_result_path,filesep,'connISC_results.mat']);
    group1_rois = load([group1_result_path,filesep,'connISC_rois.mat']);
    group2_res = load([group2_result_path,filesep,'connISC_results.mat']);
    group2_rois = load([group2_result_path,filesep,'connISC_rois.mat']);

    selected_rois = intersect(group2_rois.selected_rois,group1_rois.selected_rois);
    selected_rois = sort(selected_rois,'ascend');
    
    for i=1:length(group1_res.roits)
        ind = ismember(group1_rois.selected_rois,selected_rois);
        group1_res.roits{i}=group1_res.roits{i}(:,ind);
    end
    
    for i=1:length(group2_res.roits)
        ind = ismember(group2_rois.selected_rois,selected_rois);
        group2_res.roits{i}=group2_res.roits{i}(:,ind);
    end
    
    group1_sel_ind=[];
    group2_sel_ind=[];
    for i=selected_rois
        ind = find(i==group1_rois.selected_rois);
        if ~isempty(ind)
            group1_sel_ind(end+1) = ind;
        end
        ind = find(i==group2_rois.selected_rois);
        if ~isempty(ind)
            group2_sel_ind(end+1) = ind;
        end
    end
    
    full_rois = group1_rois.original_rois;
    internal_rois = full_rois(selected_rois);
    
    group1_sel=group1_rois.selected_rois(group1_sel_ind);
    group2_sel=group2_rois.selected_rois(group2_sel_ind);
    
    if max(abs(group1_sel-selected_rois))>0 || max(abs(group2_sel-selected_rois))>0
        error('selection list not the same!')
    end
    
    pois1 = ~ismember(group2_rois.selected_rois(group2_res.edgelist(:,1)),selected_rois);
    pois2 = ~ismember(group2_rois.selected_rois(group2_res.edgelist(:,2)),selected_rois);
    pois=pois1 | pois2;
    
    group2_res.edgelist(pois,:)=[];
    group2_res.conn_mat(pois,:)=[];
    group2_res.mean_conn(pois)=[];
    
    group2_res.mean_conn_internal(pois)=[];
    group2_res.conn_mat_internal(pois,:)=[];    
    
    pois1 = ~ismember(group1_rois.selected_rois(group1_res.edgelist(:,1)),selected_rois);
    pois2 = ~ismember(group1_rois.selected_rois(group1_res.edgelist(:,2)),selected_rois);
    pois=pois1 | pois2;
    
    group1_res.edgelist(pois,:)=[];
    group1_res.conn_mat(pois,:)=[];
    group1_res.mean_conn(pois)=[];
    
    group1_res.mean_conn_internal(pois)=[];
    group1_res.conn_mat_internal(pois,:)=[];   
    
    for i=1:size(group1_res.edgelist,1)
        group1_res.edgelist(i,1)=group1_rois.selected_rois(group1_res.edgelist(i,1));
        group1_res.edgelist(i,2)=group1_rois.selected_rois(group1_res.edgelist(i,2));
    end
    
    for i=1:size(group2_res.edgelist,1)
        group2_res.edgelist(i,1)=group2_rois.selected_rois(group2_res.edgelist(i,1));
        group2_res.edgelist(i,2)=group2_rois.selected_rois(group2_res.edgelist(i,2));
    end    
    
    if max(abs(group2_res.edgelist(i,1)-group1_res.edgelist(i,1)))>0 || max(abs(group2_res.edgelist(i,2)-group1_res.edgelist(i,2)))>0
        error('edgelist not the same!')
    end
    
    full_edgelist = group1_res.edgelist;
    
    full2internal = zeros(1,length(full_rois));
    full2internal(selected_rois) = 1:length(selected_rois);          
        
    save([comparison_result_path,filesep,'connISC_comparison_data.mat'],'group2_res','group1_res','selected_rois','full_edgelist','full2internal','full_rois','selected_rois','internal_rois','-v7.3');
    
end

%nulldist = group2_stats.nulldist - group1_stats.nulldist;

N_group1 = size(group1_res.conn_mat,2);
N_group2 = size(group2_res.conn_mat,2);

N_group1_internal = size(group1_res.conn_mat_internal,2);
N_group2_internal = size(group2_res.conn_mat_internal,2);

%N_edges = size(data,1);
N_nodes = length(internal_rois);

%adj_ind = find(triu(ones(N_nodes,N_nodes),1));
internal_edgelist = full_edgelist;
internal_edgelist(:,1)=full2internal(internal_edgelist(:,1));
internal_edgelist(:,2)=full2internal(internal_edgelist(:,2));

data = atanh([group1_res.conn_mat,group2_res.conn_mat]);
fprintf('\n----starting graph computation (isc)----\n');
kk=0;
for dens = link_density
    kk=kk+1;
    
    fprintf('...converting data into adjacency matrices (density %.1f)\n',dens);
    adj_group1 = zeros(N_nodes,N_nodes,N_group1,'int8');
    for i=1:N_group1
        th = prctile(data(:,i),100-dens);       
        edges = internal_edgelist(data(:,i)>th,:);        
        ind = sub2ind([N_nodes,N_nodes],edges(:,1),edges(:,2));
        adj = zeros(N_nodes,N_nodes,'int8');
        adj(ind)=1;
        adj=adj+adj';
        adj_group1(:,:,i)=adj;
    end
    
    adj_group2 = zeros(N_nodes,N_nodes,N_group2,'int8');
    for i=1:N_group2
        th = prctile(data(:,i+N_group1),100-dens);       
        edges = internal_edgelist(data(:,i+N_group1)>th,:);        
        ind = sub2ind([N_nodes,N_nodes],edges(:,1),edges(:,2));
        adj = zeros(N_nodes,N_nodes,'int8');
        adj(ind)=1;
        adj=adj+adj';
        adj_group2(:,:,i)=adj;
    end  
        
    fprintf('...computing node degree (density %.1f)\n',dens);
    degreedata=zeros(N_nodes,N_group1+N_group2);
    for i=1:N_nodes
        nodedegrees1 = squeeze(sum(adj_group1(i,:,:)))';
        nodedegrees2 = squeeze(sum(adj_group2(i,:,:)))';
        degreedata(i,:)=[nodedegrees1,nodedegrees2];
    end
    
    fprintf('...starting two-sample test (density %.1f)\n',dens);
    a = bramila_two_sample_test(degreedata,[ones(1,N_group1),2*ones(1,N_group2)],10000,CPUs,0,1);
    node_stats{kk}.differences = a.tvals;
    node_stats{kk}.pvals = 2*min(a.pvals,[],2);  
    node_stats{kk}.pvals_corrected = mafdr(node_stats{kk}.pvals,'BHFDR',true);

    [~,p]=ttest2(degreedata(:,1:N_group1)',degreedata(:,N_group1 + (1:N_group2))',0.05,'both','unequal');
    node_stats{kk}.pvals_param = p;
    node_stats{kk}.pvals_param_corrected = mafdr(p,'BHFDR',true);

    node_stats{kk}.density=dens;
    node_stats{kk}.degree=degreedata;
end
save([comparison_result_path,filesep,'connISC_comparison_stats.mat'],'node_stats','internal_rois','internal_edgelist','-v7.3')

data_internal = atanh([group1_res.conn_mat_internal,group2_res.conn_mat_internal]);
fprintf('\n----starting graph computation (internal)----\n');
kk=0;
for dens = link_density
    kk=kk+1;
    
    fprintf('...converting data into adjacency matrices (density %.1f)\n',dens);
    adj_group1 = zeros(N_nodes,N_nodes,N_group1_internal,'int8');
    for i=1:N_group1_internal
        th = prctile(data_internal(:,i),100-dens);       
        edges = internal_edgelist(data_internal(:,i)>th,:);        
        ind = sub2ind([N_nodes,N_nodes],edges(:,1),edges(:,2));
        adj = zeros(N_nodes,N_nodes,'int8');
        adj(ind)=1;
        adj=adj+adj';
        adj_group1(:,:,i)=adj;
    end
    
    adj_group2 = zeros(N_nodes,N_nodes,N_group2_internal,'int8');
    for i=1:N_group2_internal
        th = prctile(data_internal(:,i+N_group1_internal),100-dens);       
        edges = internal_edgelist(data_internal(:,i+N_group1_internal)>th,:);        
        ind = sub2ind([N_nodes,N_nodes],edges(:,1),edges(:,2));
        adj = zeros(N_nodes,N_nodes,'int8');
        adj(ind)=1;
        adj=adj+adj';
        adj_group2(:,:,i)=adj;
    end  
        
    fprintf('...computing node degree (density %.1f)\n',dens);
    degreedata=zeros(N_nodes,N_group1_internal+N_group2_internal);
    for i=1:N_nodes
        nodedegrees1 = squeeze(sum(adj_group1(i,:,:)))';
        nodedegrees2 = squeeze(sum(adj_group2(i,:,:)))';
        degreedata(i,:)=[nodedegrees1,nodedegrees2];
    end
    
    fprintf('...starting two-sample test (density %.1f)\n',dens);
    a = bramila_two_sample_test(degreedata,[ones(1,N_group1_internal),2*ones(1,N_group2_internal)],10000,CPUs,0,1);
    node_stats_internal{kk}.differences = a.tvals;
    node_stats_internal{kk}.pvals = 2*min(a.pvals,[],2);  
    node_stats_internal{kk}.pvals_corrected = mafdr(node_stats_internal{kk}.pvals,'BHFDR',true);

    [~,p]=ttest2(degreedata(:,1:N_group1_internal)',degreedata(:,N_group1_internal + (1:N_group2_internal))',0.05,'both','unequal');
    node_stats_internal{kk}.pvals_param = p;
    node_stats_internal{kk}.pvals_param_corrected = mafdr(p,'BHFDR',true);

    node_stats_internal{kk}.density=dens;
    node_stats_internal{kk}.degree=degreedata;
end
save([comparison_result_path,filesep,'connISC_comparison_stats.mat'],'node_stats','internal_rois','node_stats_internal','internal_edgelist','-v7.3')

fprintf('\n----starting two-sample test for links (internal)----\n')
link_stats_internal = bramila_two_sample_test(data_internal,[ones(1,N_group1_internal),2*ones(1,N_group2_internal)],8000,CPUs,0,0);
link_stats_internal.pvals=2*min(link_stats_internal.pvals,[],2);
[~,p]=ttest2(data_internal(:,1:N_group1_internal)',data_internal(:,N_group1_internal + (1:N_group2_internal))',0.05,'both','unequal');
link_stats_internal.pvals_param=p;
if length(link_stats_internal.pvals)>100000
    [link_stats_internal.pvals_corrected,link_stats_internal.qvals] = mafdr(link_stats_internal.pvals);
    [link_stats_internal.pvals_param_corrected,link_stats_internal.qvals_param] = mafdr(link_stats_internal.pvals_param);
else    
    link_stats_internal.pvals_corrected = mafdr(link_stats_internal.pvals,'BHFDR',true);
    link_stats_internal.pvals_param_corrected = mafdr(link_stats_internal.pvals_param,'BHFDR',true);
end

save([comparison_result_path,filesep,'connISC_comparison_stats.mat'],'link_stats_internal','node_stats','internal_rois','node_stats_internal','internal_edgelist','-v7.3')

fprintf('\n----starting two-sample test for links (isc)----\n')
link_stats = bramila_two_sample_test(data,[ones(1,N_group1),2*ones(1,N_group2)],8000,CPUs,0,0);
link_stats.pvals=2*min(link_stats.pvals,[],2);
[~,p]=ttest2(data(:,1:N_group1)',data(:,N_group1 + (1:N_group2))',0.05,'both','unequal');
link_stats.pvals_param=p;
if length(link_stats.pvals)>100000
    [link_stats.pvals_corrected,link_stats.qvals] = mafdr(link_stats.pvals);
    [link_stats.pvals_param_corrected,link_stats.qvals_param] = mafdr(link_stats.pvals_param);
else    
    link_stats.pvals_corrected = mafdr(link_stats.pvals,'BHFDR',true);
    link_stats.pvals_param_corrected = mafdr(link_stats.pvals_param,'BHFDR',true);
end

save([comparison_result_path,filesep,'connISC_comparison_stats.mat'],'link_stats','link_stats_internal','node_stats','internal_rois','node_stats_internal','internal_edgelist','-v7.3')

fprintf('\n----All done (before Mantel)!----\n')

if DO_MANTEL==1
            
    if length(selected_rois)>400
        warning('Too large ROI set (>400), cannot do mantel comparison!')
        return;
    end
             
    if size(group1_res.roits{1},1)~=size(group1_res.roits{1},1)
        warning('groups have different number of volumes, cannot do mantel comparison')
        return;                
    end    
    
    N1=length(group1_res.roits);
    N2=length(group2_res.roits);
    
    %---- TESTINGGGGGGG
%         TTT = size(group1_res.roits{1},1);
%     
%         t11 = detrend(cumsum(randn(TTT,1)));
%         t12 = t11 + 0.2*detrend(cumsum(randn(TTT,1)));        
%         t21 = detrend(cumsum(randn(TTT,1)));
%         t22 = t21 + 0.2*detrend(cumsum(randn(TTT,1)));
%         
%         for i=1:N1
%             for j=1:50
%                 group1_res.roits{i}(:,j)=zscore(t11.*(1+0.05*rand(TTT,1)));
%             end
%             for j=51:127
%                 group1_res.roits{i}(:,j)=zscore(t12.*(1+0.05*rand(TTT,1)));
%             end                        
%         end
%         
%         for i=1:N2
%             for j=1:50
%                 group2_res.roits{i}(:,j)=zscore(t21.*(1+0.05*rand(TTT,1)));
%             end
%             for j=51:127
%                 group2_res.roits{i}(:,j)=zscore(t22.*(1+0.05*rand(TTT,1)));
%             end                        
%         end        
    
    %---- TESTINGGGGGGG
    
    roits=[group1_res.roits,group2_res.roits];
    [~,~,conn_mat,~,~,isc_inds] = connISC_compute(roits);
    
    NN = N1+N2;
    N_pairs = NN*(NN-1);
    if N_pairs~=size(conn_mat,2)
        error('incorrect data size! (bug!)')
    end
    
    testdata = conn_mat';
    testdata = 1 - testdata; % dissimilarity measure!
    
    mat = nan(length(roits),length(roits));
    mat(1:N1,1:N1)=1;
    mat(N1+(1:N2),N1+(1:N2))=2;
    mat((1:N1),N1+(1:N2))=3;
    mat = triu(mat,1);
    mat = mat+mat';
    if nnz(isnan(mat))>0
        error('Perspective matrix has NaNs!')
    end
    
    perspective_model_mat = 0*mat;
    perspective_model_mat(mat==1)=0; % high similarity between group1 viewers
    perspective_model_mat(mat==2)=0; % high similarity between group2 viewers
    perspective_model_mat(mat==3)=1; % low similarity between mixed viewers
    perspective_model = perspective_model_mat(sub2ind(size(mat),isc_inds(:,1),isc_inds(:,2)));
    
    % testdata = pairs x voxels
    % perspective_model = pairs x 1
    [isc_connectivity_perspective_corr,isc_connectivity_perspective_pval]=bramila_mantel_vector(testdata,perspective_model,3e+5,'spearman',CPUs,1);
    isc_connectivity_perspective_pval_cor = mafdr(isc_connectivity_perspective_pval,'BHFDR',true);
    
    perspective_rois = full_rois(selected_rois);
    save([comparison_result_path,filesep,'connISC_comparison_stats_perspective.mat'],'perspective_model_mat','perspective_model','isc_inds','perspective_model','perspective_rois','isc_connectivity_perspective_corr','isc_connectivity_perspective_pval','isc_connectivity_perspective_pval_cor','conn_mat','internal_edgelist','-v7.3')    
    
end

fprintf('\n----All done!----\n')

end

