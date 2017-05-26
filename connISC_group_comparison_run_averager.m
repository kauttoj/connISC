function connISC_group_comparison_run_averager(comparison_result_path,datafiles,CPUs,MODE,WEIGHTS)
% datafiles = files saved by 'connISC_group_comparison_run' function

if nargin<5
    WEIGHTS=[];
end

if strcmp(MODE,'perspective')
    MODE=2;
elseif strcmp(MODE,'simple')
    MODE=1;
else
    error('unknown comparison type!')
end

N = length(datafiles);
mean_conn_mat = 0;

if ~isempty(WEIGHTS)
   if length(WEIGHTS)~=N
      error('Number of weights must match to datafiles!') 
   end
   frintf(' Doing WEIGHTED analysis\n ')
else
   WEIGHTS = ones(1,N); 
   frintf(' Doing UN-WEIGHTED analysis\n ')
end

if MODE==2
    
    load(datafiles{1},'perspective_model_mat','perspective_model','isc_inds','perspective_rois','internal_edgelist');
    
    for i=1:N
        load(datafiles{i},'conn_mat');
        if i==1
            mean_conn_mat = WEIGHTS(i)*conn_mat;
        else
            if nnz(size(conn_mat)-size(mean_conn_mat))>0
                error('conn_mat files are not of same size!')
            end
            mean_conn_mat = mean_conn_mat + WEIGHTS(i)*conn_mat;
        end
    end
    mean_conn_mat = mean_conn_mat/sum(WEIGHTS);
    clear conn_mat;
    
    %[~,~,conn_mat,~,~,subject_id] = connISC_compute(roits);
    
    testdata = mean_conn_mat';
    testdata = 1 - testdata; % dissimilarity measure!
    
    % testdata = pairs x voxels
    % perspective_model = pairs x 1
    [isc_connectivity_perspective_corr,isc_connectivity_perspective_pval]=bramila_mantel_vector(testdata,perspective_model,1e+6,'spearman',CPUs,1);
    isc_connectivity_perspective_pval_cor = mafdr(isc_connectivity_perspective_pval,'BHFDR',true);
    
    save([comparison_result_path,filesep,'connISC_comparison_stats_perspective_averaged.mat'],'perspective_model_mat','isc_inds','perspective_model','perspective_rois','isc_connectivity_perspective_corr','isc_connectivity_perspective_pval','isc_connectivity_perspective_pval_cor','mean_conn_mat','datafiles','internal_edgelist','WEIGHTS','-v7.3')
    
elseif MODE==1
    
    load(datafiles{1});
    
    N1=size(group1_res.conn_mat,2);
    N2=size(group2_res.conn_mat,2);
    
    for i=1:N
        load(datafiles{i},'group1_res','group2_res');
        if i==1
            mean_conn_mat = WEIGHTS(i)*[group1_res.conn_mat,group2_res.conn_mat];
        else
            if N1~=size(group1_res.conn_mat,2) || N2~=size(group2_res.conn_mat,2)
                error('group size not consistent!')
            end
            mean_conn_mat = mean_conn_mat + WEIGHTS(i)*[group1_res.conn_mat,group2_res.conn_mat];
        end
    end
    mean_conn_mat = mean_conn_mat/sum(WEIGHTS);
    
    %[~,~,conn_mat,~,~,subject_id] = connISC_compute(roits);
    mean_conn_mat = atanh(mean_conn_mat);
    
    group_model = [ones(1,N1),2*ones(1,N2)];
    
    a = bramila_two_sample_test(mean_conn_mat,group_model,10000,CPUs,0,1);   
    
    isc_connectivity_diff = a.tvals;
    isc_connectivity_pval = 2*min(a.pvals,[],2);
    isc_connectivity_pval_cor = mafdr(isc_connectivity_pval,'BHFDR',true);    
    
    internal_rois=full_rois(selected_rois);
    internal_edgelist=[full2internal(full_edgelist(:,1))',full2internal(full_edgelist(:,2))'];
        
    save([comparison_result_path,filesep,'connISC_comparison_stats_averaged.mat'],'datafiles','isc_connectivity_diff','isc_connectivity_pval','isc_connectivity_pval_cor','mean_conn_mat','selected_rois','internal_rois','internal_edgelist','group_model','WEIGHTS','-v7.3');
    
end

end

