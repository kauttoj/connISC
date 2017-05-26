function [roits,rois,blacklist] = connISC_create_data(niifiles,rois,group_mask,cut_edge_percentage,DROP_BAD,interpolation_cfg)
%CONNISC_CREATE_DATA Summary of this function goes here
%   Detailed explanation goes here    

    if nargin < 5
        DROP_BAD=1;
    end
    if nargin < 6
        interpolation_cfg=[];
    end
    
    parfor i=1:length(niifiles)

        fprintf('...subject %i\n',i)
        nii=load_nii(niifiles{i});
        vol = nii.img;
        temp_cfg=[];
        temp_cfg.infile='';
        temp_cfg.vol=vol;
        temp_cfg.mask = group_mask;
        temp_cfg.rois = rois;
    
        roits{i} = bramila_roiextract(temp_cfg);                
                   
        bad_nodes{i}=find(sum(isnan(roits{i}),1)>0 | var(roits{i})==0);
        
        if length(bad_nodes{i})==0 && ~isempty(interpolation_cfg)
            if length(interpolation_cfg.volume_time{i})~=size(vol,4)
                error('Number of volumes does not match given volume times! (interpolation module)')
            end
            if nnz(isnan(roits{i}))>0
               error('NaN values found! (should not happen here)') 
            end
            if interpolation_cfg.volume_time{i}(1)<=interpolation_cfg.requested_time{i}(1) ...
                    && interpolation_cfg.volume_time{i}(end)>=interpolation_cfg.requested_time{i}(end)
                roits{i} = interp1(interpolation_cfg.volume_time{i},roits{i},interpolation_cfg.requested_time{i},interpolation_cfg.INTERPOLATION_METHOD);
                if nnz(isnan(roits{i}))>0
                   error('Interpolation resulted in NaN values!') ;
                end
            else
               error('Bad interpolation limits!');
            end
        end        
    end        
    
    bad_nodes_total = [];
    for i=1:length(niifiles)
        bad_nodes_total=cat(2,bad_nodes_total,bad_nodes{i});
    end
    bad_nodes_total=unique(bad_nodes_total);
    blacklist = bad_nodes_total;
    if DROP_BAD==1
        rois(blacklist)=[];
    end
    
    T = size(roits{1},1);
    if cut_edge_percentage>1
        edge_vols = floor(cut_edge_percentage/2);
    else
        edge_vols = floor(cut_edge_percentage*T/2);
    end
    
    if edge_vols > 0
        fprintf('Cutting out %i volumes from both edges (%.1f%%)\n\n',edge_vols,(2*edge_vols/T)*100);
    end
    
    for i=1:length(niifiles)
        if DROP_BAD==1
            roits{i}(:,blacklist)=[];
        end
        
        roits{i}(1:edge_vols,:)=[];
        roits{i}((end-edge_vols+1):end,:)=[];
        
        roits{i} = zscore(roits{i});
    end

end

