function connISC_run(roifile,maskfile,niifiles,iterations,output_path,CPUs,cut_edge_percentage,interpolation_cfg)

    if nargin<7
        cut_edge_percentage=0;%0.05;
    end
    
    if nargin<8
        interpolation_cfg=[];
    end
    
    if ~isempty(interpolation_cfg)
       if length(interpolation_cfg.volume_time)~=length(niifiles) || length(interpolation_cfg.requested_time)~=length(niifiles)
           error('Interpolation data does not match niftidata!')
       end                
    end
    
    if cut_edge_percentage<1 && cut_edge_percentage>0.5
       warning('!!Very high cut_edge_percentage!!!');
    end
    
% create local cluster
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        if nargin>5
            myCluster.NumWorkers=CPUs;
        end
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
        if nargin>5
            eval(['matlabpool local ',num2str(CPUs)]);
        else
            matlabpool local
        end
    end
    N_workers = matlabpool('size');
end

fprintf('\n-----Loading ROI data-------\n')

A = load(roifile);
rois = A.rois;
original_rois = rois;
N_nodes = length(rois);
selected_rois = 1:length(rois);
clear A;

nii=load_nii(maskfile);
mask=nii.img;
clear nii;

N_subjects = length(niifiles);

[blacklist1,roimask] = bramila_maskrois(rois,mask,0.25);
%[rois,blacklist1,roimask] = bramila_maskrois(rois,mask,0.25);
rois(blacklist1)=[];
fprintf('...%i rois dropped at first stage\n',length(blacklist1));

a = (N_subjects*(N_subjects-1))*(length(rois)*(length(rois)-1)/2); 
b = a*4/(1e+9);
if b>15
    fprintf('Results data size ~%fGB\n',b);
    error('Too large file, reduce ROIs!')
end

fprintf('\n-----Loading data-------\n')
[roits,rois,blacklist2] = connISC_create_data(niifiles,rois,mask,cut_edge_percentage,1,interpolation_cfg);
fprintf('...%i rois dropped at second stage\n',length(blacklist2));

selected_rois(blacklist1)=[];
blacklist = cat(2,blacklist1,selected_rois(blacklist2));
selected_rois(blacklist2)=[];

if N_nodes ~= length(selected_rois)+length(blacklist)
	error('roi-set sizes do not match!')
end

save([output_path,filesep,'connISC_rois.mat'],'rois','blacklist','roimask','selected_rois','original_rois','-v7.3');

fprintf('\n-----Computing connectivies-------\n')
[mean_conn,edgelist,conn_mat,mean_conn_internal,conn_mat_internal,isc_inds] = connISC_compute(roits);

save([output_path,filesep,'connISC_results.mat'],'mean_conn','edgelist','conn_mat','mean_conn_internal','conn_mat_internal','roits','isc_inds','-v7.3');

fprintf('\n-----Computing simple statistics-------\n')
nulldist = connISC_stats(roits,iterations);

save([output_path,filesep,'connISC_nulldist.mat'],'nulldist','-v7.3');

fprintf('\n-----All done!-------\n')