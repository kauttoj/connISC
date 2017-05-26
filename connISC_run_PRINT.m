function blacklist = connISC_run_PRINT(roifile,maskfile,niifiles,output_path,CPUs,cut_edge_percentage,FILEPRE)

    if nargin<7
        cut_edge_percentage=0;%0.05;
    end
    
    if cut_edge_percentage>0.5
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
fprintf('...%i rois dropped at first stage\n',length(blacklist1));

a = (N_subjects*(N_subjects-1))*(length(rois)*(length(rois)-1)/2); 
b = a*4/(1e+9);
if b>15
    fprintf('Results data size ~%fGB\n',b);
    error('Too large file, reduce ROIs!')
end

fprintf('\n-----Loading data-------\n')
[roits,~,blacklist2] = connISC_create_data(niifiles,rois,mask,cut_edge_percentage,0);
fprintf('...%i rois dropped at second stage\n',length(blacklist2));

%selected_rois(blacklist1)=[];
blacklist = unique(cat(2,blacklist1,blacklist2));
blacklist_ind = blacklist;
%selected_rois(blacklist2)=[];
a=true(1,length(rois));
a(blacklist)=false;
blacklist_bool=a;

blacklist=blacklist_ind;
save([output_path,filesep,'blacklist_ind.mat'],'blacklist','-v7.3');
blacklist=blacklist_bool;
save([output_path,filesep,'blacklist_bool.mat'],'blacklist','-v7.3');

%save([output_path,filesep,'connISC_rois.mat'],'rois','blacklist','roimask','selected_rois','original_rois','-v7.3');

fprintf('\n-----Computing connectivies-------\n')
[~,edgelist,conn_mat,~,conn_mat_internal] = connISC_compute(roits,blacklist_ind);

M = length(rois);
MM = size(edgelist,1);
for i=1:size(conn_mat,2)
    f=sprintf('%s_pair_%i.mat',FILEPRE,i);
    adj = zeros(M,M);
    for j=1:MM
        adj(edgelist(j,1),edgelist(j,2))=conn_mat(j,i);
        adj(edgelist(j,2),edgelist(j,1))=conn_mat(j,i);
    end
    save([output_path,filesep,f],'adj');
    fprintf('...pair %i matrix saved: %s\n',i,f);    
end
for i=1:size(conn_mat_internal,2)
    f=sprintf('%s_subj_%i.mat',FILEPRE,i);
    adj = zeros(M,M);
    for j=1:MM
        adj(edgelist(j,1),edgelist(j,2))=conn_mat_internal(j,i);
        adj(edgelist(j,2),edgelist(j,1))=conn_mat_internal(j,i);
    end
    save([output_path,filesep,f],'adj');
    fprintf('...individual %i matrix saved: %s\n',i,f);    
end



fprintf('\n-----All done!-------\n')