
clc;
clear all;
close all;

BRAMILAPATH = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila';
addpath(BRAMILAPATH);
addpath(pwd);

home = pwd;

% for group comparison, use common ROIs and mask!!
roifile = [BRAMILAPATH,filesep,'external',filesep,'rois_Power264.mat'];

iterations = 1e+5;

for run = {'run1','run2'}
    
    maskfile = [pwd,filesep,'testdata',filesep,run{1},filesep,'group_mask.nii'];
    
    k=0;
    k=k+1;niifiles1{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group1',filesep,'S1.nii'];
    k=k+1;niifiles1{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group1',filesep,'S2.nii'];
    k=k+1;niifiles1{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group1',filesep,'S3.nii'];
    k=0;
    k=k+1;niifiles2{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group2',filesep,'S1.nii'];
    k=k+1;niifiles2{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group2',filesep,'S2.nii'];
    k=k+1;niifiles2{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group2',filesep,'S3.nii'];
    k=k+1;niifiles2{k}=[pwd,filesep,'testdata',filesep,run{1},filesep,'group2',filesep,'S4.nii'];
    
    output_path1 =[pwd,filesep,'testdata',filesep,run{1},filesep,'group1'];
    output_path2 =[pwd,filesep,'testdata',filesep,run{1},filesep,'group2'];
    
    %--------------------------------------------------
    T=50;
    interpolation_cfg1.INTERPOLATION_METHOD='linear';
    for i=1:length(niifiles1)
        interpolation_cfg1.volume_time{i}=1:T;
        interpolation_cfg1.requested_time{i}=linspace(1+1e-3,T-1e-3,T*2);
    end
    
    interpolation_cfg2.INTERPOLATION_METHOD='linear';
    for i=1:length(niifiles2)
        interpolation_cfg2.volume_time{i}=1:T;
        interpolation_cfg2.requested_time{i}=linspace(1+1e-3,T-1e-3,T*2);
    end
    
    connISC_run(roifile,maskfile,niifiles1,iterations,output_path1,3,8,interpolation_cfg1);
    connISC_run(roifile,maskfile,niifiles2,iterations,output_path2,3,8,interpolation_cfg2);
    
end