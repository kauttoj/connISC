clc;
clear all;
close all;

%addpath('/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/july14');
%addpath('/triton/becs/scratch/braindata/kauttoj2/code/connISC_toolbox');

for run = {'run1','run2'}
    
    group1_result_path = [pwd,filesep,'testdata',filesep,run{1},filesep,'group1'];
    group2_result_path = [pwd,filesep,'testdata',filesep,run{1},filesep,'group2'];
    comparison_result_path = [pwd,filesep,'testdata',filesep,run{1}];
    link_density = [5,10]; % in percentages
    
    %------------------------------------
    
    connISC_group_comparison_run(group1_result_path,group2_result_path,comparison_result_path,link_density,3,1);
    
end


