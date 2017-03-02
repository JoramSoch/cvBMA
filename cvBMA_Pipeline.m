% cvBMA Toolkit: Pipeline Template
% _
% Pipeline Template for the cvBMA Toolkit
% For details, see the Manual, page 2-3.
% This script is written for SPM12.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 24/02/2017, 11:00


%%% Step 0: Study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project directories
stat_dir = 'C:\Joram\projects\MACS\DataSets\Meyer_et_al_2015\analyses';
work_dir = 'C:\Joram\projects\MACS\DataSets\Meyer_et_al_2015\subjects';
% Note: This was the data (Meyer & Haynes, in prep.) analyzed in the paper (Soch et al., 2017).

% list of subjects
subj_ids = {'sub01' 'sub02' 'sub03' 'sub04' 'sub05' 'sub06' 'sub07' 'sub08' 'sub09' 'sub10' ...
            'sub11' 'sub12' 'sub13' 'sub14' 'sub15' 'sub16' 'sub17' 'sub18' 'sub19' 'sub20' ...
            'sub21' 'sub22' 'sub23' 'sub24' 'sub25'};

% list of models
mod_nams = {'mod01' 'mod02' 'mod03' 'mod04' 'mod05' 'mod06' 'mod07' 'mod08' 'mod09' 'mod10'};

% list of parameters
para_mat = [1 2 3 4;
            1 2 3 4;
            1 3 4 5;
            1 3 4 5;
            1 2 4 5;
            1 2 4 5;
            1 2 3 5;
            1 2 3 5;
            1 2 3 4;
            1 2 3 4];

% model space details
ms_name  =  'MS01';
ms_suff  =  'BMA1';

% study dimensions
N = numel(subj_ids);
M = numel(mod_nams);
P = size(para_mat,2);


%%% Step 1: First-level model assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate model evidences
for i = 1:N                     % all subjects
    for j = 1:M                 % all models
        load(strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/','SPM.mat'));
        MA_cvLME_multi(SPM);
    end;
end;

% get cvLME filename
load(strcat(work_dir,'/',subj_ids{1},'/',mod_nams{1},'/','SPM.mat'));
S = numel(SPM.Sess);
if S == 1                       % single-session
    data = 'Ky';
    disc = 10 + mod(size(SPM.xX.X,1),10);
    AnC  = false;
    LME_map = strcat('MA_cvLME_',data,'_',int2str(disc));
end;
if S > 1                        % multi-session
    data = 'Ky';
    mode = 'N-1';
    AnC  = false;
    LME_map = strcat('MA_cvLME_',data,'_',mode);
end;
% Note: This only works if the standard settings from "MA_cvLME_single" and
% "MA_cvLME_multi" are used. If these are modified, the cvLME filename will
% also change.


%%% Step 2: First-level model averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare SPM batch
BMA_dir = strcat(stat_dir,'/','model_averaging','/',ms_name,'_',int2str(M),'mods','_',int2str(N),'subj','_',ms_suff,'/');
clear sess_map
for i = 1:N                     % all subjects
    % clear sess_map
    for k = 1:S                 % all sessions
        clear mod_map
        for j = 1:M             % all models
            if S == 1
                mod_map{j,1} = strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/',LME_map,'.nii');
            end;
            if S > 1
                mod_map{j,1} = strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/',LME_map,'_S',int2str(k),'.nii');
            end;
        end;
        sess_map{i}(k).mod_map = mod_map;
    end;
end;

% create SPM batch
matlabbatch{1}.spm.stats.bms_map.inference.dir{1}      = BMA_dir;
matlabbatch{1}.spm.stats.bms_map.inference.sess_map    = sess_map;
matlabbatch{1}.spm.stats.bms_map.inference.mod_name    = mod_nams';
matlabbatch{1}.spm.stats.bms_map.inference.method_maps = 'RFX';
matlabbatch{1}.spm.stats.bms_map.inference.out_file    = 2;
matlabbatch{1}.spm.stats.bms_map.inference.mask        = {''};
matlabbatch{1}.spm.stats.bms_map.inference.nsamp       = '1e6';
filename = strcat(BMA_dir,'design.mat');
save(filename,'matlabbatch');
% Note: Except from sess_map, nothing else matters and in
% principle, all other fields can be left as default or empty.

% load SPM batch
load(strcat(BMA_dir,'design.mat'));
params = para_mat;
method = 'as';

% perform SPM batch
MS_BMA_group(matlabbatch,para_mat,method);
BMA_dir = strcat('MS_BMA_subject_',method);
% Note: MS_BMA_group(matlabbatch,para_mat,'bmwras') invokes additional
% calculation of best/median/worst/random model's parameter estimates.


%%% Step 3: Second-level statistical analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display welcome
fprintf('\n\n');
fprintf('Use these beta image pathes for second-level analyses:\n');

% display beta image pathes
for i = 1:P
    fprintf('\n\n');
    for i = 1:N
        subj_dir = strcat(work_dir,'/',subj_ids{i},'/',BMA_dir);
        beta_img = strcat(subj_dir,'/','beta_',MF_int2str0(para(1,k),4),'_BMA','.nii');
        fprintf('%s\n',beta_img);
    end;
end;

% display farewell
fprintf('\n\n');
fprintf('Good luck and good night!');
fprintf('\n\n');