function MS_BMA_group(matlabbatch, params, method) % Release cvBMA
% _
% Bayesian Model Averaging over Model Parameters (subject group)
% FORMAT [] = MS_BMA_group(matlabbatch, method, params)
%     matlabbatch - an SPM batch editor job specifying BMS maps
%     params      - an M x P matrix of model parameters to be averaged
%     method      - a string indicating analysis options (see below)
% 
% FORMAT [] = MS_BMA_group(matlabbatch, method, params) performs Bayesian
% model averaging according to the batch editor job matlabbatch and using
% the method indicated by method with parameters of interest indexed by
% params.
% 
% The input variable "method" is a string containing up to five letters:
%     If it contains 'b', then the best model's parameters are extracted.
%     If it contains 'm', then the median model's parameters are extracted.
%     If it contains 'w', then the worst model's parameters are extracted.
%     If it contains 'r', then a random model's parameters are extracted.
%     If it contains 'a', then the averaged model parameters are estimated.
% 
% The input variable "params" can be an M x P matrix indexing parameters of
% interest for all models separately or a 1 x P vector if parameter indices
% are the same for all models, where P is the number of parameters.
% 
% Further information:
%     help MS_BMA_subject
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 03/03/2016, 18:35 (V0.4/V13)
%  Last edit: 11/04/2017, 12:15 (V0.9b/V13b) [retro-spective edit]


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;

% Get matlabbatch if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    design_mat = spm_select(1,'^*.mat','Select Batch Editor Job!');
    load(design_mat);
    MS_BMA_group(matlabbatch);
    return
else
    if isfield(matlabbatch{1}.spm.stats,'bms')  % catch SPM8 batch version
        matlabbatch{1}.spm.stats.bms_map.inference = matlabbatch{1}.spm.stats.bms.bms_map_inf;
        matlabbatch{1}.spm.stats = rmfield(matlabbatch{1}.spm.stats,'bms');
    end;
end;

% Set parameter indices if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(params), params = spm_input('parameter indices:',1,'r','[]'); end;

% Set analysis method if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(method), method = 'ba'; end;


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Get model parameters
%-------------------------------------------------------------------------%
I.sess_map = matlabbatch{1}.spm.stats.bms_map.inference.sess_map;
N = length(I.sess_map);                     % number of subjects
S = length(I.sess_map{1}); % = 1;           % number of sessions
M = length(I.sess_map{1}(1).mod_map);       % number of models

% Get log model evidence maps
%-------------------------------------------------------------------------%
LMEs = cell(N,M);
for i = 1:N
    for j = 1:M
        LMEs{i,j} = I.sess_map{i}(1).mod_map{j};
    end;
end;

% Get model space directories
%-------------------------------------------------------------------------%
BMA_dirs = cell(N,1);
for i = 1:N
    GLM1_dir = fileparts(LMEs{i,1});
    cd(strcat(GLM1_dir,'/../'));
    GLMs_dir = pwd;
    BMA_dirs{i,1} = strcat(GLMs_dir,'/','MS_BMA_subject','_',method);
end;
cd(orig_dir);

% Perform Bayesian model averaging
%-------------------------------------------------------------------------%
fprintf('\n');
for i = 1:N
    fprintf('-> Bayesian model averaging for Subject %d out of %d using Method "%s"... ', i, N, method);
    MS_BMA_subject(LMEs(i,:)', params, method, BMA_dirs{i})
    fprintf('succesful!\n');
end;
fprintf('\n');