function MS_BMA_subject(LMEs, para, opts, folder) % Release cvBMA [3]
% _
% Bayesian Model Averaging over Model Parameters (single subject)
% FORMAT [] = MS_BMA_subject(LMEs, para, opts, folder)
%     LMEs   - an M x S cell array specifying log model evidence maps
%     para   - an M x P matrix of model parameters to be averaged
%     opts   - a string indicating analysis options (see below)
%     folder - a string indicating the folder for the output images
% 
% FORMAT [] = MS_BMA_subject(LMEs, para, opts, folder) performs Bayesian
% model averaging based on the log model evidence maps LMEs using the
% settings indicated by opts and saves averaged parameters indexed by
% para to sub-directories of the directory folder.
% 
% The model evidence is the probability of the data y (e.g. fMRI signal),
% just given the model m (e.g. a GLM), regardless of any particular
% parameter values b (e.g. regression weights) ([1], eq. 3):
%     p(y|m) = INT p(y|b,m) p(b|m) db
% Then, the LME can be used to calculate posterior probabilties (PP) using
% Bayes' Theorem ([1], eq. 2):
%     p(m_i|y) = [p(y|m_i) * p(m_i)] / SUM_j^M [p(y|m_j) p(m_j)]
% Then, these PPs can be used to calculate averaged model parameters as
% weighted averages ([1], eq. 1):
%     b_BMA = SUM_i^M [b_i * p(m_i|y)]
% These are the three steps of Bayesian model averaging (BMA). In fMRI,
% BMA has been implemented for dynamic causal models (DCMs) [2] and now
% also for general linear models (GLMs) [3]. In our implementation, we
% use the cross-validated log model evidence (cvLME) to calculate PPs.
% 
% The input variable "LMEs" is an M x S cell array where M is the number of
% models, S is the number of sessions (and P is the number of parameters).
% 
% The input variable "para" can be an M x P matrix indexing parameters of
% interest for all models separately or a 1 x P vector if parameter indices
% are the same for all models.
% 
% The input variable "opts" is a string which contains up to six letters:
%     If it contains 'b', then the best model's parameters are extracted.
%     If it contains 'm', then the median model's parameters are extracted.
%     If it contains 'w', then the worst model's parameters are extracted.
%     If it contains 'r', then a random model's parameters are extracted.
%     If it contains 'a', then the averaged model parameters are estimated.
%     If it contains 's', then the session-wise parameters are summed up.
% 
% Further information:
%     help MS_BMA_group
% 
% References:
% [1] Hoeting JA, Madigan D, Raftery AE, Volinsky CT (1999):
%     "Bayesian Model Averaging: A Tutorial".
%     Statistical Science, vol. 14, no. 4, pp. 382-417.
% [2] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [3] Soch J, Meyer AP, Haynes JD, Allefeld C (2017):
%     "How to improve parameter estimates in GLM-based fMRI data analysis:
%      cross-validated Bayesian model averaging". NeuroImage, in review.
%      URL: http://biorxiv.org/content/early/2016/12/20/095778
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 03/03/2016, 13:30 (V0.4/V13)
%  Last edit: 24/02/2017, 15:45 (V0.9b/V13b)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;

% Get model parameters
%-------------------------------------------------------------------------%
M = size(LMEs,1);               % number of models
S = size(LMEs,2);               % number of sessions
P = size(para,2);               % number of parameters

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1,1});         % LME image header
V = prod(H.dim);                % number of voxels

% Expand parameters if necessary
%-------------------------------------------------------------------------%
if size(para,1) == 1
    para = repmat(para,[M 1]);
end;

% Expand options if appropriate
%-------------------------------------------------------------------------%
if S == 1 && ~ismember('s',opts)
    opts = strcat(opts,'s');
end;

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: load');

% Load log model evidences
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');
LME = zeros(S,M,V);             % S x M x V array of LMEs
for i = 1:M                     % models
    for j = 1:S                 % sessions
        lme_hdr = spm_vol(LMEs{i,j});
        lme_img = spm_read_vols(lme_hdr);
        lme_img = reshape(lme_img,[1 1 V]);
        LME(j,i,:) = lme_img;
        spm_progress_bar('Set',(((i-1)*S+j)/(M*S))*100);
    end;
end;
clear lme_hdr lme_img

% Load parameter estimates
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Load parameter estimates...' , '');
B = zeros(S*P,M,V);             % S*P x M x V array of LMEs
for i = 1:M                     % models
    SPM_dir = fileparts(LMEs{i,1});
    SPM_mat = strcat(SPM_dir,'/','SPM.mat');
    load(SPM_mat);
    cd(SPM.swd);
    for j = 1:S                 % sessions
        for k = 1:P             % parameters
            beta_hdr = SPM.Vbeta(SPM.Sess(j).col(para(i,k)));
            beta_img = spm_read_vols(beta_hdr);
            beta_img = reshape(beta_img,[1 1 V]);
            B((j-1)*P+k,i,:) = beta_img;
            spm_progress_bar('Set',(((i-1)*S*P+(j-1)*P+k)/(M*S*P))*100);
        end;
    end;
end;
clear SPM_dir SPM_mat beta_hdr beta_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
LME_1 = squeeze(LME(:,1,:));    % S x V matrix of LMEs for 1st model
if size(LME_1,2) == S, LME_1 = LME_1'; end;
[m_img m_hdr m_ind] = MS_create_mask(LME_1, H);
clear LME_1


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Get best model's parameters
%-------------------------------------------------------------------------%
if ismember('b',opts)
    spm_progress_bar('Init', 100, 'Get best model''s parameters...' , '');
    Bb = NaN(S*P,V);            % best model's betas
    % uncomment for using the cvLME
    % LMEb = squeeze(sum(LME,1));
    % [c,i] = max(LMEb,[],1);
    for j = 1:S                 % for each session
        % uncomment for using the oosLME
        [c,i] = max(squeeze(LME(j,:,:)),[],1);
        for k = 1:P
            for l = 1:v
                Bb((j-1)*P+k,m_ind(l)) = B((j-1)*P+k,i(m_ind(l)),m_ind(l));
            end;
            spm_progress_bar('Set',(((j-1)*P+k)/(S*P))*100);
        end;
    end;
end;

% Get median model's parameters
%-------------------------------------------------------------------------%
if ismember('m',opts)
    spm_progress_bar('Init', 100, 'Get median model''s parameters...' , '');
    Bm = NaN(S*P,V);            % median model's betas
    for j = 1:S                 % for each session
        for k = 1:P
            for l = 1:v
                [c,i] = min(abs(LME(j,:,m_ind(l))-median(LME(j,:,m_ind(l)))));
                Bm((j-1)*P+k,m_ind(l)) = B((j-1)*P+k,i(1),m_ind(l));
            end;
            spm_progress_bar('Set',(((j-1)*P+k)/(S*P))*100);
        end;
    end;
end;

% Get worst model's parameters
%-------------------------------------------------------------------------%
if ismember('w',opts)
    spm_progress_bar('Init', 100, 'Get worst model''s parameters...' , '');
    Bw = NaN(S*P,V);            % worst model's betas
    for j = 1:S                 % for each session
        [c,i] = min(squeeze(LME(j,:,:)),[],1);
        for k = 1:P
            for l = 1:v
                Bw((j-1)*P+k,m_ind(l)) = B((j-1)*P+k,i(m_ind(l)),m_ind(l));
            end;
            spm_progress_bar('Set',(((j-1)*P+k)/(S*P))*100);
        end;
    end;
end;

% Get random model's parameters
%-------------------------------------------------------------------------%
if ismember('r',opts)
    spm_progress_bar('Init', 100, 'Get random model''s parameters...' , '');
    Br = NaN(S*P,V);            % random model's betas
    for j = 1:S                 % for each session
        i = randi([1 M],[1 V]);
        for k = 1:P
            for l = 1:v
                Br((j-1)*P+k,m_ind(l)) = B((j-1)*P+k,i(m_ind(l)),m_ind(l));
            end;
            spm_progress_bar('Set',(((j-1)*P+k)/(S*P))*100);
        end;
    end;
end;

% Get averaged model parameters
%-------------------------------------------------------------------------%
if ismember('a',opts)
    spm_progress_bar('Init', 100, 'Get averaged model parameters...' , '');
    Ba = NaN(S*P,V);            % averaged betas
    for j = 1:S                 % for each session
        % obtain posterior probabilities
        prior = 1/M * ones(M,1);
        LMEj  = squeeze(LME(j,:,:));
        LMEj  = LMEj - repmat(mean(LMEj,1),[M 1]);
        LMEj  = exp(LMEj) .* repmat(prior,[1 V]);
        post  = LMEj ./ repmat(sum(LMEj,1),[M 1]);
        % obtain averaged parameters
        for k = 1:P
            Ba((j-1)*P+k,:) = sum(squeeze(B((j-1)*P+k,:,:)).*post,1);
            spm_progress_bar('Set',(((j-1)*P+k)/(S*P))*100);
        end;
    end;
end;

% Sum up model parameters
%-------------------------------------------------------------------------%
if ismember('s',opts)
    if ismember('b',opts)       % best model's betas
        Bbs = NaN(P,V);         % sum over sessions
        for k = 1:P, Bbs(k,:) = sum(Bb([k:P:S*P],:),1); end;
    end;
    if ismember('m',opts)       % median model's betas
        Bms = NaN(P,V);         % sum over sessions
        for k = 1:P, Bms(k,:) = sum(Bm([k:P:S*P],:),1); end;
    end;
    if ismember('w',opts)       % worst model's betas
        Bws = NaN(P,V);         % sum over sessions
        for k = 1:P, Bws(k,:) = sum(Bw([k:P:S*P],:),1); end;
    end;
    if ismember('r',opts)       % random model's betas
        Brs = NaN(P,V);         % sum over sessions        
        for k = 1:P, Brs(k,:) = sum(Br([k:P:S*P],:),1); end;
    end;    
    if ismember('a',opts)       % averaged betas
        Bas = NaN(P,V);         % sum over sessions
        for k = 1:P, Bas(k,:) = sum(Ba([k:P:S*P],:),1); end;
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1,1});
if ~exist(folder,'dir')
    mkdir(folder);
end;

% Save best model's parameters
%-------------------------------------------------------------------------%
if ismember('b',opts)
    cd(folder);             % create folder
    mkdir('best_model');
    cd('best_model');
    if S ~= 1               % session betas
        for k = 1:P
            for j = 1:S
                H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_S',int2str(j),'.nii');
                H.descrip = 'MS_BMA_subject: session-wise best model''s parameters';
                spm_write_vol(H,reshape(Bb((j-1)*P+k,:),H.dim));
            end;
        end;
    end;
    if ismember('s',opts)   % sum over sessions
        for k = 1:P
            H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'.nii');
            H.descrip = 'MS_BMA_subject: across-session best model''s parameters';
            spm_write_vol(H,reshape(Bbs(k,:),H.dim));
        end;
    end;
end;

% Save median model's parameters
%-------------------------------------------------------------------------%
if ismember('m',opts)
    cd(folder);             % create folder
    mkdir('median_model');
    cd('median_model');
    if S ~= 1               % session betas
        for k = 1:P
            for j = 1:S
                H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_S',int2str(j),'.nii');
                H.descrip = 'MS_BMA_subject: session-wise median model''s parameters';
                spm_write_vol(H,reshape(Bm((j-1)*P+k,:),H.dim));
            end;
        end;
    end;
    if ismember('s',opts)   % sum over sessions
        for k = 1:P
            H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'.nii');
            H.descrip = 'MS_BMA_subject: across-session median model''s parameters';
            spm_write_vol(H,reshape(Bms(k,:),H.dim));
        end;
    end;
end;

% Save worst model's parameters
%-------------------------------------------------------------------------%
if ismember('w',opts)
    cd(folder);             % create folder
    mkdir('worst_model');
    cd('worst_model');
    if S ~= 1               % session betas
        for k = 1:P
            for j = 1:S
                H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_S',int2str(j),'.nii');
                H.descrip = 'MS_BMA_subject: session-wise worst model''s parameters';
                spm_write_vol(H,reshape(Bw((j-1)*P+k,:),H.dim));
            end;
        end;
    end;
    if ismember('s',opts)   % sum over sessions
        for k = 1:P
            H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'.nii');
            H.descrip = 'MS_BMA_subject: across-session worst model''s parameters';
            spm_write_vol(H,reshape(Bws(k,:),H.dim));
        end;
    end;
end;

% Save random model's parameters
%-------------------------------------------------------------------------%
if ismember('r',opts)
    cd(folder);             % create folder
    mkdir('random_model');
    cd('random_model');
    if S ~= 1               % session betas
        for k = 1:P
            for j = 1:S
                H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_S',int2str(j),'.nii');
                H.descrip = 'MS_BMA_subject: session-wise random model''s parameters';
                spm_write_vol(H,reshape(Br((j-1)*P+k,:),H.dim));
            end;
        end;
    end;
    if ismember('s',opts)   % sum over sessions
        for k = 1:P
            H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'.nii');
            H.descrip = 'MS_BMA_subject: across-session random model''s parameters';
            spm_write_vol(H,reshape(Brs(k,:),H.dim));
        end;
    end;
end;

% Save averaged model parameters
%-------------------------------------------------------------------------%
if ismember('a',opts)
    cd(folder);             % change folder
    if S ~= 1               % session betas
        for k = 1:P
            for j = 1:S
                H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_BMA','_S',int2str(j),'.nii');
                H.descrip = 'MS_BMA_subject: session-wise averaged model parameters';
                spm_write_vol(H,reshape(Ba((j-1)*P+k,:),H.dim));
            end;
        end;
    end;
    if ismember('s',opts)   % sum over sessions
        for k = 1:P
            H.fname   = strcat('beta_',MF_int2str0(para(1,k),4),'_BMA','.nii');
            H.descrip = 'MS_BMA_subject: across-session averaged model parameters';
            spm_write_vol(H,reshape(Bas(k,:),H.dim));
        end;
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);