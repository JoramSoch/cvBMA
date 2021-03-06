function MA_cvLME_single(SPM, data, disc, AnC) % Release cvBMA [1]
% _
% Cross-Validated Log Model Evidence for General Linear Model (single-session)
% FORMAT MA_cvLME_single(SPM, data, disc, AnC)
%     SPM  - a structure specifying an estimated GLM
%     data - a string indicating which data to use (see below)
%     disc - an integer indicating how many volumes to discard (see below)
%     AnC  - a logical indicating accuracy and complexity computation
% 
% FORMAT MA_cvLME_single(SPM, data, disc, AnC) generates a cross-validated
% log model evidence map for a single-session GLM specified by SPM, using
% data indicated by data and discarding a number of volumes indicated by
% disc.
% 
% The present procedure splits the data set into two parts and uses the
% first (second) one to calculate parameter priors for calculating the
% Bayesian log model evidence on the second (first) one. Assumming
% independence between the two parts, the total (cross-validated) log model
% evidence is then equal to the sum of the individual log model evidences.
% 
% The input variable "data" is a string indicating which data to use:
%     If data is 'y',   then the raw data are used.
%     If data is 'Wy',  then the whitened data are used.
%     If data is 'Ky',  then the filtered data are used.
%     If data is 'KWy', then the whitened and filtered data are used.
%     If data is 'WKy', then the filtered and whitened data are used.
% 
% The default for this variable is 'Ky' which means that the analysis
% operates on the filtered data and uses the whitening matrix for
% non-sphericity correction. This is recommended if high-pass filter
% settings are invariant across all models in the model space and
% accomodates for different AR assumptions due to SPM's ReML algorithm.
% 
% The input variable "disc" is an integer indicating how much volumes are
% left out ("discarded") in the middle of the fMRI data set. The default
% for this variable is set such that at least 10 scans are left out and the
% number of remaining scans is divisible by 10, i.e. disc = 10 + mod(n,10).
% 
% The input variable "AnC" is a logical indicating whether model accuracy
% and model complexity are calculated and written to images. The log model
% evidence is the difference of accuracy and complexity: LME = Acc - Com.
% The default for this variable is false.
% 
% Further information:
%     help ME_GLM_NG
%     help ME_GLM_NG_LME
%     help ME_GLM_NG_AnC
% 
% Exemplary usage:
%     MA_cvLME_single(SPM, 'Ky', 15, true);
% 
% References:
% [1] Soch J, Meyer AP, Haynes JD, Allefeld C (2017):
%     "How to improve parameter estimates in GLM-based fMRI data analysis:
%      cross-validated Bayesian model averaging". NeuroImage, in review.
%      URL: http://biorxiv.org/content/early/2016/12/20/095778
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/03/2014, 19:00 (V0.1/V1)
%  Last edit: 24/02/2017, 02:05 (V0.9b/V13b)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    MA_cvLME_single(SPM);
    return
end;

% Estimate model if necessary
%-------------------------------------------------------------------------%
if ~isfield(SPM.xVi,'V')
    SPM_mat = strcat(SPM.swd,'/','SPM.mat');
    MA_GLM_AR_only(SPM_mat); load(SPM_mat);
    MA_cvLME_single(SPM);
    return
end;

% Set data flag if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(data), data = 'Ky'; end;

% Set disc number if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(disc), disc = 10 + mod(size(SPM.xX.X,1),10); end;

% Inactivate AnC if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(AnC), AnC = false; end;

% Call other function if multi-run
%-------------------------------------------------------------------------%
if length(SPM.Sess) > 1
    mode = 'N-1';
    MA_cvLME_multi(SPM,data,mode,AnC);
    return
end;

% Change to SPM.swd if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Get model parameters
%-------------------------------------------------------------------------%
X = SPM.xX.X;                   % design matrix
K = SPM.xX.K;                   % filtering matrix
W = SPM.xX.W;                   % whitening matrix
V = SPM.xVi.V;                  % non-sphericity
n = size(X,1);                  % number of observations
p = size(X,2);                  % number of regressors

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load time series
%-------------------------------------------------------------------------%
% Y = MA_load_data(SPM);
Y = MA_load_data_im(SPM,m_ind);
v = length(m_ind);


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   P A R T I T I O N                       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (1)');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
% Y = Y(:,m_ind);
% v = length(m_ind);

% Preprocess data if required
%-------------------------------------------------------------------------%
if strcmp(data,'y') | strcmp(data,'Wy')
    if size(K.X0,2) > 1         % Delete the lowest frequency in filter to
        K.X0 = K.X0(:,2:end);   % prevent Ln1/2 from being non-invertible
    end;                        % which can happen due to the partition.
    p = p + size(K.X0,2);       % Design matrix is augmented by filter.
end;                            
if strcmp(data,'y')
  % Y = Y;                      % RAW data are used
    X = [X K.X0];               % design matrix must have filter added
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'Wy')
    Y = W*Y;                    % WHITENED data are used
    X = W*[X K.X0];             % design must have filter and be whitened
    P = speye(n);               % precision is equal to identity matrix
end;
if strcmp(data,'Ky')
    Y = spm_filter(K,Y);        % FILTERED data are used
    X = spm_filter(K,X);        % design matrix must be filtered
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'KWy')
    Y = spm_filter(K,W*Y);      % WHITENED and FILTERED data are used
    X = spm_filter(K,W*X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
end;
if strcmp(data,'WKy')
    Y = W*spm_filter(K,Y);      % FILTERED and WHITENED data are used
    X = W*spm_filter(K,X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
end;

% Partition data into two parts (1)
%-------------------------------------------------------------------------%
if mod(n-disc,2) == 0           % EVEN number of scans
    s1 = [1:(n-disc)/2];        % leave out d volumes in the middle
    s2 = [(n-disc)/2+disc+1:n];
end;
if mod(n-disc,2) == 1           % UNEVEN number of scans
    s1 = [1:(n-disc-1)/2];      % leave out d volumes and the last one
    s2 = [(n-disc-1)/2+disc+1:(n-1)];
end;

% Partition data into two parts (2)
%-------------------------------------------------------------------------%
if length(s1) == length(s2)
    Y1 = Y(s1,:);               % time series
    Y2 = Y(s2,:);
    X1 = X(s1,:);               % design matrix
    X2 = X(s2,:);
    P1 = P(s1,s1);              % precision matrix
    P2 = P(s2,s2);
    n1 = numel(s1);             % data points
    n2 = numel(s2);
end;
clear Y X P


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (2)');

% Set (non-informative) priors for both parts
%-------------------------------------------------------------------------%
m0 = zeros(p,1);                % flat Gaussian
L0 = zeros(p,p);
a0 = 0;                         % Jeffrey's prior
b0 = 0;

% Estimate posteriors for 1st part (as priors for 2nd part)
%-------------------------------------------------------------------------%
[mn1, Ln1, an1, bn1] = ME_GLM_NG(Y1, X1, P1, m0, L0, a0, b0, 'Estimate 1st part posteriors (as 2nd part priors)');

% Estimate posteriors for 2nd part (as priors for 1st part)
%-------------------------------------------------------------------------%
[mn2, Ln2, an2, bn2] = ME_GLM_NG(Y2, X2, P2, m0, L0, a0, b0, 'Estimate 2nd part posteriors (as 1st part priors)');

% Set (informative) priors for both parts
%-------------------------------------------------------------------------%
m01 = mn2;  m02 = mn1;          % normal distribution
L01 = Ln2;  L02 = Ln1;
a01 = an2;  a02 = an1;          % gamma distribution
b01 = bn2;  b02 = bn1;
clear mn* Ln* an* bn*

% Estimate posteriors for 1st part (with priors from 2nd part)
%-------------------------------------------------------------------------%
[mn1, Ln1, an1, bn1] = ME_GLM_NG(Y1, X1, P1, m01, L01, a01, b01, 'Estimate 1st part posteriors (with 2nd part priors)');

% Estimate posteriors for 2nd part (with priors from 1st part)
%-------------------------------------------------------------------------%
[mn2, Ln2, an2, bn2] = ME_GLM_NG(Y2, X2, P2, m02, L02, a02, b02, 'Estimate 2nd part posteriors (with 1st part priors)');


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   L O G   M O D E L   E V I D E N C E     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (3)');

% Preallocate images
%-------------------------------------------------------------------------%
log_py1 = NaN(size(M));
log_py2 = NaN(size(M));
if AnC
    Acc_ym1 = NaN(size(M));
    Acc_ym2 = NaN(size(M));
    Com_ym1 = NaN(size(M));
    Com_ym2 = NaN(size(M));
end;

% Calculate evidence for 1st part (using estimates from 2nd part)
%-------------------------------------------------------------------------%
log_py1(m_ind) = ME_GLM_NG_LME(P1, L01, a01, b01, Ln1, an1, bn1);

% Calculate accuracy and complexity for 1st part (using 2nd part)
%-------------------------------------------------------------------------%
if AnC, [Acc_ym1(m_ind), Com_ym1(m_ind)] = ME_GLM_NG_AnC(X1, P1, m01, L01, a01, b01, mn1, Ln1, an1, bn1, 'Compute accuracy and complexity for 1st part'); end;

% Calculate evidence for 2nd part (using estimates from 1st part)
%-------------------------------------------------------------------------%
log_py2(m_ind) = ME_GLM_NG_LME(P2, L02, a02, b02, Ln2, an2, bn2);

% Calculate accuracy and complexity for 2nd part (using 1st part)
%-------------------------------------------------------------------------%
if AnC, [Acc_ym2(m_ind), Com_ym2(m_ind)] = ME_GLM_NG_AnC(X2, P2, m02, L02, a02, b02, mn2, Ln2, an2, bn2, 'Compute accuracy and complexity for 2nd part'); end;

% Calculate total model evidence (assuming independence between parts)
%-------------------------------------------------------------------------%
cvLME = log_py1 + log_py2;
if AnC
    cvAcc = Acc_ym1 + Acc_ym2;
    cvCom = Com_ym1 + Com_ym2;
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);

% Write log model evidence
%-------------------------------------------------------------------------%
H.fname   = strcat('MA_cvLME_',data,'_',int2str(disc),'.nii');
H.descrip = 'MA_cvLME_single: cross-validated log model evidence for general linear model with normal-gamma priors (GLM-NG)';
spm_write_vol(H,reshape(cvLME,m_dim));
H.fname   = strcat('MA_cvLME_',data,'_',int2str(disc),'_P1.nii');
H.descrip = 'MA_cvLME_single: log model evidence for 1st part based on priors estimated from 2nd part';
spm_write_vol(H,reshape(log_py1,m_dim));
H.fname   = strcat('MA_cvLME_',data,'_',int2str(disc),'_P2.nii');
H.descrip = 'MA_cvLME_single: log model evidence for 2nd part based on priors estimated from 1st part';
spm_write_vol(H,reshape(log_py2,m_dim));

% Write accuracy and complexity
%-------------------------------------------------------------------------%
if AnC
    H.fname   = strcat('MA_cvAcc_',data,'_',int2str(disc),'.nii');
    H.descrip = 'MA_cvLME_single: cross-validated model accuracy for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvAcc,m_dim));
    H.fname   = strcat('MA_cvCom_',data,'_',int2str(disc),'.nii');
    H.descrip = 'MA_cvLME_single: cross-validated model complexity for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvCom,m_dim));
    H.fname   = strcat('MA_cvAcc_',data,'_',int2str(disc),'_P1.nii');
    H.descrip = 'MA_cvLME_single: model accuracy for 1st part based on priors estimated from 2nd part';
    spm_write_vol(H,reshape(Acc_ym1,m_dim));
    H.fname   = strcat('MA_cvCom_',data,'_',int2str(disc),'_P1.nii');
    H.descrip = 'MA_cvLME_single: model complexity for 1st part based on priors estimated from 2nd part';
    spm_write_vol(H,reshape(Com_ym1,m_dim));
    H.fname   = strcat('MA_cvAcc_',data,'_',int2str(disc),'_P2.nii');
    H.descrip = 'MA_cvLME_single: model accuracy for 2nd part based on priors estimated from 1st part';
    spm_write_vol(H,reshape(Acc_ym2,m_dim));
    H.fname   = strcat('MA_cvCom_',data,'_',int2str(disc),'_P2.nii');
    H.descrip = 'MA_cvLME_single: model complexity for 2nd part based on priors estimated from 1st part';
    spm_write_vol(H,reshape(Com_ym2,m_dim));
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);