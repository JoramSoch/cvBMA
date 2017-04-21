% cvBMA Toolkit: Hands-On Example
% _
% Hands-On Example for the cvBMA Toolkit
% For details, see the Manual, page 4-7.
% This script is written for SPM8.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 21/04/2017, 15:15


%%% Step 0: Pre-processing and model estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please download the SPM8 Manual or SPM12 Manual and work through Chapter
% 29/31 until the end of Section 29.3/31.3. At this point, you will have
% pre-processed the data and estimated two first-level GLMs, a "categorical
% model" and a "parametric model". [...] After these preliminary analyses,
% the categorical model is located in DIR/categorical and the parametric
% model is located in DIR/parametric where DIR is some folder on your PC.

% analysis directory
DIR = 'C:\Joram\projects\MACS\DataSets\Henson_et_al_2002\stats';


%%% Step 1: First-level model assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% categorical model
load(strcat(DIR,'/classical_categorical/SPM.mat'));
MA_cvLME_single(SPM)

% parametric model
load(strcat(DIR,'/classical_parametric/SPM.mat'));
MA_cvLME_single(SPM)


%%% Step 2: First-level model averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bayesian model selection
params = [1, 4, 7, 10; 1, 2, 5, 6];
load(strcat(DIR,'/BMA_design.mat'));
MS_BMA_group(matlabbatch,params)
% Note: MS_BMA_group(matlabbatch,params,'bmwra') invokes additional
% calculation of best/median/worst/random model's parameter estimates.


%%% Step 3: Ad-hoc parameter estimate check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image calculator
filename = strcat(DIR,'/BMA_check.mat');
spm_jobman('run',filename);

% Prepare visualization
xVis.data     = {strcat(DIR,'/classical_categorical/beta_0001.img');
                 strcat(DIR,'/classical_parametric/beta_0001.img');
                 strcat(DIR,'/MS_BMA_subject_ba/beta_0001_BMA.nii')};
xVis.overlay  = {strcat(DIR,'/MS_BMA_subject_ba/beta_0001_diff.nii')};
xVis.thresh   = '>1';
xVis.PlotType = 'bar';
xVis.LineSpec = 'b';
xVis.XTicks   = {'cat', 'para', 'avg'};
xVis.YLimits  = [];
xVis.Title    = 'categorical vs. parametric vs. average';

% Display visualization
MF_visualize('Setup',xVis);