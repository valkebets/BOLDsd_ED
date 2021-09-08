% This script computes a Partial Least Squares (PLS) correlation analysis
% between BOLD signal variability and clinical variables
%


clear

addpath(genpath('./myPLS'));


%% Specify variables/parameters for PLS analysis (in myPLS toolbox)

% Open myPLS_inputs (in myPLS toolbox) and change the following variables:


% Comment lines 84-86: load('example_data.mat'); X0 = brain_data; Y0behav = [age,sex,FSIQ]; % sex - 0=female/1=male
input.brain_data = % enter matrix containing BOLD signal variability data (Subjects x Voxels)
input.behav_data = % enter matrix containing clinical data (Subjects x Clinical variables)

%%%% Regress out potential confounds before running PLS analysis !!

input.group_names = {'ADHD','Bipolar','Borderline','Control'}; 
input.behav_names = {'ALS','MADRS','YMRS'}; 
% Comment line 120: for ii = 1:20; input.img_names{ii,1} = ['img ' num2str(ii)]; end
pls_opts.nPerms = 1000;
pls_opts.nBootstraps = 1000;
pls_opts.normalization_img = 1;
pls_opts.normalization_behav = 1;
pls_opts.grouped_PLS = 0; 
pls_opts.grouped_perm = 1;
pls_opts.grouped_boot = 1;
pls_opts.boot_procrustes_mod = 1;
pls_opts.save_boot_resampling = 0;
save_opts.output_path = './results';
save_opts.alpha = 0.05;
% Uncomment line 203: % save_opts.img_type = 'volume'; 
% Comment lines 210-212: input.brain_data = input.brain_data(:,1:20); save_opts.img_type = 'barPlot'; %save_opts.fig_pos_img = [440   606   560   192];
save_opts.BSR_thres = [-3 3];
save_opts.mask_file = % enter filename of binary group template that will constrain analysis
save_opts.errorbar_mode = 'std'; 

%% Run PLS analysis (myPLS toolbox)

myPLS_main
