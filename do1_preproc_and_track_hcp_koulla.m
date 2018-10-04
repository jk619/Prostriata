clc
clear all

% This script is used to initialize dt6 with common pipeline (AC-PC
% registration, dtiInit etc...), alliging ROIS 2 acpc and tracking using streamtrack
%%
project = 'Koulla'


if isempty(strfind(cd,'Volumes')) == 0
    project_dir = sprintf('/Volumes/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'Users')) == 0
    
    
    project_dir = sprintf('/Users/jankurzawski/Dropbox/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'media')) == 0
    
    addpath(genpath('/home/fmridti/Documents/software'))
    [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/media/fmridti/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    setenv('PATH','/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/freesurfer/bin/')
    addpath(genpath('/home/fmridti/Documents/software'))
    
end



cmap = [round([0 189 61]/255,2);round([237 26 55]/255,2);round([151 189 61]/255,2);round([128 0 128]/255,2);round([0 0 255]/255,2)]

hemi = {'lh','rh'};

hcp = 1; %1 with lines 49, 51 normalize bvals
preproc_anat = 0; %transformation of t1 anatomy to acpc; creates a folder 't1'
do_2acpc = 0; %creates the acpc matrix from line above
wmmask = 0; %takes the wm mask from freesurfer dir and transforms it into acpc space
ecc = 0;
do_regrois = 0; % transforming rois into acpc space
dtiInitalize = 0 ; % initalize dt6
track_roi = 1; %1





subjects = {'101107';'101309';'101410';'101915';'102008';'100206';'100307';'100408';'100610';'101006'}; %done 101006  100610


lmax_list = [8];
tic
%%
for s = 1:  length(subjects);% 3 7 8 9 ran just before Jan left March2018
    
    
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life  = [config.subject_dir 'life/'];
    
    config.anat = sprintf('%s%s_t1.nii.gz',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    config.anat_acpc = sprintf('%st1/t1_acpc.nii',config.subject_dir);
    config.anat_acpc_ref = sprintf('%st1/t1.nii',config.subject_dir);
    config.acqparams = sprintf('%s/params/acqparams.txt',project_dir);
    config.index = sprintf('%s/params/index.txt',project_dir);
    [~,config.freesurfer_dir] = system('echo $SUBJECTS_DIR');
    
    config.rois_acpc = [config.rois_dir 'acpc/'];
    
    
    mkdir(sprintf('%st1',config.subject_dir));
    mkdir(config.subject_dir_life);
    
    disp('0. Configuring directories')
    disp(sprintf('-----------------------------------\n\n'))
    
    params.thresholds.b0_normalize    = 200;
    params.thresholds.bvals_normalize = 100;
    
    %% Normalize HCP files to the VISTASOFT environment
    if hcp
        bvals.val = dlmread(config.bvals);
        
        % Round the numbers to the closest thousand
        % This is necessary because the VISTASOFT software does not handle the B0
        % when they are not round numbers
        [bvals.unique, ~, bvals.uindex] = unique(bvals.val);
        
        bvals.unique(bvals.unique <= params.thresholds.b0_normalize) = 0;
        bvals.unique  = round(bvals.unique./params.thresholds.bvals_normalize) ...
            *params.thresholds.bvals_normalize;
        bvals.valnorm = bvals.unique( bvals.uindex );
        dlmwrite(config.bvals,bvals.valnorm,'delimiter',' ');
    end
    
    %% AC-PC allignment before the recon-all to have right information in the header
    
    if preproc_anat
        
        
        % convert orig from Fs to nifti
        system(sprintf('mri_convert %s/%s/mri/orig.mgz %s',config.freesurfer_dir(1:end-1),subjects{s},config.anat_acpc_ref))
        
        % Load the file from disk
        disp('1. Registering to AC-PC')
        disp(sprintf('-----------------------------------\n\n'))
        ni = niftiRead(config.anat);
        
        % Make sure the file is aligned properly
        ni = niftiApplyCannonicalXform(ni);
        
        % Load a standard template from vistasoft
        MNI_template =  fullfile(mrDiffusionDir, 'templates', 'MNI_T1.nii.gz');
        
        % Compute the spatial normalization to align the current raw data to the template
        SpatialNormalization = mrAnatComputeSpmSpatialNorm(ni.data, ni.qto_xyz, MNI_template);
        
        % Assume that the AC-PC coordinates in the template are in a specific location:
        % X, Y, Z = [0,0,0; 0,-16,0; 0,-8,40]
        % Use this assumption and the spatial normalization to extract the corresponding AC-PC location on the raw data
        coords = [0,0,0; 0,-16,0; 0,-8,40]; % Coordinats of AC an PC
        
        ImageCoords = mrAnatGetImageCoordsFromSn(SpatialNormalization, tal2mni(coords)', true)';
        
        % Now we assume that ImageCoords contains the AC-PC coordinates that we need for the Raw data.
        % We will use them to compute the AC_PC alignement automatically. The new file will be saved to disk.
        % Check the alignement.
        mrAnatAverageAcpcNifti(ni, config.anat_acpc, ImageCoords, [], [], [], false); % saves a file AC-PC alligned
        system(['gunzip -f ' config.anat_acpc '.gz'])
    end
    
    if do_2acpc %allign ROIS
        
        mat = spm_coreg(config.anat_acpc,config.anat_acpc_ref);
        transformation = spm_matrix(mat);
        roi2acpc = worldmat2flirtmat(transformation,config.anat_acpc_ref,config.anat_acpc);
        save(sprintf('%st1/roi2acpc',config.subject_dir),'roi2acpc','-ascii');
        system(sprintf('flirt -applyxfm -init %st1/roi2acpc -in %s -ref %s -out %st1/%s',config.subject_dir,config.anat_acpc_ref,config.anat_acpc,config.subject_dir,'t12t_acpc.nii'));
        
    end
    
    
    
    
    %% ECC fsl
    
    if ecc
        
        disp('4. Running eddy current correction on DTI data')
        disp(sprintf('-----------------------------------\n\n'))
        mkdir(sprintf('%s/ecc',config.subject_dir))
        copyfile(sprintf('%s',config.dwi),sprintf('%s/ecc/ecc.nii.gz',config.subject_dir))
        system(sprintf('fslroi %secc/ecc.nii.gz %secc/%s 0 1',config.subject_dir,config.subject_dir,'b0'));
        system(sprintf('bet %secc/%s %secc/%s -f 0.05 -m',config.subject_dir,'b0',config.subject_dir,'brain'));
        system(sprintf('eddy --imain=%secc/ecc --out=%secc/%s --mask=%secc/brain_mask.nii.gz --acqp=%s --index=%s --bvecs=%s --bvals=%s'...
            ,config.subject_dir,config.subject_dir,'data_ecc_corr',config.subject_dir,config.acqparams,config.index,config.bvecs,config.bvals));
        
        movefile(config.dwi,[config.dwi(1:end-7) '_uncorr.nii.gz'])
        movefile(config.bvecs,[config.bvecs(1:end-6) '_uncorr.bvecs'])
        movefile(sprintf('%secc/data_ecc_corr.nii.gz',config.subject_dir),config.dwi)
        movefile(sprintf('%secc/data_ecc_corr.eddy_rotated_bvecs',config.subject_dir),config.bvecs)
        
        
    end
    
    %% dtiInit (creates dt6.mat and all necessary files; also makes it perfectly isotropic)
    
    if dtiInitalize
        
        disp('5. Creating dt6.mat file')
        disp(sprintf('-----------------------------------\n\n'))
        
        dwi = niftiRead(config.dwi);
        res = dwi.pixdim(1:3);
        dwParams = dtiInitParams;
        dwParams.clobber           =  1; %if 1 and processed already delete and rerun ; 0 not ovewrite
        dwParams.eddyCorrect       = -1; % if -1 only align dwi to T1 otherwise 1 = all correction 0 = motion correction
        dwParams.phaseEncodeDir    = 2;
        dwParams.rotateBvecsWithRx = 1; % allign bvecs to diffusion space
        dwParams.rotateBvecsWithCanXform = 1; % allign to bvecs to T1 space
        dwParams.bvecsFile  = config.bvecs; %paths to bvecs
        dwParams.bvalsFile  = config.bvals; %paths to bvals
        dwParams.dt6BaseName = sprintf('%sdti/',config.subject_dir); %folder name
        dwParams.dt6BaseName = sprintf('%sdti/',config.subject_dir); %folder name
        dwParams.outDir     = [dwParams.dt6BaseName] ; %path to outdir
        dwParams.dwOutMm    = [1.5 1.5 1.5]; % make sure its isotropic!
        
        % dtiInit(config.dwi, './t1/t1_acpc.nii.gz', dwParams)
        dtiInit(config.dwi, sprintf('%st1/t1_acpc.nii',config.subject_dir), dwParams)
        
    end
    
    
    if wmmask % this creates a wm mask from freesurfer segmentation and uses it to constrain tracking
        
        mkdir(sprintf('%s/dti/bin/',config.subject_dir))
        system(sprintf('mri_convert %s/%s/mri/wm.mgz %s/%s/mri/white.nii.gz',config.freesurfer_dir(1:end-1),subjects{s},config.freesurfer_dir(1:end-1),subjects{s}));
        white = niftiRead(sprintf('/usr/local/freesurfer/subjects/%s/mri/white.nii.gz',subjects{s}))
        wm = find(white.data>0 & white.data<255);
        newwm = zeros(size(white.data));
        newwm(wm) = 1;
        white.data = newwm;
        wm_freesurfer_nii = sprintf('%sdti/bin/wm_freesurfer.nii.gz', config.subject_dir)
        niftiWrite(white,wm_freesurfer_nii);
        system(sprintf('flirt -applyxfm -init %st1/roi2acpc -in %s -ref %s -out %s',config.subject_dir,wm_freesurfer_nii,config.anat_acpc,wm_freesurfer_nii))
        wm_bin = load_nifti(wm_freesurfer_nii)
        wm_bin.vol = wm_bin.vol > 0;
        save_nifti(wm_bin,wm_freesurfer_nii)
        wm_freesurfer_mif = sprintf('%sdti/bin/wm_freesurfer.mif', config.subject_dir)
        %         system(sprintf('mrconvert %s %s',wm_freesurfer_nii,wm_freesurfer_mif))
        
    end
    
    wm_freesurfer_mif = sprintf('%sdti/bin/wm_freesurfer.mif', config.subject_dir)
    
    
    %% Loop through ROIS and apply to each of them the acpc matrix and save as nifti
    for l = 1: length(lmax_list)
        for h = 1 :length(hemi)
            
            hem = hemi{h};
            
            rois = dir([config.rois_dir, sprintf('*%s*.nii',hem)]);
            
            
            if do_regrois
                for r = 1 : length(rois)
                    
                    system(sprintf('flirt -applyxfm -init %st1/roi2acpc -in %s%s -ref %s -out %s%s',config.subject_dir,config.rois_dir,rois(r).name,config.anat_acpc,config.rois_acpc,[rois(r).name(1:end-4) '_acpc.nii']));
                    system(sprintf('gunzip -f %s%s',config.rois_acpc,[rois(r).name(1:end-4) '_acpc.nii.gz']))
                    myroi = load_nifti(sprintf('%s%s',config.rois_acpc,[rois(r).name(1:end-4) '_acpc.nii']))
                    myroi.vol = myroi.vol > 0;
                    save_nifti(myroi,(sprintf('%s%s',config.rois_acpc,[rois(r).name(1:end-4) '_acpc.nii'])))
                end
            end
            rois = dir([config.rois_acpc, sprintf('*%s*.nii',hem)]);
            %% Make all possible combinations of Rois
            
            
            combinations = nchoosek(1:length(rois),2);
            
            
            %% mrtrix tracking (https://github.com/francopestilli/life_scripts/blob/master/mrtrix_track_between_rois.m)
            
            if track_roi
                
                
                
                dtFile = fullfile(config.subject_dir, '/dti/dt6.mat');
                
                
                
                %         refImg = fullfile(config.subject_dir, '/t1/t1_acpc.nii');
                fibersFolder = fullfile(config.subject_dir, '/dti/fibers/');
                % Set upt the MRtrix trakign parameters
                trackingAlgorithm = {'prob'};
                lmax    = lmax_list(l);
                % The appropriate value depends on # of directions. For 32, use lower #'s like 4 or 6. For, 6 or 10 is good [10];
                %http://jdtournier.github.io/mrtrix-0.2/tractography/preprocess.html
                nSeeds  = 10000; % 10000; Total number of fibers that we want to get between 2 rois
                nFibers = 1000000; %1000000; %Maximum number of attempts to get above
                wmMaskName= sprintf('%sdti/bin/wmMask',config.subject_dir);
                
                % Then transform the niftis into .mif (mrTrix compatible format)
                [pp,f,~] = fileparts(wmMaskName);
                wmMaskMifName    = fullfile(pp,sprintf('%s.mif',f));
                wmMaskNiftiName  = sprintf('%s.nii.gz',wmMaskName);
                mrtrix_mrconvert(wmMaskNiftiName, wmMaskMifName);
                
                % This first step initializes all the files necessary for mrtrix.
                % This can take a long time.
                
                files = mrtrix_init(dtFile,lmax,wmMaskName);
                
                for c = 1: size(combinations,1)
                    
                    
                    %% Define the ROIs
                    % We want to track the cortical pathway (LGN -> V1/V2 and V1/V2 -> MT)
                    fromRois = {sprintf('%s%s',config.rois_acpc,rois(combinations(c,1)).name)};
                    toRois   = {sprintf('%s%s',config.rois_acpc,rois(combinations(c,2)).name)};
                    
                    % Some of the following steps only need to be done once for each ROI,
                    % so we want to do some sort of unique operation on the from/toRois
                    
                    individualRois = unique([fromRois, toRois]); % put together all the pairs of ROIs
                    
                    % Convert the ROIs from .mat or .nii.gz to .mif format.
                    for i_roi = 1:length(individualRois)
                        
                        
                        [f,p,e] = fileparts(individualRois{i_roi});
                        
                        roi_trk{i_roi} = [f '/' p '.mif'];
                        
                        mrtrix_mrconvert(individualRois{i_roi},roi_trk{i_roi});
                        
                    end
                    
                    %%
                    % Create joint from/to Rois to use as a mask
                    
                    % MRTRIX tracking between 2 ROIs template.
                    
                    roi1 = dtiRoiFromNifti([individualRois{1}],[],[],'.mat');
                    roi2 = dtiRoiFromNifti([individualRois{2}],[],[],'.mat');
                    
                    % Make a union ROI to use as a seed mask:
                    % We will generate as many seeds as requested but only inside the voume
                    % defined by the Union ROI.
                    %
                    % The union ROI is used as seed, fibers will be generated starting ONLy
                    % within this union ROI.
                    
                    floor1 = strfind(roi1.name,'_');
                    floor2 = strfind(roi2.name,'_');
                    
                    roi1.name = roi1.name(1:floor1(2)-1);
                    roi2.name = roi2.name(1:floor2(2)-1);
                    
                    
                    
                    roi1data = load_nifti(individualRois{1});
                    roi2data = load_nifti(individualRois{2});
                    
                    mkdir(sprintf('%s%s-%s',config.fibers_dir,roi1.name,roi2.name))
                    temp_dir = sprintf('%s%s-%s',config.fibers_dir,roi1.name,roi2.name);
                    roiUnion        = roi1data; % seed union roi with roi1 info
                    roiUnion.name   = ['union of ' roi1.name ' and ' roi2.name]; % r lgn calcarine';
                    roiUnion.vol = roiUnion.vol + roi2data.vol;
                    roiName         = fullfile(temp_dir,[roi1.name '-' roi2.name '_union']);
                    save_nifti(roiUnion,[roiName '.nii.gz']);
                    
                    seedRoiNiftiName= sprintf('%s.nii.gz',roiName);
                    seedRoiMifName  = sprintf('%s.mif',roiName);
                    
                    % Transform the niftis into .mif
                    mrtrix_mrconvert(seedRoiNiftiName, seedRoiMifName);
                    
                    % We cd into the folder where we want to sae the fibers.
                    
                    %% tractography
                    cd(temp_dir);
                    curv = [0.25 0.5 1 2 4];
                    for x = 1:length(curv)
                        %We geenrate and save the fibers in the current
                        %folder. Ensemble roi2roi just runs the same
                        %streamtrack command with different curvatures and
                        %saves them as fib files
                        [fibersPDB, status, results] = mrtrix_track_roi2roi_ensamble(files, [roi_trk{1}], [roi_trk{2}], ...
                            seedRoiMifName, wm_freesurfer_mif, trackingAlgorithm{1}, ...
                            nSeeds, nFibers,[],[],num2str(curv(x)));
                        
                        myfib = fgRead(fibersPDB);
                        fgWrite(myfib,sprintf('%s-%s_l%i_c%s',roi1.name,roi2.name,lmax_list(l),num2str(curv(x))),'mat');
                    end
                    cd(project_dir)
                    
                end
            end
        end
        cd(project_dir)
    end
end


toc
