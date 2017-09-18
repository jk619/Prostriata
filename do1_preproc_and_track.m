clc
clear all
project_dir = '/Users/jankurzawski/Dropbox/mrTrix_Prostriata';
freesurfer_init;

hemi = {'lh','rh'};
% rois = {'ProS','_V1'} % for autorois
preproc_anat = 0;
ecc = 0;
dtiInitalize = 0;
register_ROIS = 0;
track_roi = 0;
whole_brain_tractography = 1;
single_rois = 0;
auto_rois = 0;

subjects = {'pro_AB';'pro_CL';'pro_EG';'pro_FM';'pro_GA';'pro_LV';'pro_MC';'pro_IS'};
%%
for s = 5 %: length(subject)
    
    
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life = [config.subject_dir '/life/'];

    config.anat = sprintf('%s%s_t1.nii.gz',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    config.anat_bv = sprintf('%s%s_t1_BV.vmr',config.subject_dir,config.subject);
    config.anat_bv_nii = sprintf('%s%s_t1_BV.nii',config.subject_dir,config.subject);
    config.anat_acpc = sprintf('%st1/t1_acpc.nii',config.subject_dir);
    config.acqparams = sprintf('%s/params/acqparams.txt',project_dir);
    config.index = sprintf('%s/params/index.txt',project_dir);
    
    
    
    mkdir(sprintf('%st1',config.subject_dir));
    mkdir(config.subject_dir_life);

    disp('0. Configuring directories')
    disp(sprintf('-----------------------------------\n\n'))
    
    %% AC-PC allignment before the recon-all to have right information in the header
    
    if preproc_anat
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
        system(['gunzip ' config.anat_acpc '.gz'])
        
        %% Register Rois from BV to acpc
        disp('2. Registering anatomies BV - mrTRix')
        disp(sprintf('-----------------------------------\n\n'))
        anat_bv = xff(config.anat_bv);
        anat_bv.ExportNifti(config.anat_bv_nii);
        mat = spm_coreg(config.anat_acpc,config.anat_bv_nii);
        transformation = spm_matrix(mat);
        system(sprintf('gunzip %s',config.anat_acpc));
        bv2acpc = worldmat2flirtmat(transformation,config.anat_bv_nii,config.anat_acpc);
        save(sprintf('%st1/bv2acpc',config.subject_dir),'bv2acpc','-ascii');
        system(sprintf('flirt -applyxfm -init %st1/bv2acpc -in %s -ref %s -out %st1/%s',config.subject_dir,config.anat_bv_nii,config.anat_acpc,config.subject_dir,'Anatomy_BV_reg.nii'));
    end
    
    if register_ROIS
        
        disp('3. Registering ROIs to AC-PC')
        disp(sprintf('-----------------------------------\n\n'))

        rois_bv = dir([config.subject_dir 'dti/ROIs/bv/*.voi']);
        mkdir([config.subject_dir 'dti/ROIs/anat/']);
        for r = 1 : length(rois_bv)
            anat_bv = xff(config.anat_bv);
            
            anat_bv.VMRData = zeros(256,256,256);
            roi_bv = xff([config.subject_dir 'dti/ROIs/bv/' rois_bv(r).name]);
            IND = sub2ind(size(anat_bv.VMRData),roi_bv.VOI.Voxels(:,1),roi_bv.VOI.Voxels(:,2),roi_bv.VOI.Voxels(:,3));
            anat_bv.VMRData(IND) = 1;
            anat_bv.ExportNifti(sprintf('%s/dti/ROIs/bv/%s.nii',config.subject_dir,rois_bv(r).name(1:end-4)));
            system(sprintf('flirt -applyxfm -init %st1/bv2acpc -in %sdti/ROIs/bv/%s -ref %s -out %sdti/ROIs/anat/%s',config.subject_dir,config.subject_dir,[rois_bv(r).name(1:end-4) '.nii'],config.anat_acpc,config.subject_dir,[rois_bv(r).name(1:end-4) '_anat.nii.gz']));
            system(sprintf('fslmaths %sdti/ROIs/anat/%s -bin %sdti/ROIs/anat/%s',config.subject_dir,[rois_bv(r).name(1:end-4) '_anat.nii.gz'],config.subject_dir,[rois_bv(r).name(1:end-4) '_anat.nii.gz']));
            
        end
        system(sprintf('gunzip %s*',config.rois_dir))
        
    end
    
    %% Run the recon-all -all on the AC-PC aligned T1 image.
    
    
    
    
    %% convert brain.mgz to brain.nii.gz
    %     system(sprintf('mri_convert /Applications/freesurfer/subjects/%s/mri/brain.mgz /Applications/freesurfer/subjects/%s/brain.nii.gz',config.subject,config.subject))
    
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
        movefile(config.bvecs,[config.bvecs(1:end-5) '_uncorr.bvecs'])
        movefile(sprintf('%secc/data_ecc_corr.nii.gz',config.subject_dir),config.dwi)
        movefile(sprintf('%secc/data_ecc_corr.eddy_rotated_bvecs',config.subject_dir),config.bvecs)
        
        
    end
    
    %% dtiInit (creates dt6.mat)
    
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
        dwParams.outDir     = [dwParams.dt6BaseName] ; %path to outdir
        dwParams.dwOutMm    = [2 2 2]; % make sure its isotropic!
        
        % dtiInit(config.dwi, './t1/t1_acpc.nii.gz', dwParams)
        dtiInit(config.dwi, sprintf('%st1/t1_acpc.nii',config.subject_dir), dwParams)
        
    end
    %% Imports the Prostriata ROIS from Glasser template to Subjects space
    % This requires Freesurfer and Glasser template (hemi.HCP-MMP1.annot) inside the fsaverage/label  + Subjects processed with
    % recon-all -all
    
    % breaking the annotation file into sepeate labels with freesurfer
    
    if auto_rois
        
        if ~exist(sprintf('%slabels',config.subject_dir),'dir')
            mkdir(sprintf('%slabels',config.subject_dir))
        end
        
        for h = 1 : length(hemi)
            system(sprintf('%smri_annotation2label --subject fsaverage --hemi %s --annotation HCP-MMP1 --outdir %slabels',freesurfer_string,hemi{h},config.subject_dir))
        end
        
        for r = 1 : length(rois)
            
            roi = rois{r};
            
            
            
            % registering the label to HCP subject
            d = dir(sprintf('%s/labels/*%s*',config.subject_dir,roi));
            
            for p = 1 : length(d)
                
                system(sprintf('%smri_label2label --srclabel %slabels/%s --srcsubject fsaverage --trglabel %s.%s.hcp.label --trgsubject %s --regmethod surface --hemi %s',freesurfer_string,config.subject_dir,d(p).name,d(p).name(1:2),roi,subject,d(p).name(1:2)))
                system(sprintf('%smri_label2vol  --label /Applications/freesurfer/subjects/%s/label/%s.%s.hcp.label --temp /Applications/freesurfer/subjects/%s/brain.nii.gz --o %sdti/ROIs/%s_%s_vol.nii --identity --subject %s',freesurfer_string,config.subject,d(p).name(1:2),roi,config.subject,config.subject_dir,d(p).name(1:2),roi,config.subject))
                
            end
            
        end
        % make ./mri/brain.nii.gz in the subject directory using mri_convert brain.mgz brain.nii.gz and copy it to subject directory
        %converting the surface label into volume
    else
        
        rois = dir([config.rois_dir, '*.nii']);
        
    end
    
    
    %% Make all possible combinations of Rois
    
    
    combinations = nchoosek(1:length(rois),2);
    
    
    %% mrtrix tracking (https://github.com/francopestilli/life_scripts/blob/master/mrtrix_track_between_rois.m)
    
    
    dtFile = fullfile(config.subject_dir, '/dti/dt6.mat');
    
    refImg = fullfile(config.subject_dir, '/t1/t1_acpc.nii');
    fibersFolder = fullfile(config.subject_dir, '/dti/fibers/');
    % Set upt the MRtrix trakign parameters
    trackingAlgorithm = {'prob'};
    lmax    = [4]; % The appropriate value depends on # of directions. For 32, use lower #'s like 4 or 6. For, 6 or 10 is good [10];
    %http://jdtournier.github.io/mrtrix-0.2/tractography/preprocess.html
%     nSeeds  = 60000; % 10000; Total number of fibers that we want to get between 2 rois
%     nFibers = 100000; %1000000; %Maximum number of attempts to get above
    nSeeds  = 1000000; % 10000; Total number of fibers that we want to get between 2 rois
    nFibers = 10000000; %1000000; %Maximum number of attempts to get above
    wmMask  = []; % white matter mask (wmparc.nii)
    wmMaskName= sprintf('%sdti/bin/wmMask',config.subject_dir);
    % [~, wmMaskName] = dtiRoiNiftiFromMat(wmMaskName,refImg,wmMaskName,1); %save the roi when last param equal to 1
    %%
    
    % Then transform the niftis into .mif (mrTrix compatible format)
    [pp,f,~] = fileparts(wmMaskName);
    wmMaskMifName    = fullfile(pp,sprintf('%s.mif',f));
    wmMaskNiftiName  = sprintf('%s.nii.gz',wmMaskName);
    mrtrix_mrconvert(wmMaskNiftiName, wmMaskMifName);
    
    % This first step initializes all the files necessary for mrtrix.
    % This can take a long time.
    files = mrtrix_init(dtFile,lmax,wmMaskName);
    
    if track_roi
        
        for c = 1 : size(combinations,1)
            
            
            %define the ROIs
            % We want to track the cortical pathway (LGN -> V1/V2 and V1/V2 -> MT)
            fromRois = {sprintf('%s%s',config.rois_dir,rois(combinations(c,1)).name)};
            toRois   = {sprintf('%s%s',config.rois_dir,rois(combinations(c,2)).name)};
            
            % Make an (include) white matter mask ROI. This mask is the smallest
            % set of white matter that contains both ROIS (fromRois and toRois)
            %
            % We use a nifti ROi to select the portion of the White matter to use for
            % seeding
            
            % wmMaskName      = fullfile(baseDir,  subjectDir, '/ROIs/front_mask');
            %% Setup the WM mask from Freesurfer
            %     system(sprintf('mri_convert /Applications/freesurfer/subjects/%s/mri/wm.mgz /Applications/freesurfer/subjects/%s/mri/white.nii.gz',subject,subject))
            %     white = niftiRead(sprintf('/Applications/freesurfer/subjects/%s/mri/white.nii.gz',subject))
            %     wm = find(white.data>0 & white.data<255);
            %     newwm = zeros(size(white.data));
            %     newwm(wm) = 1;
            %     white.data = newwm;
            %     niftiWrite(white,'./ROIs/wm_mask.nii.gz')
            
            
            
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
            roipairs{c,1} = sprintf('%s',roi1.name);
            roipairs{c,2} = sprintf('%s',roi2.name);

            
            mkdir(sprintf('%s%s-%s',config.fibers_dir,roi1.name,roi2.name))
            temp_dir = sprintf('%s%s-%s',config.fibers_dir,roi1.name,roi2.name);
            roiUnion        = roi1; % seed union roi with roi1 info
            roiUnion.name   = ['union of ' roi1.name ' and ' roi2.name]; % r lgn calcarine';
            roiUnion.coords = vertcat(roiUnion.coords,roi2.coords);
            roiName         = fullfile(temp_dir,[roi1.name '-' roi2.name '_union']);
            [~, seedMask]   = dtiRoiNiftiFromMat(roiUnion,refImg,[roiName '.nii.gz'],1);
            seedRoiNiftiName= sprintf('%s.nii.gz',roiName);
            seedRoiMifName  = sprintf('%s.mif',roiName);
            
            % Transform the niftis into .mif
            mrtrix_mrconvert(seedRoiNiftiName, seedRoiMifName);
            
            % We cd into the folder where we want to sae the fibers.
            
            %% tractography
            cd(temp_dir);
            
            % We geenrate and save the fibers in the current folder.
            [fibersPDB, status, results] = mrtrix_track_roi2roi(files, [roi_trk{1}], [roi_trk{2}], ...
                seedRoiMifName, wmMaskMifName, trackingAlgorithm{1}, ...
                nSeeds, nFibers);
            myfib = fgRead(fibersPDB);
            fgWrite(myfib,[roi1.name '_' roi2.name '_track'],'mat');
            
            cd(project_dir)
            %% Run Life
            
            % Create two fibers
            %         1. Whole brain without the tract interest
            %         2. Tract of interest
            %         3. To view mrDiffusion
            %         4. To make it pretty MBA toolbox (www.github.com/francopestilli/mba
            %
            
        end
        save([config.subject_dir 'roipairs.mat'],'roipairs')

    end
    %% Full brain tractography
    
    if whole_brain_tractography
        
        mkdir([config.subject_dir '/dti/fibers/wholebrain/']);
        cd([config.subject_dir '/dti/fibers/wholebrain']);
        [fibersPDB_wb, status, results] = mrtrix_track_wholebrain(files, files.wm, ...
            wmMaskMifName, trackingAlgorithm{1},nSeeds, nFibers);
        
        myfib = fgRead(fibersPDB_wb);
        fgWrite(myfib,['wbt_track'],'mat');
        cd(project_dir);
    end
    
    
    %% Single ROIS
    
    rois_mif = dir([config.rois_dir '*.mif']);
    
    if single_rois
        
        for mif = 2 : length(rois_mif)
            
            floor1 = strfind(rois_mif(mif).name,'_');
            roi_name = rois_mif(mif).name(1:floor1(2)-1);
            
            
            mkdir([config.fibers_dir roi_name]);
            cd([config.fibers_dir roi_name]);
            [fibersPDB_wb, status, results] = mrtrix_track_wholebrain(files, [config.rois_dir rois_mif(mif).name]  , ...
                wmMaskMifName, trackingAlgorithm{1},nSeeds, nFibers);
            
            myfib = fgRead(fibersPDB_wb);
            fgWrite(myfib,[roi_name],'mat');
            cd(project_dir);
        end
    end
    
end