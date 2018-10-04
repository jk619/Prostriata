
clear all
close all

project = 'Koulla'


%setup paths

if isempty(strfind(cd,'Volumes')) == 0
    
    project_dir = sprintf('/Volumes/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'Users')) == 0
    
    
    project_dir = sprintf('/Users/jankurzawski/Dropbox/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'media')) == 0
    
    setenv('PATH''/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/home/fmridti/Documents/mrtrix3')
    addpath(genpath('/home/fmridti/Documents/software'))
    [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/media/fmridti/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    setenv('PATH','/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
    
    addpath(genpath('/home/fmridti/Documents/software'))
    
end


[status,cmdout] = system('echo $HOME')



%%
freesurfer_init;


subjects = {'101107';'101309';'101410';'101915';'102008';'100206';'100307';'100408';'100610';'101006'}; %done 101006  100610

rois_fs = {'RSC';'p24_';} % names of ROIS from Glasser


rois_fsl = {'lgn_lh';'lgn_rh';'pulv_rh';'pulv_lh'} % strings for MNI template
rois_fsl_num = [497 500 623 622] % numbers that Coresspond to above structures in atlas MNI from FSL


hemi = {'lh','rh'};
% hemi = hemi(2)
cortical_rois = 1; %  if 1 mapping cortical ROIs from Glasser template
subcortical_rois = 0; % if 1 mapping subcortical ROIs from Talaraich
reg2mni = 1;
%%
for s = 1:length(subjects)
    %setup config
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life = [config.subject_dir 'life'];
    config.anat = sprintf('%s%s_t1.nii',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.anat_dir = sprintf('%st1/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    config.anat_acpc = sprintf('%st1/t1_acpc.nii',config.subject_dir);
    config.anat_acpc_freesurf = sprintf('%st1/t1.nii',config.subject_dir);
    dwiFile = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    dwiFileRepeat = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    [~,config.freesurfer_dir] = system('echo $SUBJECTS_DIR');
    config.fsl_atlas = '/usr/local/fsl/data/'
    config.rois_acpc = [config.rois_dir 'acpc/'];

    t1File = config.anat_acpc;
    Niter = 500;
    
    mkdir(config.rois_dir)
    system(sprintf('%smri_convert %s/%s/mri/brain.mgz %s/%s/brain.nii.gz ',freesurfer_string,config.freesurfer_dir(1:end-1),subjects{s},config.freesurfer_dir(1:end-1),subjects{s}))
    
    %%
    if cortical_rois
        
        if exist(sprintf('%slabels',config.subject_dir),'dir') 
           
            
        else
            mkdir(sprintf('%slabels_PALS',config.subject_dir))
            % first break Glasser template into ROIS for each hemi
            
            for h = 1 : length(hemi)
                system(sprintf('%smri_annotation2label --subject fsaverage --hemi %s --annotation HCPMMP1 --outdir %slabels_PALS',freesurfer_string,hemi{h},config.subject_dir))
            end
            
        end
        
        for r = 1 : length(rois_fs)
            
            roi = rois_fs{r};
            
            for h = 1 : length(hemi)
                
                % find labels of interest from the folder
                d = dir(sprintf('%slabels/%s.*%s*',config.subject_dir,(hemi{h}),roi));
                
                % 1) map the ROI from native space to freesurfer space,
                % 2) From subject native surface 2 volume
                for p = 1: length(d)
                    
                    system(sprintf('%smri_label2label --srclabel %slabels/%s --srcsubject fsaverage --trglabel %s.%s.hcp.label --trgsubject %s --regmethod surface --hemi %s',freesurfer_string,config.subject_dir,d(p).name,d(p).name(1:2),roi,config.subject,d(p).name(1:2)))
                    system(sprintf('%smri_label2vol  --label %s/%s/label/%s.%s.hcp.label --temp %s/%s/brain.nii.gz --o %sdti/ROIs/anat/%s_%s_anat.nii --identity --subject %s',freesurfer_string,config.freesurfer_dir(1:end-1),config.subject,d(p).name(1:2),roi,config.freesurfer_dir(1:end-1),config.subject,config.subject_dir,lower(roi),d(p).name(1:2),config.subject))
                end        
            end
        end
        
        
        %%
        if subcortical_rois
            
            if reg2mni
                % 1) register Anatomy to MNI and get the affine matrix
                % 2) Get the reverse matrix to register MNI ROIS 2 native
                system(sprintf('flirt -ref %sstandard/MNI152_T1_1mm.nii.gz -in %st1.nii -omat %st12mni.mat -out %st12mni.nii.gz',config.fsl_atlas,config.anat_dir,config.anat_dir,config.anat_dir))
                system(sprintf('convert_xfm -inverse -omat %s/reverse.mat %s/t12mni.mat',config.anat_dir,config.anat_dir))
                
            end
            
            
            
            for rr = 1 : length(rois_fsl)

                
                    % 1) Load the Talaraich atlas
                    % 2) Fin the ROI based on the number it has in atlas
                    % 3) Apply the reverse affine matrix
                    % 4) gunzip the ROI
                    tal = load_nifti(sprintf('%s/atlases/Talairach/Talairach-labels-1mm.nii.gz',config.fsl_atlas));
                    temp = find(tal.vol==rois_fsl_num(rr));
                    template = zeros(size(tal.vol));
                    template(temp) = 1;
                    tal.vol = template;
                    save_nifti(tal,[config.rois_dir sprintf('%s_anat.nii.gz',rois_fsl{rr})]);
                    system(sprintf('flirt -applyxfm -init %sreverse.mat -in %s%s_anat.nii.gz -ref %s -out %s%s_anat.nii.gz',config.anat_dir,config.rois_dir,rois_fsl{rr},config.anat_acpc_freesurf,config.rois_dir,rois_fsl{rr}))
                    system(sprintf('gunzip %s%s_anat.nii.gz',config.rois_dir,rois_fsl{rr}))
                    
          
                    %% Do some dilation end erosion depending on the ROI
                
                if ~isempty(strfind(rois_fsl{rr},'pulv'))
                    
                    system(sprintf('/home/fmridti/Documents/mrtrix3/bin/maskfilter %s%s_anat.nii erode %s%s_anat_dil.nii -npass 1',config.rois_dir,rois_fsl{rr},config.rois_dir,rois_fsl{rr}))
                    system(sprintf('rm %s%s_anat.nii',config.rois_dir,rois_fsl{rr}))
                    system(sprintf('mv %s%s_anat_dil.nii %s%s_anat.nii',config.rois_dir,rois_fsl{rr},config.rois_dir,rois_fsl{rr}))
                    
   
                    
                elseif ~isempty(strfind(rois_fsl{rr},'lgn'))
                    
                    system(sprintf('/home/fmridti/Documents/mrtrix3/bin/maskfilter %s%s_anat.nii dilate %s%s_anat_dil.nii -force -npass 2',config.rois_dir,rois_fsl{rr},config.rois_dir,rois_fsl{rr}))
                    system(sprintf('rm %s%s_anat.nii',config.rois_dir,rois_fsl{rr}))
                    system(sprintf('mv %s%s_anat_dil.nii %s%s_anat.nii',config.rois_dir,rois_fsl{rr},config.rois_dir,rois_fsl{rr}))
         
                end
            end 
        end
    end
end
