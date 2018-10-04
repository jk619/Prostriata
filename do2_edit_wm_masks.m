%% This script just opens the anatomy T1 and corresponding WM mask for manual correction
clear all
close all

project = 'Koulla'
addpath('./notuse')

if isempty(strfind(cd,'Volumes')) == 0
    %     [status,cmdout] = system('echo $HOME')
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

% cmap = {round([151 189 61]/255,2),round([237 26 55]/255,2),round([51 51 255]/255,2)}
cmap = [round([0 189 61]/255,2);round([237 26 55]/255,2);round([151 189 61]/255,2);round([128 0 128]/255,2);round([0 0 255]/255,2)]




%subjects = {'101107';'101309';'101410';'101915';'102008';'100206';'100307';'100408';'100610';'101006'}; %done 101006  100610
subjects = {'100307'}



for s = 1 : length(subjects)
    
    
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
    wm_freesurfer_mif = sprintf('%sdti/bin/wm_freesurfer.mif', config.subject_dir)
    delete(wm_freesurfer_mif)
    wm_freesurfer_nii = sprintf('%sdti/bin/wm_freesurfer.nii.gz', config.subject_dir)
    system(sprintf('fsleyes %s %s &',wm_freesurfer_nii,config.anat_acpc))
    
    
end

