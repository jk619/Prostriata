%% This script combines results from tracking with different curvatures together.

clc
clear all
project = 'Koulla'
addpath('./notuse')

if isempty(strfind(cd,'Volumes')) == 0
    %     [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/Volumes/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'Users')) == 0
    
    
    project_dir = sprintf('/Users/jankurzawski/Dropbox/mrTrix/mrTrix_%s',project);
    
elseif ~contains(cd,'media') == 0
    
    setenv('PATH','/usr/lib/fsl/5.0/:/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
    addpath(genpath('/home/fmridti/Documents/software'))
    [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/media/fmridti/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    setenv('PATH','/usr/lib/fsl/5.0/:/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
    
    addpath(genpath('/home/fmridti/Documents/software'))
    
end


[status,cmdout] = system('echo $HOME')

subjects = {'101107';'101309';'101410';'101915';'102008';'100206';'100307';'100408';'100610';'101006'}; %done 101006  100610

%%
for s = 1:length(subjects)
    %setup config
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life = [config.subject_dir 'life/'];
    config.anat = sprintf('%s%s_t1.nii',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    config.anat_acpc = sprintf('%st1/t1_acpc2.nii',config.subject_dir);
    dwiFile = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    dwiFileRepeat = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    t1File = config.anat_acpc;
    Niter = 500;
    
    
    
    
    load(sprintf('%s/roipairs.mat',project_dir)); % here we have all the pairs of ROIs saved
    
    
    
    for l = 1 : length(roipairs)
        
        track{l}(1) =  roipairs(l,1)
        track{l}(2) =  roipairs(l,2)
        
    end
    rois = find_pairs(track,roipairs);
    rois = rois(end-4:end,:);
    
    %%
    
    for r = 1 : size(rois,1)
        clear total
        
        c = 1;
        
        fgl_sep = dir(sprintf('%sdti/fibers/%s-%s/*.tck',config.subject_dir,rois{r,1},rois{r,2}));
        %         fgl_sep = fgl_sep(3:5:end)
        for rr = 1 : length(fgl_sep)
            
            temp =  dtiImportFibersMrtrix(sprintf('%sdti/fibers/%s-%s/%s',config.subject_dir,rois{r,1},rois{r,2},fgl_sep(rr).name));
            
            if isempty(temp.fibers{1})
                
                disp('empty')
                
            else
                
                total{c} = temp;
                
                c = c+1;
                
            end
            
        end
        
        if ~exist('total','var')
            total{1}.name = Inf;
            
        end
        %%
        
        
        for len = 1 : length(total)
            
            if len == 1 && sum(isfinite(total{1}.name)) > 0
                
                allfib  = total{len};
                
            elseif len == 1 &&  sum(isfinite(total{1}.name)) == 0
                
                allfib.fibers = {};
                
                
            else
                
                allfib = fgMerge(allfib,total{len});
                
                
                
            end
            
        end
        save(sprintf('%s%s-%s/allfib',config.fibers_dir,rois{r,1},rois{r,2}),'allfib')
        fgWrite(allfib,sprintf('%s%s-%s/allfib_mr',config.fibers_dir,rois{r,1},rois{r,2}),'mat')
        save(sprintf('%s%s-%s/config',config.fibers_dir,rois{r,1},rois{r,2}),'config')
        
        
        %% save combined ROIS as one file allfib, allfib_mr (readable by mrDiffusion)
    end
end

