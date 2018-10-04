

% This script loads the merged tracktography results (different curvatures)
% and runs LiFe on the concatanated fiber. It also saved cleaned version
% without weights = 0

clc
clear all
project = 'Koulla'

if isempty(strfind(cd,'Volumes')) == 0
    %     [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/Volumes/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'Users')) == 0
    
    
    project_dir = sprintf('/Users/jankurzawski/Dropbox/mrTrix/mrTrix_%s',project);
    
elseif isempty(strfind(cd,'media')) == 0
    
    setenv('PATH','/usr/lib/fsl/5.0/:/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
    addpath(genpath('/home/fmridti/Documents/software'))
    [status,cmdout] = system('echo $HOME')
    project_dir = sprintf('/media/fmridti/fMRI_DTI/DTI/mrTrix/mrTrix_%s',project);
    setenv('PATH','/usr/lib/fsl/5.0/:/usr/local/fsl/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin')
    
    addpath(genpath('/home/fmridti/Documents/software'))
    
end


[status,cmdout] = system('echo $HOME')

cmap = [round([0 189 61]/255,2);round([237 26 55]/255,2);round([151 189 61]/255,2);round([128 0 128]/255,2);round([0 0 255]/255,2)]

subjects = {'101107';'101309';'101410';'101915';'102008';'100206';'100307';'100408';'100610';'101006'}; %done 101006  100610

life = 1;

%%
for 1 = 2:length(subjects)
    tic
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
    fe_rois = {}
    
    %%
    
    load(sprintf('%s/roipairs.mat',project_dir)); % load roipiars
    
    
    
    for l = 1 : length(roipairs)
        
        track{l}(1) =  roipairs(l,1)
        track{l}(2) =  roipairs(l,2)
        
    end
    
    rois = find_pairs(track,roipairs);
    %%
    
    
    parfor r = 1 : size(rois,1) % run in parallel LiFe on selected track
        
        
        allfib = load(sprintf('%sdti/fibers/%s-%s/allfib',config.subject_dir,rois{r,1},rois{r,2}),'allfib')
        allfib = allfib.allfib
        
        %%
        
        if numel(fieldnames(allfib)) > 1 && size(allfib.fibers,1) > 1
            
            
            if exist(sprintf('%s%s-%s.mat',config.subject_dir_life,rois{r,1},rois{r,2}),'file')
                
                
                
            else
                
                
                feFileName = sprintf('%s-%s',rois{r,1},rois{r,2});
                
                fe = feConnectomeInit(dwiFile,allfib,feFileName,[],dwiFileRepeat,t1File);
                fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));
                
                
                fe_rois{r} = fe;
                
                name{r} = sprintf('%s%s-%s.mat',config.subject_dir_life,rois{r,1},rois{r,2})
                fib_name{r} = sprintf('%s%s-%s_afterlife.mat',config.subject_dir_life,rois{r,1},rois{r,2})
                
                
            end
        else
            
            fe_rois{r} = {}
            
        end
    end
    
    
    not_empty = find(~cellfun(@isempty,fe_rois)) % save the clenaed track in a structured and save them
    
    for f = not_empty
        fe = fe_rois{f}
        w  = feGet(fe_rois{f},'fiber weights');
        positive_w_all = w > 0;
        fgl_pos = feGet(fe,'fibers acpc');
        fgl_pos = fgExtract(fgl_pos, positive_w_all, 'keep');
        save(name{f},'fe')
        fgWrite(fgl_pos,fib_name{f},'mat')
        
        
    end
    
    toc
end
