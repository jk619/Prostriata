% Build the file names for the diffusion data, the anatomical MRI.
cls
project_dir = '/Users/jankurzawski/Dropbox/mrTrix_Prostriata';
subjects = {'pro_AB';'pro_CL';'pro_EG';'pro_FM';'pro_GA';'pro_LV';'pro_MC';'pro_IS'};

%%
for s = 5
    % set up dir + files
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life = [config.subject_dir 'life'];
    config.anat = sprintf('%s%s_t1.nii',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    dwiFile = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    dwiFileRepeat = sprintf('%sdti/%s_dwi_aligned_trilin_noMEC.nii.gz',config.subject_dir,config.subject);
    t1File = sprintf('%st1/t1_acpc.nii',config.subject_dir);
    
    % load all roipairs
    load(sprintf('%sroipairs.mat',config.subject_dir));
    
    %set rois
    
    roi1 = roipairs{1,1};
    roi2 = roipairs{1,2};
    
    sprintf('Track %s-%s',roi1,roi2)
    if  exist(sprintf('%s/%s-%s/fg_life.mat',config.subject_dir_life,roi1,roi2),'file') ~=2 || exist(sprintf('%s/%s-%s/fe_life.mat',config.subject_dir_life,roi1,roi2),'file') ~=2
        
        % load wholebrain tck
        wb = sprintf('%sdti/fibers/wholebrain/%s_dwi_aligned_trilin_noMEC_csd_lmax4_%s_dwi_aligned_trilin_noMEC_wm_wmMask_prob.tck',config.subject_dir,subjects{s},subjects{s});
        % make wholebrain tck excluding roi2roi and save
        system(sprintf('filter_tracks -invert %s -include %s%s_anat.mif -include %s%s_anat.mif %s/%s-%s/connectome_ex_roi2roi.tck',wb,config.rois_dir,roi1,config.rois_dir,roi2,config.fibers_dir,roi1,roi2))
        cut_wb = dtiImportFibersMrtrix(sprintf('%s%s-%s/connectome_ex_roi2roi.tck',config.fibers_dir,roi1,roi2));
        fgWrite(cut_wb,sprintf('%s%s-%s/connectome_ex_roi2roi.mat',config.fibers_dir,roi1,roi2));
        
        
        % combine wholebrain without tck with roi2roi tck
        roi2roi = dtiImportFibersMrtrix(sprintf('/Users/jankurzawski/Dropbox/mrTrix_Prostriata/pro_GA/dti/fibers/%s-%s/pro_GA_dwi_aligned_trilin_noMEC_csd_lmax4_%s_anat_%s_anat_%s-%s_union_wmMask_prob.tck',roi1,roi2,roi1,roi2,roi1,roi2));
        
        wb_roi2roi = fgMerge(roi2roi,cut_wb);
        wb_roi2roi_name = sprintf('%s%s-%s/connectome_in_roi2roi.mat',config.fibers_dir,roi1,roi2);
        fgWrite(wb_roi2roi,wb_roi2roi_name);
        
        feFileName    = sprintf('%s-%s',roi1,roi2);
        
        % run life on combined
        
        mkdir(sprintf('%s/%s-%s/',config.subject_dir_life,roi1,roi2))
        
        fe = feConnectomeInit(dwiFile,wb_roi2roi_name,feFileName,[],dwiFileRepeat,t1File);
        
        
        Niter = 500;
        % Fit LiFE model
        fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));
        %     fe = feSet(fe,'fit',feFitModel(feGet(fe,'mfiber'),feGet(fe,'dsigdemeaned'),'bbnnls'));
        % Remove non zero weighted fascicles
        fg = feGet(fe,'fibers acpc');
        w = feGet(fe,'fiber weights');
        positive_w = w > 0;
        fg = fgExtract(fg, positive_w, 'keep');
        
        save(sprintf('%s/%s-%s/fg_life.mat',config.subject_dir_life,roi1,roi2),'fg')
        save(sprintf('%s/%s-%s/fe_life.mat',config.subject_dir_life,roi1,roi2),'fe')
        %
        %
    else
        %
        load(sprintf('%s/%s-%s/fg_life.mat',config.subject_dir_life,roi1,roi2))
        load(sprintf('%s/%s-%s/fe_life.mat',config.subject_dir_life,roi1,roi2))
        %
    end
    return
    %%
    
    %     config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir)
    %
    %     rois = subdir([config.rois_dir, '*.mat']);
    %     roi1 = dtiReadRoi(rois(1).name);
    %     roi2 = dtiReadRoi(rois(2).name);
    %     [fg_rois fascicleIndices] = feSegmentFascicleFromConnectome(fg, {rois(1).name,rois(2).name},{'and','and'},'v1-lgn')
    
    [fglClean, fibersToKeep] = mbaComputeFibersOutliers(roi2roi, 3, 3);
    fascicleIndices  = (1:length(fglClean.fibers));
    
    
    % set up boolean variable inputs for virtual lesion function
    display.tract = true;
    display.distributions = true;
    display.evidence = true;
    
    [se, feLesion,  newfe, newFascicleIndices, indicesFibersKept, commonCoords] = feVirtualLesion(fe, fascicleIndices, display);
    
    %
    
    % [sr_el_fg, sr_el_cmatrix] = dtiConnectMultipleROIs([], 'and', [], [roi1 roi2], fg)
    
    %%
    % Use MBA to visualize
    colors     = {[.75 .25 .1]};
    viewCoords = [90,0];
    proportion_to_show = .05;
    threshold_length   = 5;
    slices     = {[10 0 0],[0 -40 0],[0 0 -14]}; % Values in mm from AC
    
    % Prepare the plot of tract
    %     fg_rois = rmfield(fg_rois,'coordspace');
    
    % Pick a percentage of fascicles to display (the PN can be too dense for visualization).
    %fibs_indx = RandSample(1:length(fg.fibers),round(length(fg.fibers)*proportion_to_show));
    %     fg.fibers = fg.fibers(1:2000);
    
    % Display fascicles and anatomy
    fig_h = figure('name','Whole Brain','color','k');
    hold on
    
    % display anatomy
    t1 = niftiRead(t1File);
    %     h  = mbaDisplayBrainSlice(t1, slices{2});
    h  = mbaDisplayBrainSlice(t1, slices{1});
    
    % Disdplay fasciles
    set(gca,'visible','off','ylim',[-108 69],'xlim',[-75 75],'zlim',[-45 78],'Color','w')
    [~, light_h] = mbaDisplayConnectome(fg_rois.fibers,fig_h,colors{1},'single');
    delete(light_h)
    view(viewCoords(1),viewCoords(2))
    light_h = camlight('right');
    lighting phong;
    %set(fig_h,'Units','normalized', 'Position',[0.5 .2 .4 .8]);
    set(gcf,'Color',[1 1 1])
    drawnow
    
    %     quit
    
    % save figure to disk
    print(fig_h,'my_figure','-dpng','-r600');%'-dpsc2')
end