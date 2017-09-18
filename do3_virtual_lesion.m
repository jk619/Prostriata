cls
project_dir = '/Users/jankurzawski/Dropbox/mrTrix_Prostriata';
subjects = {'pro_AB';'pro_CL';'pro_EG';'pro_FM';'pro_GA';'pro_LV';'pro_MC';'pro_IS'};

draw = 0;
std_clean = 0;
weight_clean = 1;
vl = 1;
%%
for s = 5
    %setup config
    config.subject = subjects{s};
    config.subject_dir = sprintf('%s/%s/',project_dir,config.subject);
    config.subject_dir_life = [config.subject_dir 'life'];
    config.anat = sprintf('%s%s_t1.nii',config.subject_dir,config.subject);
    config.dwi = sprintf('%s%s_dwi.nii.gz',config.subject_dir,config.subject);
    config.bvals = sprintf('%s%s_dwi.bvals',config.subject_dir,config.subject);
    config.bvecs = sprintf('%s%s_dwi.bvecs',config.subject_dir,config.subject);
    config.rois_dir = sprintf('%sdti/ROIs/anat/',config.subject_dir);
    config.fibers_dir = sprintf('%sdti/fibers/',config.subject_dir);
    config.anat_acpc = sprintf('%st1/t1_acpc.nii',config.subject_dir);
    
    
    % load all roipairs
    load(sprintf('%sroipairs.mat',config.subject_dir));
    
    %set rois
    
    roi1 = roipairs{1,1};
    roi2 = roipairs{1,2};
    
    %load fe
    load(sprintf('%s/%s-%s/fe_life.mat',config.subject_dir_life,roi1,roi2))
    load(sprintf('%s/%s-%s/fg_life.mat',config.subject_dir_life,roi1,roi2))

    
    %%
    %load roi2roi
    
    %     ff2 = roipairs{8,1};
    %     fend2 = roipairs{8,2};
    %     fgl = fgRead(sprintf('/Users/jankurzawski/Dropbox/mrTrix_Prostriata/pro_GA/dti/fibers/pros_lh/pro_GA_dwi_aligned_trilin_noMEC_csd_lmax4_pros_lh_anat_wmMask_prob.tck');
    
    fgl = fgRead(sprintf('/Users/jankurzawski/Dropbox/mrTrix_Prostriata/pro_GA/dti/fibers/%s-%s/pro_GA_dwi_aligned_trilin_noMEC_csd_lmax4_%s_anat_%s_anat_%s-%s_union_wmMask_prob.tck',roi1,roi2,roi1,roi2,roi1,roi2));
    %     fgl2 = fgRead(sprintf('/Users/jankurzawski/Dropbox/mrTrix_Prostriata/pro_GA/dti/fibers/%s-%s/pro_GA_dwi_aligned_trilin_noMEC_csd_lmax4_%s_anat_%s_anat_%s-%s_union_wmMask_prob.tck',ff2,fend2,ff2,fend2,ff2,fend2))
    
    %clean roi2roi
    
    if weight_clean
        w = feGet(fe,'fiber weights');
        w = w(1:length(fgl.fibers));
        positive_w = w > 0;
        fglClean = fgExtract(fgl, positive_w, 'keep');
        
        
    elseif std_clean
        [fglClean, fibersToKeep] = mbaComputeFibersOutliers(fgl, 5,2);
        
    else
        
        fglClean = fgl;
        
    end
    %     [fglClean2, fibersToKeep] = mbaComputeFibersOutliers(fgl2, 3,3);
    
    
    
    if draw
        draw_fg(fglClean,config.anat_acpc)
    end
    
    %%
  if vl
    %%
    % Save indices of the roi2roi
    fascicleIndices  = (1:length(fglClean.fibers));
    %     fascicleIndices2  = (1:length(fglClean2.fibers));
    
    % Setup display
    display.tract = false;
    display.distributions = true;
    display.evidence = true;
    
    
   
    %clean wb connectome 
    
%     w_all = feGet(fe,'fiber weights');
%     positive_w_all = w_all > 0;
%     fg = feGet(fe,'fibers acpc');
%     fg = fgExtract(fg, positive_w_all, 'keep');
%     fe.fg = fg;
%     
    
     % run Virtual lesion
    [se, feLesion,  feNoLesion, newFascicleIndices, indicesFibersKept, commonCoords] = ...
        feVirtualLesion(fe, fascicleIndices, display);
    
    
  end
  return
    
    %% AFQ
    fg_clean = fglClean;
    numNodes = 100;
    dt = dtiLoadDt6(sprintf('%s/dti/dt6.mat',config.subject_dir));
    
    [fa md rd ad] = AFQ_ComputeTractProperties(fg_clean, dt, numNodes,0);
    
    %%
    
    % Open a new figure window
    figure; hold('on');
    % Set the coloring of each plot
    set(gca,'ColorOrder',jet(20));
    % Plot the Tract FA Profiles for each tract
    plot(fa,'linewidth',2);
    xlabel('Node');
    ylabel('Fractional Anisotropy');
    title('Tract Profiles');
    % Add a legend with the names of each fiber tract. Here are the names of
    % each fiber group.
    fgNames={'Optic Radiation'};
    legend(fgNames,'Location','EastOutside' );
    
    
    %%
    % Render the Tract FA Profile for the left arcuate fasciculus. When the
    % argument 'dt' is passed in follwed by a variable containing the dt6
    % structure then the tract profile is added to the plot. The colormap
    % denotes the FA value at each point along the tract core.
    AFQ_RenderFibers(fg_clean,'dt',dt);

end