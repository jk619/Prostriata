Processing Pipeline for Prostriata ROI2ROI tractography 
------------------------------------------------
1. do0_hcp_rois_koulla.m

This maps the cortical ROIs from Glasser template and saves them as a nifti volume in the native space of the subjects. It also uses the Talaraich Deamon Atlas for subcortical ROIs (i.e LGN). Registers the anatomy to the MNI template to obtain the affine registration matrix. The registration mat is inversed to allow mapping of the template ROIs to the native space
The end result is a set of ROIs (nifti) alligned to the native space

------------------------------------------------


2. do1_preproc_and_track_hcp_koulla.m 
This does the ac-pc alignemnt, preprocessing, dtiInit, mapping ROIs to acpc space, ensemble tractography with different curvatures

The end result is a set of tracts (5 files for 5 curvatures) for each ROI

-----------------------------------------------

3. do2_edit_wm_masks.m 
This loads anatomy + wm for manual correction
The end results are the corrected WM masks further used to constrain tractography

-----------------------------------------------


4. do3_merge_allfib_koulla.m 
This merges the seperate tracts with different curvature thresholds to one big file
The end results is one bigfile per pair of ROIS containing all discovered tracts

-----------------------------------------------


5. do4_life_on_merged.m 
This runs LiFe on merged fibers and saves a clean version.
-----------------------------------------------


The end results is one clean tract per set of ROIs

