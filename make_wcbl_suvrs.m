clear
clc

% CHOOSE A TEST SUBJECT
subj = 'LDS1770351';
tp = 'Timepoint1';
tracer = 'FDG';
acq_date = '2021-05-07';

% ------------------------------
leads_proc_dir = '/mnt/coredata/Projects/LEADS/data_f7p1/processed';
subj_fdg_dir = strcat('/', subj, '/', tp, '/', tracer, '_', acq_date);
fdg_proc_dir = {strcat(leads_proc_dir, '/', subj_fdg_dir)};

for i = 1:size(fdg_proc_dir, 1)
    % Load the files that we need to work with

    % Native aparc
    aparc_f = dir(strcat(fdg_proc_dir{i,1},'/LDS*raparc+aseg.nii'));
    aparc_f = strcat(aparc_f.folder, '/', aparc_f.name);
    aparc = spm_read_vols(spm_vol(aparc_f));

    % Native FDR SUVR (Pons)
    fdg_pons_f = dir(strcat(fdg_proc_dir{i,1},'/LDS*_suvr_pons.nii'));
    fdg_pons_f = strcat(fdg_pons_f.folder,'/',fdg_pons_f.name);
    fdg_pons = spm_read_vols(spm_vol(fdg_pons_f));

    % Warped FDR SUVR (Pons)
    wfdg_pons_f=dir(strcat(fdg_proc_dir{i,1},'/wLDS*_suvr_pons.nii'));
    wfdg_pons_f=strcat(wfdg_pons_f.folder,'/',wfdg_pons_f.name);
    wfdg_pons = spm_read_vols(spm_vol(wfdg_pons_f));

    % Store the file names we need to create

    % Native FDG SUVR (Whole Cerebellum)
    fdg_wcbl_f = strcat(fdg_proc_dir{i,1}, '/', subj, '_', tracer, '_', acq_date, '_suvr_wcbl.nii');

    % Warped FDG SUVR (Whole Cerebellum)
    wfdg_wcbl_f = strcat(fdg_proc_dir{i,1}, '/w', subj, '_', tracer, '_', acq_date, '_suvr_wcbl.nii');

    % Create the native space whole cerebellum mask and save it to disc
    wcbl_mask_f = strcat(fdg_proc_dir{i,1},'/wcbl_ref_mask.nii');
    wcbl_idx = [8; 47; 7; 46];
    wcbl_mask = zeros(size(aparc));
    for i = 1:length(wcbl_idx)
        idx = wcbl_idx(i);
        wcbl_mask = wcbl_mask + (aparc == idx);
    end
    V = spm_vol(aparc_f);
    V.dt = [spm_type('uint8'), 0];
    V.pinfo = [1, 0, 0]';
    V.fname = wcbl_mask_f;
    spm_write_vol(V, wcbl_mask);
    clear V
  
    % Calculate and save the whole cerebellum FDG SUVRs to disc
    wcbl_mean = mean(fdg_pons(logical(wcbl_mask)));
    
    fdg_wcbl = fdg_pons ./ wcbl_mean;
    V = spm_vol(fdg_pons_f);
    V.dt = [16, 0]; % 16 = float32; saving the data as float32
    V.fname = fdg_wcbl_f;
    spm_write_vol(V, fdg_wcbl);
    clear V

    wfdg_wcbl = wfdg_pons / wcbl_mean;
    V = spm_vol(wfdg_pons_f);
    V.dt = [16, 0]; % 16 = float32; saving the data as float32
    V.fname = wfdg_wcbl_f;
    spm_write_vol(V, wfdg_wcbl);
    clear V

end
