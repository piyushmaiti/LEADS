function generate_tau_pet_multislice_pdf(path_to_ftp, overwrite)
    % Script to create multislice PDFs for the TAU PET images. The multislice PDFs are created within the FTP directory
    %
    % Parameters:
    % -----------
    % path_to_ftp : str
    %     The path to the FTP directory
    % overwrite : bool, optional
    %     Flag to overwrite the file if it already exists. Default is false.
    % Output:
    % -------
    % pdf file saved in the FTP directory containing the multislice PET images

    if nargin < 2
        overwrite = [];
    end

    % Check if the path to the FTP directory exists
    if ~exist(path_to_ftp, 'dir')
        error('The path to the FTP directory does not exist.');
    end
    
    % Add the path to SPM12
    addpath('/home/mac/username/code/matlab/spm12/')
    % Extract the ID
    [~, subject_id, ~] = fileparts(fileparts(path_to_ftp));
    
    % Extract the directory name and date
    [~, dir_name] = fileparts(path_to_ftp);
    scan_date = extractAfter(dir_name, 'FTP_');

    % Creating the full path to the PDF file to be saved to check if it exists
    crid = char(cellstr(subject_id));
    crdate = char(scan_date);
    affine_slovname = fullfile(path_to_ftp, [crid, '_Tau-PET_FTP_', crdate, '.pdf']);
    
    % Check if file exists 
    if exist(affine_slovname, 'file')
        if isempty(overwrite)
            error('Multislice File \n %s already exists. --> Set overwrite flag to "true" to overwrite. \nExiting Code ....', affine_slovname);
        elseif ~overwrite
            error('Multislice File \n %s already exists and overwrite is set to false.', affine_slovname);
        end
        % If overwrite is true, proceed with the processing
    end

    % Display the results
    fprintf('Processing ID: %s for %s\n', subject_id, dir_name);
    
    % Gather the affine suvr image
    fprintf('\nSearching for the affine suvr image and the affine nu image through the MRI folder containing the link....\n');
    affine_suvr_img = char(strcat(path_to_ftp, '/ar', subject_id, '_FTP_', scan_date, '_suvr-infcblgm.nii'));
    fprintf('* affine suvr image found -->  %s\n', affine_suvr_img);

    %% Searching for the affine nu image and the mri folder containing the link
    mri_containing_the_link = char(strcat(path_to_ftp, '/mri'));
    [~, mri_folder] = system(char(strcat('readlink -f',{' '}, mri_containing_the_link)));
    mri_folder = strtrim(mri_folder); % Remove any leading or trailing whitespace
    nu_img = dir(char(strcat(mri_folder,'/',cellstr(subject_id),'_MRI-T1*nu.nii')));
    

    % Gather the affine nu image
    affine_nu_img = char(strcat(mri_folder,'/a',nu_img.name));
    fprintf('* affine nu image found -->  %s\n', affine_nu_img);
    
    % affine warp, NIH with auto range
    fprintf('\n Creating the multislice PDF for the TAU PET images....\n');
    tempimgvals = spm_read_vols(spm_vol(affine_suvr_img));
    tempimgvals = nonzeros(tempimgvals);
    tempimgvals = tempimgvals(~isnan(tempimgvals));
    rangecolorscale = [min(tempimgvals) max(tempimgvals)];
    
    o = slover;
    o.cbar = 2;
    o.img(1).vol = spm_vol(affine_nu_img);
    o.img(1).type = 'structural';
    o.img(1).prop = 1;
    o.img(2).vol = spm_vol(affine_suvr_img);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = 'nih.lut';
    o.img(2).range = [0.5 rangecolorscale(2)];
    o.img(2).prop = 0.7;
    o.transform = 'axial';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults(o);
    o.slices = -30:6:58;
    o = paint(o);
    
    crdate = char(scan_date);
    crid = char(cellstr(subject_id));
    jpeglab = strcat('ID:',{' '},crid,{' '},'***',{' '},'FTP-PET',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.98 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab2 = strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Generated Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.96 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab3 = strcat('SUVR Map - Ref region: Inferior Cerebellar GM',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.94 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab4 = strcat('PET resolution is now 6mm FWHM; direct comparison with older data (8mm) is not recommended.');
    hTitAx4 = axes('Parent',o.figure,'Position',[0 0.92 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab4,'Parent',hTitAx4,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','red','FontSize',9);
    
    % Save the file
    print(affine_slovname, '-dpdf', '-r300');
    fprintf('File saved: %s\n', affine_slovname);
    
end

