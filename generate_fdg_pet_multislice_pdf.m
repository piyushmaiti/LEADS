function generate_fdg_pet_multislice_pdf(path_to_fdg, overwrite)
    % Script to create multislice PDFs for the FDG PET images. The multislice PDFs are created within the FDG directory
    %
    % Parameters:
    % -----------
    % path_to_fdg : str
    %     The path to the FDG directory
    % overwrite : bool, optional
    %     Flag to overwrite the file if it already exists. Default is false.
    % Output:
    % -------
    % pdf file saved in the FDG directory containing the multislice PET images

    if nargin < 2
        overwrite = [];
    end

    % Check if the path to the FDG directory exists
    if ~exist(path_to_fdg, 'dir')
        error('The path to the FDG directory does not exist.');
    end
    
    % Extract the ID
    [~, subject_id, ~] = fileparts(fileparts(path_to_fdg));
    
    % Extract the directory name and date
    [~, dir_name] = fileparts(path_to_fdg);
    scan_date = extractAfter(dir_name, 'FDG_');

    % Creating the full path to the PDF file to be saved to check if it exists
    crid = char(cellstr(subject_id));
    crdate = char(scan_date);
    affine_slovname = fullfile(path_to_fdg, [crid, '_FDG-PET_', crdate, '.pdf']);
    
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
    affine_suvr_img = char(strcat(path_to_fdg, '/ar', subject_id, '_FDG_', scan_date, '_suvr-pons.nii'));
    fprintf('* affine suvr image found -->  %s\n', affine_suvr_img);

    %% Searching for the affine nu image and the mri folder containing the link
    mri_containing_the_link = char(strcat(path_to_fdg, '/mri'));
    [~, mri_folder] = system(char(strcat('readlink -f',{' '}, mri_containing_the_link)));
    mri_folder = strtrim(mri_folder); % Remove any leading or trailing whitespace
    nu_img = dir(char(strcat(mri_folder,'/',cellstr(subject_id),'_MRI-T1*nu.nii')));
    

    % Gather the affine nu image
    affine_nu_img = char(strcat(mri_folder,'/a',nu_img.name));
    fprintf('* affine nu image found -->  %s\n', affine_nu_img);
    
    % affine warp, NIH with auto range
    fprintf('\n Creating the multislice PDF for the FDG PET images....\n');
    tempimgvals = spm_read_vols(spm_vol(affine_suvr_img));
    tempimgvals = nonzeros(tempimgvals);
    tempimgvals = tempimgvals(~isnan(tempimgvals));
    rangecolorscale = [min(tempimgvals) max(tempimgvals)];
    
    affine_slovname1=fullfile(path_to_fdg, [crid, '_FDG-PET_axial', crdate, '.pdf']);
    o = slover;
    o.cbar = 2;
    o.img(1).vol = spm_vol(affine_nu_img);
    o.img(1).type = 'structural';
    o.img(1).prop = 1;
    o.img(2).vol = spm_vol(affine_suvr_img);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = 'nih.lut';
    o.img(2).range = [0.1 2.2];
    o.img(2).prop = 0.7;
    o.transform = 'axial';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults(o);
    o.slices = -30:6:58;
    o = paint(o);
    
    crdate = char(scan_date);
    crid = char(cellstr(subject_id));
    jpeglab = strcat('ID:',{' '},crid,{' '},'***',{' '},'FDG-PET',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.98 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab2 = strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Generated Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.96 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab3 = strcat('SUVR Map - Ref region: Pons',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.94 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab4 = strcat('PET resolution is now 6mm FWHM; direct comparison with older data (8mm) is not recommended.');
    hTitAx4 = axes('Parent',o.figure,'Position',[0 0.92 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab4,'Parent',hTitAx4,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','red','FontSize',9);
    print(affine_slovname1, '-dpdf', '-r300');
    
    % affine warp with NIH - coronal
    affine_slovname2=fullfile(path_to_fdg, [crid, '_FDG-PET_coronal', crdate, '.pdf']);

    o = slover;
    o.cbar = 2;
    o.img(1).vol=spm_vol(affine_nu_img);
    o.img(1).type='structural';
    o.img(1).prop=1;
    o.img(2).vol = spm_vol(affine_suvr_img);
    o.img(2).type = 'truecolour';
    o.img(2).cmap = 'nih.lut';
    o.img(2).range = [0.1 2.2];
    o.img(2).prop=0.7;
    o.transform = 'coronal';
    o.figure = spm_figure('GetWin','Graphics');
    o = fill_defaults (o);
    o.slices = -65:8:50;
    o = paint(o);

    crdate = char(scan_date);
    crid = char(cellstr(subject_id));
    jpeglab = strcat('ID:',{' '},crid,{' '},'***',{' '},'FDG-PET',{' '},'***',{' '},'LEADS.PETCORE@ucsf.edu');
    hTitAx = axes('Parent',o.figure,'Position',[0 0.98 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab,'Parent',hTitAx,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab2 = strcat('Scan Date:',{' '},datestr(crdate),{' '},'***',{' '},'Generated Date:',{' '},date);
    hTitAx2 = axes('Parent',o.figure,'Position',[0 0.96 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab2,'Parent',hTitAx2,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab3 = strcat('SUVR Map - Ref region: Pons',{' '},'***',{' '},'Left is Left');
    hTitAx3 = axes('Parent',o.figure,'Position',[0 0.94 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab3,'Parent',hTitAx3,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','black','FontSize',10);
    jpeglab4 = strcat('PET resolution is now 6mm FWHM; direct comparison with older data (8mm) is not recommended.');
    hTitAx4 = axes('Parent',o.figure,'Position',[0 0.92 0.06 0.02],'Visible','off');
    text(0.5,0,jpeglab4,'Parent',hTitAx4,'HorizontalAlignment','left','VerticalAlignment','baseline','Color','red','FontSize',9);
    print(affine_slovname2,'-dpdf','-r300');

    cmdmerge=strcat('pdfunite', {' '},  affine_slovname1, {' '}, affine_slovname2, {' '}, affine_slovname);
    system(char(cmdmerge));
    fprintf('File saved: %s\n', affine_slovname);

    % Remove the temporary files
    delete(affine_slovname1);
    delete(affine_slovname2);
    
end
