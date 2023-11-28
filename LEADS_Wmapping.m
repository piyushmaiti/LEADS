% This script is used to generate the Wmaps for the LEADS data
% Updated - Novemeber 20, 2023
% This is the version 2 of the Wmapping script - This version uses the updated calibration data
% Author - Piyush Maiti
% Email - piyush.maiti@ucsf.edu
% ------------------------------------------------- %

path_processed='/mnt/coredata/Projects/LEADS/data_f7p1/processed/'; % LEADS Data Directory

wscore_path='/mnt/coredata/Projects/LEADS/script_f7p1/wscoring'; % Path to the folder containing the wscoring scripts
path_wscore_meta = fullfile(wscore_path,'meta'); % Path to the meta
wscore_calibration_path=fullfile(wscore_path, 'calibration_61AmyNegCN-v2'); % Change the Folder to the calibration data based on the version of the calibration data you want to use
fprintf(1, "For more information about the calibration data, please refer to the README.md file in the calibration folder %s\n", wscore_calibration_path);

% MRI calibration Locations and Variables
itcpt_mri = fullfile(wscore_calibration_path, 'mri', 'beta_0001.nii');
betaage_mri = fullfile(wscore_calibration_path, 'mri', 'beta_0002.nii');
betativ_mri =  fullfile(wscore_calibration_path, 'mri', 'beta_0003.nii'); % TIV : Total Intracranial Volume - Used for Wsocring gray matter probability maps
sdmap_mri = fullfile(wscore_calibration_path, 'mri', 'SDmap.nii');

% PET calibration Locations and Variables for FBB and FTP
calibration_modals = {'fbb', 'ftp'};

for i = 1:length(calibration_modals)
    eval(['itcpt_' calibration_modals{i} ' = fullfile(wscore_calibration_path, calibration_modals{i}, ''beta_0001.nii'');']);
    eval(['betaage_' calibration_modals{i} ' = fullfile(wscore_calibration_path, calibration_modals{i}, ''beta_0002.nii'');']); 
    eval(['sdmap_' calibration_modals{i} ' = fullfile(wscore_calibration_path, calibration_modals{i}, ''SDmap.nii'');']);
end

% PET calibration Locations and Variables for FDG
fdg_calibration = {'fdg_pons', 'fdg_whlcerebellum'};

for i = 1:length(fdg_calibration)
    eval(['itcpt_' fdg_calibration{i} ' = fullfile(wscore_calibration_path, fdg_calibration{i}, ''beta_0001.nii'');']);
    eval(['betaage_' fdg_calibration{i} ' = fullfile(wscore_calibration_path, fdg_calibration{i}, ''beta_0002.nii'');']);
    eval(['betasex_' fdg_calibration{i} ' = fullfile(wscore_calibration_path, fdg_calibration{i}, ''beta_0003.nii'');']) 
    eval(['sdmap_' fdg_calibration{i} ' = fullfile(wscore_calibration_path, fdg_calibration{i}, ''SDmap.nii'');']);
end

emask='/mnt/coredata/Projects/LEADS/script_f7p1/templates/EM_final.nii';

% Check if the calibration files exist
for i = 1:length(calibration_modals)
    if ~isfile(eval(['itcpt_' calibration_modals{i} ]))
        fprintf('File %s not found.\n', eval(['itcpt_' calibration_modals{i} ]));
    end
    if ~isfile(eval(['betaage_' calibration_modals{i} ]))
        fprintf('File %s not found.\n', eval(['betaage_' calibration_modals{i} ]));
    end
    if ~isfile(eval(['sdmap_' calibration_modals{i} ]))
        fprintf('File %s not found.\n', eval(['sdmap_' calibration_modals{i} ]));
    end
end
fprintf(1,"All calibration files found.\n");
% ------------------------------------------------- % 
%% Reminder Message to the user of the approach

fprintf(1,'The Wmapping script now looks for the LEADS_PTDEMOG.csv file in the petcore drive,\nin /shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/\nPress enter to acknowledge and continue:\n');
pause;

% new approach, look for the demographic database in the shared/petcore drive, keeping the approach to grab the latest if multiple are present

dbage=dir('/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/LEADS_PTDEMOG*');
dbage=struct2cell(dbage)';
dbage=dbage(:,[1 3]);
dbage=array2table(dbage);
dbage.Properties.VariableNames(1) = cellstr(strcat('Filename'));
dbage.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
dbage.DateCreated=datetime(dbage.DateCreated);
filt=max(dbage.DateCreated); % Store most recent date
dbage=dbage(dbage.DateCreated==filt,:); % Subset to select the most recent file 
clear filt;

% Creating a table to store the Date of Birth of the participants
dbagef=readtable(strcat('/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/',char(dbage{1,1}))); %% ready to check differences now
dbagef=dbagef(:,[1 7]);
dbagef.Properties.VariableNames(1) = cellstr(strcat('ID'));
dbagef.Properties.VariableNames(2) = cellstr(strcat('dob'));
dbagef.ID=cellfun(@upper,dbagef.ID,'UniformOutput',false);
dbagef.dob=datetime(dbagef.dob);

% Create a table to store the Sex of the participants
dbsex=readtable(strcat('/shared/petcore/Projects/LEADS/data_f7p1/LONI_uploads/service/',char(dbage{1,1}))); %% ready to check differences now
dbsex=dbsex(:,[1 6]);
dbsex.Properties.VariableNames(1) = cellstr(strcat('ID'));
dbsex.Properties.VariableNames(2) = cellstr(strcat('sex'));
dbsex.ID=cellfun(@upper,dbsex.ID,'UniformOutput',false);
dbsex.sex = cellfun(@str2double, dbsex.sex);
dbsex.sex(dbsex.sex == 2) = 0;

% Get the list of available seg8.mat files
srcseg8s=dir(strcat(path_processed,'LDS*/*/MRI*/*seg8.mat'));
allseg8s = strcat({srcseg8s.folder}','/',{srcseg8s.name}');

% Prepare for later
allseg8s=array2table(horzcat(allseg8s, regexp(allseg8s, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once')));
allseg8s.Properties.VariableNames(1) = cellstr(strcat('File'));
allseg8s.Properties.VariableNames(2) = cellstr(strcat('Identifier'));
allseg8s=sortrows(allseg8s,'Identifier');

fprintf("Database with MRI seg8.mat ready, Count: %d\n", size(allseg8s,1));

% ------------------------------------------------- % 
olddbtiv= dir(fullfile(path_wscore_meta,'LEADS_TIV*'));
olddbtiv=struct2cell(olddbtiv)';
olddbtiv=olddbtiv(:,[1 3]);
olddbtiv=array2table(olddbtiv);
olddbtiv.Properties.VariableNames(1) = cellstr(strcat('Filename'));
olddbtiv.Properties.VariableNames(2) = cellstr(strcat('DateCreated'));
olddbtiv.DateCreated=datetime(olddbtiv.DateCreated);
filt=max(olddbtiv.DateCreated); % Store most recent date
olddbtiv=olddbtiv(olddbtiv.DateCreated==filt,:); % Subset to select the most recent file
oldinfotiv=readtable(fullfile(path_wscore_meta,char(table2cell(olddbtiv(1,1)))),'Delimiter',','); %% ready to check differences now
oldinfotiv=sortrows(oldinfotiv,'Identifier');

% Check new MRIs need to be processed and added to the database
% Comparing the list of nu_seg8.mat files with the database of TIVs and searching for nu_seg8.mat files using the 'Identifier' column.
% Example Identifier : LDS9410679_MRI_T1_2023-09-08

fprintf(1,'Checking if new MRIs needs to be processed:  \n\n'); 

checkinfo=isequal(allseg8s(:,{'Identifier'}),oldinfotiv(:,{'Identifier'})); % Check for overlap between identifiers 

if checkinfo==0
    
    Indexa=not(ismember(allseg8s(:,{'Identifier'}),oldinfotiv(:,{'Identifier'}))); %% logical indexing for the new cases. 
    newcases=allseg8s(Indexa,:); %% Creates a table with the New MRI seg8 scans
    
    if size(newcases,1)>0

        fprintf(1,'**New MRI seg8 files were found:\n\n');
        disp(newcases.Identifier);
        newcases=table2array(newcases);
        newcases=newcases(:,1); % prepare for SPM

        clear matlabbatch
        spm('defaults','PET');
        matlabbatch{1}.spm.util.tvol.matfiles = newcases;
        matlabbatch{1}.spm.util.tvol.tmax = 3;
        matlabbatch{1}.spm.util.tvol.mask = {'/mnt/neuroimaging/SPM/spm12/tpm/mask_ICV.nii,1'};
        matlabbatch{1}.spm.util.tvol.outf = char(strcat(path_wscore_meta,'tempTIV.csv'));
        spm_jobman('run',matlabbatch); clear matlabbatch;
        
        % read the new values, compute the TIV
        dbtiv=readtable(char(fullfile(path_wscore_meta,'tempTIV.csv')),'Delimiter',',');
        dbtiv.TIV=sum(dbtiv{:,2:end},2);
        dbtiv.Identifier=regexp(dbtiv.File, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once');
        
        %append to previous database, remove temporary file
        dbtiv=vertcat(oldinfotiv,dbtiv);
        dbtiv=sortrows(dbtiv,'Identifier');
        filenamedbtivs = sprintf('LEADS_TIVs_%s.csv', datestr(now,'mm-dd-yyyy_HH-MM-SS'));
        writetable(dbtiv,char(strcat(path_wscore_meta,filenamedbtivs)),'WriteRowNames',false);
        delete(char(fullfile(path_wscore_meta,'tempTIV.csv')));
        
    end % end if condition there are new MRIs for which we need TIVs


elseif checkinfo==1
    
    fprintf(2,'**There are no MRI scans to update\n\n');
    dbtiv=oldinfotiv;
    
end % end if condition there are no differences between old and new seg8 file lists

fprintf(" --------- Wmapping --------- \n");
% NOTE : Add the FDG for Whole Cerebellum 
listimgs=[dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/MRI*/s8iso*LDS*nu.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FBB*/wLDS*suvr_cbl.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FTP*/wLDS*suvr_infcblg.nii');dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FDG*/wLDS*suvr_pons.nii')]; % dir('/mnt/coredata/Projects/LEADS/data_f7p1/processed/*/*/FDG*/wLDS*suvr_wcbl.nii')];
listimgs = strcat({listimgs.folder}','/',{listimgs.name}');

% testpurposes only
listimgs = {'/mnt/coredata/Projects/LEADS/data_f7p1/processed/LDS0990416/Timepoint2/FDG_2022-12-20/wLDS0990416_FDG_2022-12-20_suvr_pons.nii'};

for v=1:size(listimgs,1)
    
    tempid=listimgs{v,1}(50:59);
    temptp=listimgs{v,1}(61:70);
    tempmod=listimgs{v,1}(72:74);
  
    if isequal(tempmod,'MRI')
        tempdate=datetime(listimgs{v,1}(79:88));
    else
        tempdate=datetime(listimgs{v,1}(76:85));
    end % end if condition what modality to grab date

    [p,f,e]=spm_fileparts(listimgs{v,1});

    fprintf(1,"Checking if the images already exist\n   ")

    if exist(strcat(p,'/Wmap-v2_',f,e),'file')==0
        
        fprintf(1,"Wmap-v2 does not exist for : %s\n", strcat(f));

        fprintf(1,"Looking up DOB/Age \n   ")
        qcdob=datetime(dbagef(ismember(dbagef.ID, tempid),:).dob, 'Format','yyyy-MM-dd');
        if size(qcdob,1)==1 && qcdob<datetime('1978-01-01') && qcdob>datetime('1953-01-01')

            tempimg=listimgs{v,1};
            tempage=calyears(between(qcdob,tempdate,'years'));
            tempsex=dbsex(ismember(dbsex.ID,tempid),:).sex;
            tempwmapname=strcat(p,'/Wmap-v2_',f,e);
            fprintf(1,"Generating Wmaps: %s\n", tempwmapname);

            if isequal(tempmod,'MRI')
            tempidentifier=cellstr(regexp(tempimg, 'LDS\d{7}_\w{3}_\w{2}_\d{4}-\d{2}-\d{2}','match','once'));
            temptiv=dbtiv(ismember(dbtiv.Identifier,tempidentifier),:).TIV;
            tempvols=vertcat(cellstr(tempimg),itcpt_mri,betaage_mri,betativ_mri,sdmap_mri);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')+(i4*',num2str(temptiv),')))./i5');
            elseif isequal(tempmod,'FBB')
            tempvols=vertcat(cellstr(tempimg),itcpt_fbb,betaage_fbb,sdmap_fbb);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')))./i4'); 
            elseif isequal(tempmod,'FTP')
            tempvols=vertcat(cellstr(tempimg),itcpt_ftp,betaage_ftp,sdmap_ftp);
            tempexpr=strcat('(i1-(i2+(i3*',num2str(tempage),')))./i4'); 
            
            elseif isequal(tempmod,'FDG')
                fprintf(1,"Depending on the ROI Selecting the calibration files \n");
                if contains(tempimg, 'pons')
                    tempvols = vertcat(cellstr(tempimg), itcpt_fdg_pons, betaage_fdg_pons, betasex_fdg_pons, sdmap_fdg_pons);
                elseif contains(tempimg, 'wcbl')
                    tempvols = vertcat(cellstr(tempimg), itcpt_fdg_whlcerebellum, betaage_fdg_whlcerebellum, betasex_fdg_whlcerebellum, sdmap_fdg_whlcerebellum);
                end    
                tempexpr = strcat('(i1-(i2+(i3*', num2str(tempage), ')+(i4*', num2str(tempsex), ')))./i5');
            end % end if condition with which modality we are working

            tempvols2=vertcat(cellstr(tempwmapname),cellstr(emask));
            tempwmapname_gm=char(strcat(p,'/Wmap-v2_',f,'_GM',e));
            disp(tempwmapname_gm)

            fprintf("Computing Wmap\n")
            clear matlabbatch
            spm('defaults','PET'); 
            matlabbatch{1}.spm.util.imcalc.input = tempvols;
            matlabbatch{1}.spm.util.imcalc.output = tempwmapname;
            matlabbatch{1}.spm.util.imcalc.outdir = '';
            matlabbatch{1}.spm.util.imcalc.expression =tempexpr;
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
            matlabbatch{2}.spm.util.imcalc.input = tempvols2;
            matlabbatch{2}.spm.util.imcalc.output = tempwmapname_gm;
            matlabbatch{2}.spm.util.imcalc.outdir = '';
            matlabbatch{2}.spm.util.imcalc.expression = 'i1.*i2';
            matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{2}.spm.util.imcalc.options.mask = 0;
            matlabbatch{2}.spm.util.imcalc.options.interp = 1;
            matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run',matlabbatch); clear matlabbatch;

            % Let's create the symbolic links to the Wmaps
            fprintf("Creating symbolic links\n")
            symlink1=strcat('ln -s',{' '},tempwmapname,{' '},strcat('/mnt/coredata/Projects/LEADS/data_f7p1/mni/',temptp,'/wmap',lower(tempmod),'/Wmap-v2_',f,e)); system(char(symlink1));
            
            % Generatig Multiaxial Images
            fprintf("Generating Multiaxial Images\n")
            run LEADS_Wmap_MultiAxial_service.m

            % Generating 3D Rendered Images
            fprintf("Generating 3D Rendered Images\n")
            run LEADS_Wmap_3DRend_service.m

            fprintf("Wmapping Services Completed!\n")

        elseif size(qcdob,1)==1 && (qcdob>datetime('1978-01-01') || qcdob<datetime('1953-01-01'))
        
            fprintf(2,'Warning! For %s the dob was %s, it is not in the 40-65 LEADS range.\nI am skipping the subject. Press Enter to acknowledge and continue:\n',tempid,char(qcdob));
            pause

        end
        
    else
        fprintf(1, "Exists: %s\n", strcat(p, '/Wmap-v2_', f, e));
        
    end % end if condition to check if the Wmap already exists 

    clear tempid temptp tempmod tempdate p f e qcdob tempimg tempage tempwmapname tempidentifier temptiv tempvols tempexpr tempvols2 tempwmapname_gm symlink1
    
end

%cd (path_processed);
%%% Let's copy the all the available JPEG files to the shared petcore
%system('cp -n   **/**/FTP*/Wmap*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/ftp/wmap/.');
%system('cp -n   **/**/FBB*/Wmap*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/fbb/wmap/.');
%system('cp -n   **/**/FDG*/Wmap*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/fdg/wmap/.');
%system('cp -n   **/**/MRI*/Wmap*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/mri/wmap/.');
%
%system('cp -n   **/**/FTP*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/ftp/wmap_3d/.');
%system('cp -n   **/**/FBB*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/fbb/wmap_3d/.');
%system('cp -n   **/**/FDG*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/fdg/wmap_3d/.');
%system('cp -n   **/**/MRI*/3DRend_W*.jpg   /shared/petcore/Projects/LEADS_training/data_f7p1/summary/mri/wmap_3d/.');
%
%disp('Finished, enjoy your Wmaps!');
%clear;