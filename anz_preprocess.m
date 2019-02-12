
function anz_preprocess(subjectdir, atlasdir, species, TR, trimf, ...
    stc, sm, fil, gsr, wmcsfr, slice, PB)

%% ________________________________________________________________________
%                                                                          
%                                                                          
%                      Anzar's Preprocessing Pipeline
%                                                                          
% _________________________________________________________________________
%                                                                          
%                                                                          
%                         Authored by: Anzar Abbas
%                          anzar.abbas@emory.edu
%                         Last updated 11/01/2016
%                                                                          
% This function can preprocess functional MRI data from an individual
% subject with an anatomical scan and one or more functional scans. The
% preprocessing can be applied to data from rodents, monkeys, and humans,
% as long as a standard space MRI atlas exists for that species. 
%
% It is important to note that this function relies heavily on FSL
% (https://fsl.fmrib.ox.ac.uk), which must be downloaded and installed.
% Additionally, this function also requires the use of a NIfTI toolbox
% by Jimmy Shen, which has been provided with this function.
%
% A document accompanying this function, named anz_preprocess_manual.pdf,
% should specify the details of the preprocessing pipeline being applied
% and how to address the  paramters that need to be inputted into the
% function. A shorter description of the parameters has been given here as
% well.
%                                                                        
%                                                                          
% _____________________________ subjectdir ________________________________
%                                                                          
% A string which specifies the entire path to the folder in which all of
% the subject's data is stored. This folder should have two types of files
%                                                                          
%   1 - The anatomical image, saved as a NIfTI file, and named 't1'.
%                                                                          
%   2 - The functional image(s), saved as a NIfTI file. One  subject may
%       have more than one functional scan. The name of the functional
%       scan must begin with 'f' and be  followed by two digits specifying
%       the scan number. For example, a subject with one functional scan
%       will have only one file named 'f01'. A subject with two or more
%       functional scans will have multiple functional scans named f01,
%       f02, ..., f10, and so on.
%                                                                          
% However, in case of data from rodents, the user will also have to have
% pre-made  anatomical and functional brain masks saved in the subject
% folder:
%                                                                          
%   3 - An anatomical mask that is a binary 3D image in which everything
%       that is brain is specified with the value of 1 and everything that
%       is not brain is specified with a value of 0.  This file must be
%       named 't1_mask'.
%                                                                          
%   4 - A functional mask for each scan that is a binary 3D image  in which
%       everything that is brain is specified with the value of 1 and
%       everything that is not brain is specified with a value of 0. This
%       file must  be named 'f??_mask',  with the question  marks being the
%       scan number.
%
% It is often helpful to have the subjectdir previously saved as a string
% variable in MATLAB.
%                                                                          
% _______________________________ atlasdir ________________________________
%                                                                          
% This is a  string which specifies the entire path to the folder in which
% the standard space atlas data is stored.  This folder should have the
% following files:
%                                                                          
%   1 - The atlas T1 brain image. This is an anatomical image of the
%       brain in which only neural tissue is present. This file must be
%       named 'atlas_t1.nii.gz'
%                                                                          
%   2 - The atlas brain mask image. This is a binary image in which
%       everything that is brain is specified with the value of 1 and
%       everything that is not brain is specified with a value of 0. This
%       file must be named 'atlas_t1_brain_mask.nii.gz'
%                                                                          
%   3 - In the case of rodents only, a segmented image of the brain. This
%       is a standard image distributed with atlases in which all areas
%       that  are CSF are specified with a value of 1, all areas that are
%       gray matter are specified with a value of 2, and all areas that are
%       white matter are specified with a value of 3. This file must be
%       named 'atlas_t1_brain_seg.nii.gz'
%
% Example atlases are provided with this function.
%
% It is often helpful to have the atlasdir previously saved as a string
% variable in MATLAB.
%                                                                          
% ________________________________ species ________________________________
%                                                                          
% A 1/2 input specifiying if you are working with  rodents (1) or primates
% (2). Primates can mean both monkeys and humans. Depending on the answer
% to this question, the preprocessing pipeline will carry out the
% appropriate steps.
%                                                                          
% __________________________________ TR ___________________________________
%                                                                          
% The TR of the functional scan(s) in seconds.
%                                                                          
% _________________________________ trim __________________________________
%                                                                          
% How many  timepoints the user wants to cut out from the beginning of the
% functional scan. The user can also enter 0 here to avoid trimming the
% functional scan.
%                                                                          
% __________________________________ stc __________________________________
%                                                                          
% Slice time correction - If the  user does not want to conduct slice time
% correction, this value should be inputted as 0. If the user does want to
% conduct slice time correction, then this will be a 1x2 matrix. The first
% value in the matrix is about the slice acquisision (SA). If the SA was
% bottom up, then enter 1. If the SA was top down, then enter 2. The second
% value in the matrix is about the slice intervals.  If  the  slices  were
% collected in series, then enter 1.  If  they  were  collected  in  an
% interleaved fashion, then enter 2 here.
%                                                                          
% ___________________________________ sm __________________________________
%                                                                          
% Smoothing - The  amount the user would like to spatially filter or smooth
% the functional data in mm. For reference, this is typically 6mm for
% humans, 2mm for macaques, and 0.5mm for rats.
%                                                                          
% ___________________________________ gsr _________________________________
%                                                                          
% Global signal regression - 1/0 input describing whether or not the user
% would like to conduct global signal regression on the functional data.
%                                                            
% _________________________________ wmcsfr ________________________________
%                                                                          
% White matter and CSF signal regression - 1/0  input  describing whether
% or not the user would like to conduct white matter and CSF signal
% regression on the functional data. This value is always 0 for rodents as
% this function does not conduct WM/CSF signal regression in rodents.
%                                                                          
% ___________________________________ fil _________________________________
%                                                                          
% Temporal filter - All functional scans will undergo tempral filtering.
% fil is a 1x2 matrix in which the first value is the high pass filter and
% the second value is the low pass filter. Both values are in Hz. 
%                                                                          
% _________________________________ slice _________________________________
%                                                                          
% Reduce by one spatial dimension - 1/0 input describing whether or not the
% user would like to slice all 3D brain volumes (anatomical and functional)
% into 2D matrices with all the slices for each timepoint laid out in one
% image. This will convert 3D anatomical images to 2D anatomical images and
% 4D functional images to 3D functional images. This is helpful for
% visualization of all slices simultaneously.
%                                                                          
% ___________________________________ PB __________________________________
%                                                                          
% PushBullet - Optional input - Pushbullet is a tool that I use to  keep me
% updated on processes I know will take a long time. This preprocessing
% pipeline is one of them. Entering a pushbullet API key specific to the
% user here can keep them updated on the progress of the preprocessing
% pipeline. What's  even more useful is that they will also be notified if
% the code errs, so they can tend to it sooner rather than later. To get
% your own API key, follow the instructions on the Pushbullet website.
%                                                                          
% _________________________________________________________________________




%% ______________________________________________________________________ %
%                                                                         %
%                                   Setup                                 %
% _______________________________________________________________________ %

fprintf('\nBeginning preprocessing pipeline\n\n');
if ismac == 1
    fsldir = '/Applications/fsl/bin/';
    setenv('FSLDIR','/Applications/fsl');
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
else
    fsldir = '/usr/local/fsl/bin/';
    setenv('FSLDIR','/usr/local/fsl');
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
end
% Setting FSL environment and output filetype to '.nii.gz' Saving the
% directories in which all the FSL functions are stored. This part of the
% script is making the assumption that the FSL directory on your computer
% is in the default location that FSL uses when it is being installed. If
% that is not the case, you should change the path of the directory in the
% code above.

fprintf('Pushbullet key ')
if nargin > 12
    p = Pushbullet(PB);
    fprintf('assigned\n')
else
    fprintf('not assigned\n\n');
end
% Entering the API key for the user if they care to use pushbullet during
% the preprocessing pipeline (recommended by author but not at all
% necessary).

k = strfind(subjectdir,'subject');
sub = [subjectdir(k:end),' - ']; 
clear k
% This is simply helping with the text that will be displayed as the
% function is running. 'sub' is a string with the name of the subject that
% will be printed with every update

warning('off','all');
% There are a few warnings that MATLAB outputs depending on the type of
% computer this function that is being run on. As far as my understanding
% goes, they are irrelevant, hence the function is temporarily turning off
% warnings so that the updates being printed as the code runs look clean.

% 
% _______________________________________________________________________ %

cd (subjectdir)

switch species 
    case 1
        
        
        
%% ______________________________________________________________________ %
%                                                                         %
%                                 Rodents                                 %
% _______________________________________________________________________ %
        
%         fprintf('Working on rodent data\n\n')
%         
%         % Anatomical Image Reorientation ________________________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Anatomical image reorientation ... '])
%         % Every step in this function will print an update in the command
%         % window. That update will include the time and the name of the
%         % subject and scan being preprocessed. The above three lines of
%         % code will repeat with every step.
%         cmd = [fsldir,'fslreorient2std t1 t1_reorient'];
%         system(cmd);
%         fprintf('Done\n')
%         % Reorienting the anatomical image to standard just in case it
%         % isn't already.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Anatomical Mask Reorientation _________________________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Anatomical mask reorientation ... '])
%         cmd = [fsldir,'fslreorient2std t1_mask t1_mask'];
%         system(cmd);
%         fprintf('Done\n')
%         % Reorienting the user-created anatomical image mask to standard
%         % just in case it isn't already.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Anatomical Brain Extraction ___________________________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Anatomical brain extraction  ... '])
%         cmd = [fsldir,'fslmaths t1_reorient -mul t1_mask', ...
%             ' t1_reorient_brain'];
%         system(cmd);
%         fprintf('Done\n')
%         % Multiplying the anatomical image by the brain mask to get a brain
%         % extracted image.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Anatomical Brain Registration to Atlas ________________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Anatomical brain registration to atlas ... '])
%         cmd = [fsldir,'flirt -in t1_reorient_brain -ref ', ...
%             atlasdir,'/atlas_t1_brain -out t1_reorient_brain_regatlas', ...
%             ' -omat t1_to_atlas.mat'];
%         system(cmd);
%         system('rm t1_to_atlas.mat');
%         fprintf('Done\n')
%         % Using FSL's FLIRT tool to register the anatomcial images to the
%         % rat atlas brain image.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Anatomical Brain Segmentation _________________________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Anatomical brain segmentation ... '])
%         nii = load_untouch_nii([atlasdir,'/atlas_t1_brain_seg.nii.gz']);
%         atlas_seg = nii.img;
%         % Loading the atlas file that is all three tissues segmented.
%         csf = 1 .* (atlas_seg == 1); 
%         % Creating CSF mask
%         gm = 1 .* (atlas_seg == 2); 
%         % Creating gray matter mask
%         wm = 1 .* (atlas_seg == 3); 
%         % Creating white matter mask
%         nii.img = csf;
%         save_untouch_nii(nii,'t1_reorient_brain_regatlas_csf.nii.gz');
%         % Saving the CSF binary mask 3D image
%         nii.img = gm;
%         save_untouch_nii(nii,'t1_reorient_brain_regatlas_gm.nii.gz');
%         % Saving the GM binary mask 3D image
%         nii.img = wm;
%         save_untouch_nii(nii,'t1_reorient_brain_regatlas_wm.nii.gz');
%         % Saving the WM binary mask 3D image
%         clear nii atlas_seg 
%         fprintf('Done\n')
%         % Normally (in the case of primates), the function uses FSL's FAST
%         % tool to segment the individual's anatomical scan. However, FAST
%         % does not work reliably on rodent data. Since the data is
%         % registered to a standard space atlas, the function instead uses
%         % the tissue segmentation masks from the atlas.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Slicing anatomical images into 2D matrices ____________________ %
%         %
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,'Slicing anatomical images into 2D images ... '])
%         if slice == 1
%             nii = load_nii('t1_reorient_brain_regatlas.nii.gz');
%             anatomical_images_3d = cell(4,1);
%             % This will hold the t1 image as well as all the segmented
%             % tissue maps
%             anatomical_images_3d{1} = nii.img;
%             % The first index holds the t1 image
%             anatomical_images_3d{2} = csf;
%             anatomical_images_3d{3} = gm;
%             anatomical_images_3d{4} = wm;
%             % The other three indices hold the CSF, gray matter, and white
%             % matter tissue maps respectively
%             [~,~,z] = size(anatomical_images_3d{1});
%             % Getting dimensions of T1 image
%             num_boxes = z;
%             % Predefining the number of boxes our image is going to have
%             if rem(sqrt(z),1) >= eps
%                 % If the brain in the z dimension is not a perfect square
%                 while rem(sqrt(z),1) >= eps
%                     % Well then while that's the case
%                     num_boxes = num_boxes + 1;
%                     % Adding one to the number of boxes
%                 end
%             end
%             anatomical_image_slices = cell(4,1);
%             % This cell array will hold all the cell arrays that will hold
%             % all the slices of the brains for each of the anatomical
%             % images
%             anatomical_images_2d = cell(4,1);
%             % This is the cell array that will hold the 2D image for each
%             % of the anatomical images
%             for i = 1:4
%                 anatomical_image_slices{i} = cell(num_boxes,1);
%                 % This is the cell array that will hold the slices for one
%                 % anatomical image
%                 for j = 1:z
%                     temp = squeeze(anatomical_images_3d{i}(:,:,j));
%                     % For each slice in each anatomical image
%                     temp = rot90(temp,1);
%                     % Rotate to make anterior face upwards
%                     anatomical_image_slices{i}{j} = temp;
%                 end
%                 [x,y] = size(anatomical_image_slices{i}{1});
%                 % Getting x and y dimensions of each box
%                 empty_space = zeros(x,y);
%                 % Making empty space for the empty boxes
%                 for j = z+1:z+(num_boxes-z)
%                     % Going through the empty boxes
%                     anatomical_image_slices{i}{j} = empty_space;
%                     % And filling them with empty space
%                 end
%                 anatomical_images_2d{i} = zeros(x*sqrt(num_boxes), ...
%                     y*sqrt(num_boxes));
%                 % Predefining the 2D anatomical image
%                 count = 0;
%                 % Going to count through the slices
%                 for row = 1:sqrt(num_boxes)
%                     for column = 1:sqrt(num_boxes)
%                         count = count + 1;
%                         % Going through each row and column
%                         anatomical_images_2d{i}((row*x-x+1):(row*x), ...
%                             (column*y-y+1):(column*y)) = ...
%                             anatomical_image_slices{i}{count};
%                         % Organizing the brain slices in rows and columns
%                     end
%                 end
%             end
%             anat = anatomical_images_2d{1}; %#ok<*NASGU>
%             csf = anatomical_images_2d{2};
%             gm = anatomical_images_2d{3};
%             wm = anatomical_images_2d{4};
%             % Renaming all the sliced anatomical images that have been laid
%             % out as 2D images
%             save('t1_reorient_brain_regatlas_2D.mat','anat');
%             save('t1_reorient_brain_regatlas_csf_2D.mat','csf');
%             save('t1_reorient_brain_regatlas_gm_2D.mat','gm');
%             save('t1_reorient_brain_regatlas_wm_2D.mat','wm');
%             % Saving each of those matrices as .mat files
%             clear nii anatomical_images_3d z anatomical_images_slices
%             clear anatomical_images_2d i j temp x y empty_space count 
%             clear row column anat
%             fprintf('Done\n')
%         else
%             fprintf('Skipped\n')
%         end
%         % This step cuts the 3D anatomical images slice by slice and lays
%         % out all the slices in one plane. If the user chooses to do this,
%         % they can view all slices of the brain at the same time by loadinf
%         % the .mat files that have just been saved.
%         % _______________________________________________________________ %
%         
%         
%         
%         % Going through each functional scan ____________________________ %
%         %
%         f_scans = dir('f*');
%         f_scans = {f_scans.name}';
%         count = 0;
%         indices = 0;
%         for i = 1:length(f_scans)
%             if length(f_scans{i}) > 10
%                 count = count + 1;
%                 indices(count) = i; %#ok<*AGROW>
%             end
%         end
%         try
%             f_scans(indices) = [];
%         catch
%         end
%         for i = 1:length(f_scans)
%             f_scans{i} = f_scans{i}(1:3);
%         end
%         clear i
%         % Creating a cell array of all the raw functional scans in which
%         % just the name of the functional scans are save in a string.
%         for i = 1:length(f_scans)
%             filename = f_scans{i};
%             scan = [f_scans{i},' - '];
%             fprintf('\n')
%         % This is the file that we're currently working on at all times and
%         % it will update with each step that we conduct
%         %
%         % The code above is simply setting up a for loop for every
%         % functional scan that this particular subject has. If the code is
%         % erring at this part of the function, please make sure that all
%         % functional scans have been named in the way specified above and
%         % in the manual.
%         % _______________________________________________________________ %
%             
%             
%             
%             % Reorienting functional scan _______________________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Reorienting functional scan ... '])
%             cmd = [fsldir,'fslreorient2std ', ...
%                 filename,' ', filename, '_reorient'];
%             system(cmd);
%             filename = [filename,'_reorient'];
%             fprintf('Done\n')
%             % Reorienting the functional scan to standard, just in case it
%             % isn't already
%             % ___________________________________________________________ %
%             
%             
%             
%             % Trimming functional scan __________________________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Trimming functional scan ... '])
%             if trimf ~= 0
%                 cmd = [fsldir,'fslnvols ',filename];
%                 [~,timepoints] = system(cmd);
%                 % Finding out how many total timepoints are in the
%                 % functional scan as we will be needing that number in the
%                 % next step
%                 cmd = [fsldir,'fslroi ',filename,' ',filename,'_trim ', ...
%                     num2str(trimf),' ',num2str(str2double((timepoints)))];
%                 system(cmd);
%                 filename = [filename,'_trim'];
%                 fprintf('Done\n')
%             else
%                 fprintf('Skipped\\n')
%             end
%             clear trim timepoints
%             % Trimming the front of the functional scan, if specified by
%             % user - and updating the filename to include _trim in it if we
%             % did conduct the trimming. Reasoning for this trimming is
%             % included in the manual.
%             % ___________________________________________________________ %
%             
%             
%             
%             % Slice time correction _____________________________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Slice time correction ... '])
%             if stc ~= 0
%                 if stc(1) == 1
%                     if stc(2) == 1
%                         cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
%                             filename,'_stc -r ',num2str(TR)];
%                         system(cmd);
%                         % Running slice time correction with slices
%                         % acquired bottom up and in series
%                     else
%                         cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
%                             filename,'_stc -r ',num2str(TR),' --odd'];
%                         system(cmd);
%                         % Running slice time correction with slices
%                         % acquired bottom up and interleaved
%                     end
%                 else
%                     if stc(2) == 1
%                         cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
%                             filename,'_stc -r ',num2str(TR), ' --down'];
%                         system(cmd);
%                         % Running slice time correction with slices
%                         % acquired top down and in series
%                     else
%                         cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
%                             filename,'_stc -r ',num2str(TR), ...
%                             ' --down --odd'];
%                         system(cmd);
%                         % Running slice time correction with slices
%                         % acquired bottom down and interleaved
%                     end
%                 end
%                 filename = [filename,'_stc'];
%                 % Updating the filename only if we slice time corrected
%                 fprintf('Done\n')
%             else
%                 fprintf('Skipped\n')
%             end
%             % This whole if statement is conducting the slice time
%             % correction if it was specified by the user in the way that it
%             % was specified by the user. If there are further complications
%             % in slice time correction, the code above can be adjusted to
%             % accomodate the way the scans were collected.
%             % ___________________________________________________________ %
%             
%             
%             
%             % Motion correction of functional scans _____________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Motion correction ... '])
%             cmd = [fsldir,'mcflirt -in ',filename,' -out ', ...
%                 filename,'_mc -refvol 0 -plots'];
%             system(cmd);
%             filename = [filename,'_mc'];
%             % Conducting motion correction on the functional scans and
%             % saving all the motion paramters. Below, I will plot out the
%             % motion parameters as an output.
%             par = load([filename,'.par']);
%             fig = figure('visible','off','position',get(0,'screensize'));
%             % Creating a hidden figure that will eventually be saved as a
%             % JPEG file for the user's convenience.
%             subplot(2,1,1);
%             hold on;
%             plot(par(:,1));
%             plot(par(:,2));
%             plot(par(:,3));
%             title('Rotation in radians');
%             legend('X','Y','Z');
%             set(gca,'ylim',[-0.2 0.2],'fontsize',12);
%             ylabel('Rotations (rad)');
%             xlabel('Timepoints');
%             % Plotting the rotation paramters
%             subplot(2,1,2);
%             hold on
%             plot(par(:,4));
%             plot(par(:,5));
%             plot(par(:,6));
%             title('Transformation in mm');
%             legend('X','Y','Z');
%             set(gca,'ylim',[-3 3],'fontsize',12);
%             ylabel('Transformations (mm)');
%             xlabel('Timepoints');
%             % Plotting the transformation parameters
%             saveas(fig,[filename,'_visual.jpg']);
%             % Saving
%             fprintf('Done\n')
%             % This series of steps is using FSL's MCFLIRT tool to conduct
%             % motion correction and then plotting out and saving the motion
%             % parameters as a separate JPEG file that can be viewed later.
%             % ___________________________________________________________ %
%             
%             
%             
%             % Functional Mask Reorientation _____________________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Functional scan reorientation ... '])
%             cmd = [fsldir,'fslreorient2std ',scan(1:3),'_mask ', ...
%                 scan(1:3),'_mask'];
%             system(cmd);
%             % Reorienting the user-created functional mask to standard just
%             % in case it isn't already
%             fprintf('Done\n')
%             % ___________________________________________________________ %
%             
%             
%             
%             % Functional Brain Extraction _______________________________ %
%             % 
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Brain extraction from skull ... '])
%             cmd = [fsldir,'fslmaths ', filename, ...
%                 ' -mul ',scan(1:3),'_mask ',filename,'_brain'];
%             system(cmd);
%             filename = [filename,'_brain'];
%             % Multiplying the functional scan by the brain mask to get a
%             % brain extracted image
%             fprintf('Done\n')
%             % 
%             % ___________________________________________________________ %
%             
%             
%             
%             % Functional Brain Registration to Atlas ____________________ %
%             % 
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Registration to atlas space ... '])
%             cmd = [fsldir,'flirt -in ',filename,' -ref ', ...
%                 atlasdir,'/atlas_t1_brain -omat f_to_atlas.mat'];
%             system(cmd);
%             cmd = [fsldir,'flirt -in ',filename,' -ref ', ...
%                 atlasdir,'/atlas_t1_brain -applyxfm -init ', ...
%                 'f_to_atlas.mat -out ',filename,'_regatlas'];
%             system(cmd);
%             cmd = [fsldir,'fslmaths ',filename,'_regatlas -mas ', ...
%                 atlasdir,'/atlas_t1_brain_mask ',filename,'_regatlas'];
%             system(cmd);
%             system('rm f_to_atlas.mat');
%             filename = [filename,'_regatlas'];
%             % Using FSL's FLIRT to register ther functional scan to the
%             % atlas brain image
%             fprintf('Done\n')
%             % ___________________________________________________________ %
%             
%             
%             
%             % Smoothing of Functional Scan ______________________________ %
%             % 
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Spatial smoothing ... '])
%             cmd = [fsldir,'fslmaths ', filename, ' -s ', num2str(sm), ...
%                 ' ',filename,'_sm'];
%             system(cmd);
%             filename = [filename,'_sm'];
%             % Using FSL's fslmaths tool to smooth the functional scans by
%             % however much the user specified
%             fprintf('Done\n')
%             % ___________________________________________________________ %
%             
%             
%             
%             % Temporal Filtering of functional data _____________________ %
%             % 
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Temporal filtering ... '])
%             % Starting this step, we are no longer using FSL for
%             % preprocessing. Rather, we're loading the NIfTI files into
%             % MATLAB and performing some steps on them directly. This is
%             % why the NIfTI toolbox provided with this function is
%             % critical.
%             nii = load_untouch_nii([filename,'.nii.gz']);
%             f = nii.img;
%             [X,Y,Z,T] = size(f);
%             % Getting the dimensions of the functional scan
%             f = reshape(f,[X*Y*Z,T]); 
%             % Reshaping the scan into a 2D matrix of voxels over time
%             f = f';
%             tc = f;
%             SF = 1/TR';
%             % Calculating the sampling frequency using the TR
%             x = size(tc,1);
%             tc_mirror = tc(end:-1:1,:);
%             tc_extended = [tc_mirror; tc; tc_mirror];
%             tc_fft = fft(tc_extended);
%             high_cutoff = round(size(tc_extended,1)*(fil(1)/SF));
%             tc_fft(1:high_cutoff,:) = 0;
%             low_cutoff = round(size(tc_extended,1)*(fil(2)/SF));
%             tc_fft(low_cutoff:end,:) = 0;
%             tc_filtered = ifft(tc_fft);
%             temp = 2*real(tc_filtered);
%             tc_filtered = temp(x+1:2*x,:);
%             f = tc_filtered';
%             % The above steps are applying a FFT filter on the functional
%             % data
%             f = reshape(f,[X,Y,Z,T]);
%             nii.img = f;
%             save_untouch_nii(nii,[filename,'_fil.nii.gz']);
%             cmd = [fsldir,'fslmaths ',filename,'_fil -mas ', atlasdir, ...
%                 '/atlas_t1_brain_mask ',filename,'_fil'];
%             system(cmd);
%             % Just to cut down on file size, I'm multiplying the functional
%             % scan by the brain mask again
%             filename = [filename,'_fil'];
%             fprintf('Done\n')
%             clear f1 tc_filtered temp tc_fft low_cutoff high_cutoff
%             clear tc_extended tc_mirror x
%             % ___________________________________________________________ %
%  
%             
%             
%             % Global signal regression __________________________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Global signal regression ... '])
%             if gsr == 1
%                 nii = load_untouch_nii([filename,'.nii.gz']);
%                 f = nii.img;
%                 % Loading the functional scan into Matlab
%                 f_gsr = f;
%                 % Predefining the global signal regressed functional scan
%                 f_mean = zeros(1,size(f,4));
%                 % Predefining what will me the mean of the functional scan
%                 for t = 1:size(f,4)
%                     temp = f(:,:,:,t);
%                     temp_sum = sum(temp(:));
%                     temp_length = length(find(temp(:)));
%                     temp_mean = temp_sum/temp_length;
%                     f_mean(t) = temp_mean;
%                 end
%                 % Taking the mean of the functional timeseries. This is the
%                 % functional timecourse.
%                 f_mean_zsc = zscore(f_mean);
%                 % Z-scoring the mean
%                 for x = 1:size(f,1)
%                     for y = 1:size(f,2)
%                         for z = 1:size(f,3)
%                             v = double(squeeze(f(x,y,z,:)));
%                             beta = (f_mean_zsc * f_mean_zsc') \ ...
%                                 (f_mean_zsc * v);
%                             f_gsr(x,y,z,:) = v' - f_mean_zsc * beta;
%                         end
%                     end
%                 end
%                 % Regressing the global signal from the functional
%                 % timeseries
%                 f = f_gsr;
%                 filename = [filename,'_gsr'];
%                 nii.img = f;
%                 save_untouch_nii(nii,[filename,'.nii.gz']);
%                 % Saving the global signal regressed functional scan
%                 fprintf('Done\n')
%             else
%                 fprintf('Skipped\n')
%             end
%             % Conducting global signal regression using GLM regression on
%             % the functional data if the user specified.
%             % ___________________________________________________________ %
%             
%             
%             
%             % Z-scoring all voxels in functional scan ___________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Z-scoring ... '])
%             if gsr ~= 1
%                 nii = load_untouch_nii([filename,'.nii.gz']);
%                 f = nii.img;
%                 % Loading the functional scan in case it wasn't already
%                 % loaded in the previous step
%             end
%             mask = load_nii([atlasdir,'/atlas_t1_brain_mask.nii.gz']);
%             mask = mask.img;
%             % Loading the T1 mask
%             [x,y,z,~] = size(f);
%             % Assigning variables to the dimensions of the scan
%             for j = 1:x
%                 for k = 1:y
%                     for l = 1:z
%                         if mask(j,k,l) == 1
%                             voxel = f(j,k,l,:);
%                             % Working on one brain voxel at a time
%                             voxel_zsc = zscore(voxel);
%                             % Z-scoring the voxel signal
%                             f(j,k,l,:) = voxel_zsc;
%                             % Adding the zscored value to the zscored
%                             % functional scan matrix
%                         end
%                     end
%                 end
%             end
%             nii.img = f;
%             save_untouch_nii(nii,[filename,'_zsc.nii.gz']);
%             filename = [filename,'_zsc'];
%             clear nii x y z j k l
%             fprintf('Done\n')
%             % Saving the newly z-scored functional scan
%             % ___________________________________________________________ %
%             
%             
%             
%             % Slicing functional images into 3D matrices ________________ %
%             %
%             T = datetime('now');
%             time = whatsthetime(T);
%             fprintf([time,sub,scan,'Slicing into 3D matrices ... '])
%             if slice == 1
%                 [~,~,z,~] = size(f);
%                 % Getting dimensions of functional scan
%                 num_boxes = z;
%                 % Predefining the number of boxes the output image is going
%                 % to have.
%                 if rem(sqrt(z),1) >= eps
%                     % If the brain in the z dimension is not a perfect
%                     % square
%                     while rem(sqrt(num_boxes),1) >= eps
%                         % Well then while that's the case
%                         num_boxes = num_boxes + 1;
%                         % Adding to the number of boxes
%                     end
%                 end
%                 slices = cell(num_boxes,1);
%                 % Predefining cell array that will hold timeseries for all
%                 % slices
%                 for j = 1:z
%                     % Going through all slices
%                     temp = squeeze(f(:,:,j,:));
%                     % This is the timeseries for a slice
%                     temp = rot90(temp,1);
%                     % Rotating the image to make anterior face upwards
%                     slices{j} = temp;
%                     % Adding timeseries for each slice in each cell
%                 end
%                 [x,y,t] = size(slices{1});
%                 % Updating these values because of the rotation
%                 empty_space = zeros(x,y,t);
%                 % Empty space to go in empty boxes
%                 for j = z+1:z+(num_boxes-z)
%                     % Going through the rest of the cells in slices
%                     slices{i} = empty_space;
%                     % Adding empty space into cell
%                 end
%                 f_sliced = zeros(x*sqrt(num_boxes),y*sqrt(num_boxes),t);
%                 % Predefining the output image shape
%                 count = 0;
%                 % Will be counting through the boxes
%                 for row = 1:sqrt(num_boxes)
%                     % For every row of brain slices
%                     for column = 1:sqrt(num_boxes)
%                         % For every column of brain slices
%                         count = count + 1;
%                         % Updating the count
%                         f_sliced((row*x-x+1):(row*x),(column*y-y+1): ...
%                             (column*y),:) = slices{count};
%                         % Organizing brain slices in rows and columns
%                     end
%                 end
%                 func = f_sliced;
%                 filename = [filename,'_3D'];
%                 save([filename,'.mat'],'func');
%                 fprintf('Done\n')
%             else
%                 fprintf('Skipped\n')
%             end
%             % The abve steps are slicing the 4D functional scans into 3D
%             % functional scans by removing one of the spatial dimensions.
%             % ___________________________________________________________ %
%             
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,scan,'Finished preprocessing\n'])
%         
%         if nargin > 12
%             string = [sub,scan,'Functional scan preprocessed'];
%             p.pushNote([],string,'');
%             % Keeping the user updated with progress for each scan
%         end
%         
%         end
%         % Ending everything that we were doing to each functional scan
%         
% % _______________________________________________________________________ %



    case 2 % Primates
        

        
%% ______________________________________________________________________ %
%                                                                         %
%                                 Primates                                %
% _______________________________________________________________________ %
        


        fprintf('Working on primate/human data\n\n')
        
        
        
        % Anatomical Image Reorientation ________________________________ %
        %
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,'Anatomical brain reorientation ... '])
        % Every step in this function will print an update in the command
        % window. That update will include the time and the name of the
        % subject and scan being preprocessed. The above three lines of
        % code will repeat with every step.
        cmd = [fsldir,'fslreorient2std t1 t1_reorient'];
        system(cmd);
        fprintf('Done\n')
        % Reorienting the anatomical image to standard just in case it
        % isn't already.
        % _______________________________________________________________ %
        
        
        
        % Anatomical Brain Registration to Atlas ________________________ %
        %
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,'Anatomical brain registration to atlas ... '])
        cmd = [fsldir,'flirt -in t1 -ref ', atlasdir, ...
            '/atlas_t1 -out t1_reorient_regatlas', ...
            ' -omat t1_to_atlas.mat'];
        system(cmd);
        system('rm t1_to_atlas.mat');
        fprintf('Done\n')
        % Using FSL's flirt to register the anatomcial images to the atlas
        % anatomical image.
        % _______________________________________________________________ %
        
        
        
        % Anatomical Brain Extraction ___________________________________ %
        % 
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,'Anatomical brain extraction ...'])
        cmd = [fsldir,'fslmaths t1_reorient_regatlas -mas ', atlasdir, ...
            '/atlas_t1_brain_mask t1_reorient_regatlas_brain'];
        system(cmd);
        fprintf('Done\n')
        % Using FSL's fslmaths tool to separate out the brain from the
        % anatomical image
        % _______________________________________________________________ %
        
        
        % Anatomical Brain Segmentation _________________________________ %
        %
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,'Anatomical brain segmentation ... '])
        cmd = [fsldir,'fast -g t1_reorient_regatlas_brain'];
        system(cmd);
        system('rm t1_reorient_regatlas_brain_seg.nii.gz');
        system(['mv t1_reorient_regatlas_brain_seg_0.nii.gz ', ...
            't1_reorient_regatlas_brain_csf.nii.gz']);
        system(['mv t1_reorient_regatlas_brain_seg_1.nii.gz ', ...
            't1_reorient_regatlas_brain_gm.nii.gz']);
        system(['mv t1_reorient_regatlas_brain_seg_2.nii.gz ', ...
            't1_reorient_regatlas_brain_wm.nii.gz']);
        system('rm *pve*'); system('rm *mixeltype*');
        % Using FSL's FAST tool to separate out the three main tissue types
        % in the brain and then renaming them so that their filenames are
        % more indicative of what they actually are
        fprintf('Done\n')
        % _______________________________________________________________ %
        
        
        
        % Slicing anatomical images into 2D matrices ____________________ %
        %
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,'Anatomical brain slicing ... '])
        if slice == 1
            anatomical_images_3d = cell(4,1);
            % This will hold the t1 image as well as all the segmented
            % tissue maps
            nii = load_nii('t1_reorient_regatlas_brain.nii.gz');
            anatomical_images_3d{1} = nii.img;
            % The first index holds the t1 image
            nii = load_nii('t1_reorient_regatlas_brain_csf.nii.gz');
            csf = nii.img;
            anatomical_images_3d{2} = csf;
            % Loading and adding the CSF mask
            nii = load_nii('t1_reorient_regatlas_brain_gm.nii.gz');
            gm = nii.img;
            anatomical_images_3d{3} = gm;
            % Loading and adding the GM mask
            nii = load_nii('t1_reorient_regatlas_brain_wm.nii.gz');
            wm = nii.img;
            anatomical_images_3d{4} = wm;
            % Loading and adding the WM mask/keilholz-lab/Anzar/Tools
            [~,~,z] = size(anatomical_images_3d{1});
            % Getting dimensions of T1 image
            num_boxes = z;
            % Predefining the number of boxes our image is going to have
            if rem(sqrt(num_boxes),1) >= eps
                % If the brain in the z dimension is not a perfect square
                while rem(sqrt(num_boxes),1) >= eps
                    % Well then while that's the case
                    num_boxes = num_boxes + 1;
                    % Adding one to the number of boxes
                end
            end
            anatomical_image_slices = cell(4,1);
            % This cell array will hold all the cell arrays that will hold
            % all the slices of the brains for each of the anatomical
            % images
            anatomical_images_2d = cell(4,1);
            % This is the cell array that will hold the 2D image for each
            % of the anatomical images
            for i = 1:4
                anatomical_image_slices{i} = cell(num_boxes,1);
                % This is the cell array that will hold the slices for one
                % anatomical image
                for j = 1:z
                    temp = squeeze(anatomical_images_3d{i}(:,:,j));
                    % For each slice in each anatomical image
                    temp = rot90(temp,1);
                    % Rotate to make anterior face upwards
                    anatomical_image_slices{i}{j} = temp;
                end
                [x,y] = size(anatomical_image_slices{i}{1});
                % Getting x and y dimensions of each box
                empty_space = zeros(x,y);
                % Making empty space for the empty boxes
                for j = z+1:z+(num_boxes-z)
                    % Going through the empty boxes
                    anatomical_image_slices{i}{j} = empty_space;
                    % And filling them with empty space
                end
                anatomical_images_2d{i} = zeros(x*sqrt(num_boxes), ...
                    y*sqrt(num_boxes));
                % Predefining the 2D anatomical image
                count = 0;
                % Going to count through the slices
                for row = 1:sqrt(num_boxes)
                    for column = 1:sqrt(num_boxes)
                        count = count + 1;
                        % Going through each row and column
                        anatomical_images_2d{i}((row*x-x+1):(row*x), ...
                            (column*y-y+1):(column*y)) = ...
                            anatomical_image_slices{i}{count};
                        % Organizing the brain slices in rows and columns
                    end
                end
            end
            anat = anatomical_images_2d{1};
            csf = anatomical_images_2d{2};
            gm = anatomical_images_2d{3};
            wm = anatomical_images_2d{4};
            % Renaming all the sliced anatomical images that have been laid
            % out as 2D images
            save('t1_reorient_brain_regatlas_2D.mat','anat');
            save('t1_reorient_brain_regatlas_csf_2D.mat','csf');
            save('t1_reorient_brain_regatlas_gm_2D.mat','gm');
            save('t1_reorient_brain_regatlas_wm_2D.mat','wm');
            % Saving each of those matrices as .mat files
            clear nii anatomical_images_3d z anatomical_images_slices
            clear anatomical_images_2d i j temp x y empty_space count 
            clear row column anat
            fprintf('Done\n')
        else
            fprintf('Skipped\n')
        end
       % This step cuts the 3D anatomical images slice by slice and lays
        % out all the slices in one plane. If the user chooses to do this,
        % they can view all slices of the brain at the same time by loadinf
        % the .mat files that have just been saved.
        % _______________________________________________________________ %
        
        
        
        % Going through each functional scan ____________________________ %
        %
        f_scans = dir('f*');
        f_scans = {f_scans.name}';
        count = 0;
        indices = 0;
        for i = 1:length(f_scans)
            if length(f_scans{i}) > 10
                count = count + 1;
                indices(count) = i;
            end
        end
        try
            f_scans(indices) = [];
        catch
        end
        for i = 1:length(f_scans)
            f_scans{i} = f_scans{i}(1:3);
        end
        clear i length_filenames 
        % Creating a cell array of all the raw functional scans in which
        % just the name of the functional scans are save in a string.
        for i = 1:length(f_scans)
            filename = f_scans{i};
            scan = [f_scans{i},' - '];
            fprintf('\n')
        % This is the file that we're currently working on at all times and
        % it will update with each step that we conduct
        %
        % The code above is simply setting up a for loop for every
        % functional scan that this particular subject has. If the code is
        % erring at this part of the function, please make sure that all
        % functional scans have been named in the way specified above and
        % in the manual.        
        % _______________________________________________________________ %
            
            
            
            % Reorienting functional scan _______________________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Reorienting functional scan ... '])
            cmd = [fsldir,'fslreorient2std ', ...
                filename,' ', filename, '_reorient'];
            system(cmd);
            filename = [filename,'_reorient'];
            fprintf('Done\n')
            % Reorienting the functional scan to standard, just in case it
            % isn't already
            % ___________________________________________________________ %
            
            
            
            % Trimming functional scan __________________________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Trimming functional scan ... '])
            if trimf ~= 0
                cmd = [fsldir,'fslnvols ',filename];
                [~,timepoints] = system(cmd);
                % Finding out how many total timepoints are in the
                % functional scan as we will be needing that number in the
                % next step
                cmd = [fsldir,'fslroi ',filename,' ',filename,'_trim ', ...
                    num2str(trimf),' ',num2str(str2double((timepoints)))];
                system(cmd);
                filename = [filename,'_trim'];
                fprintf('Done\n')
            else
                fprintf('Skipped\n')
            end
            clear timepoints
            % Trimming the front of the functional scan, if specified by
            % user - and updating the filename to include _trim in it if we
            % did conduct the trimming. Reasoning for this trimming is
            % included in the manual.
            % ___________________________________________________________ %
            
            
            
            % Slice time correction _____________________________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Slice time correction ... '])
            if stc ~= 0
                if stc(1) == 1
                    if stc(2) == 1
                        cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
                            filename,'_stc -r ',num2str(TR)];
                        system(cmd);
                        % Running slice time correction with slices
                        % acquired bottom up and in series
                    else
                        cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
                            filename,'_stc -r ',num2str(TR),' --odd'];
                        system(cmd);
                        % Running slice time correction with slices
                        % acquired bottom up and interleaved
                    end
                else
                    if stc(2) == 1
                        cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
                            filename,'_stc -r ',num2str(TR), ' --down'];
                        system(cmd);
                        % Running slice time correction with slices
                        % acquired top down and in series
                    else
                        cmd = [fsldir,'slicetimer -i ',filename,' -o ', ...
                            filename,'_stc -r ',num2str(TR), ...
                            ' --down --odd'];
                        system(cmd);
                        % Running slice time correction with slices
                        % acquired bottom down and interleaved
                    end
                end
                filename = [filename,'_stc'];
                % Updating the filename only if we slice time corrected
                fprintf('Done\n')
            else
                fprintf('Skipped\n')
            end
            % This whole if statement is conducting the slice time
            % correction if it was specified by the user in the way that it
            % was specified by the user. If there are further complications
            % in slice time correction, the code above can be adjusted to
            % accomodate the way the scans were collected.
            % ___________________________________________________________ %
            
            
            
            % Motion correction of functional scans _____________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Motion correction ... '])
            cmd = [fsldir,'mcflirt -in ',filename,' -out ', ...
                filename,'_mc -refvol 0 -plots'];
            system(cmd);
            filename = [filename,'_mc'];
            % Conducting motion correction on the functional scans and
            % saving all the motion paramters. Below, I will plot out the
            % motion parameters as an output.
            par = load([filename,'.par']);
            fig = figure('visible','off','position',get(0,'screensize'));
            % Creating a hidden figure that will eventually be saved as a
            % JPEG file for the user's convenience.
            subplot(2,1,1);
            hold on;
            plot(par(:,1));
            plot(par(:,2));
            plot(par(:,3));
            title('Rotation in radians');
            legend('X','Y','Z');
            set(gca,'ylim',[-0.2 0.2],'fontsize',12);
            ylabel('Rotations (rad)');
            xlabel('Timepoints');
            % Plotting the rotation paramters
            subplot(2,1,2);
            hold on
            plot(par(:,4));
            plot(par(:,5));
            plot(par(:,6));
            title('Transformation in mm');
            legend('X','Y','Z');
            set(gca,'ylim',[-3 3],'fontsize',12);
            ylabel('Transformations (mm)');
            xlabel('Timepoints');
            % Plotting the transformation parameters
            saveas(fig,[filename,'_visual.jpg']);
            % Saving
            fprintf('Done\n')
            % This series of steps is using FSL's MCFLIRT tool to conduct
            % motion correction and then plotting out and saving the motion
            % parameters as a separate JPEG file that can be viewed later.
            % ___________________________________________________________ %
            
            
                        
            % Functional Brain Registration to Atlas ____________________ %
            % 
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Registration to atlas space ... '])
            cmd = [fsldir,'flirt -in ',filename,' -ref ', ...
                atlasdir,'/atlas_t1_brain -omat f_to_atlas.mat'];
            system(cmd);
            cmd = [fsldir,'flirt -in ',filename,' -ref ', ...
                atlasdir,'/atlas_t1_brain -applyxfm -init ', ...
                'f_to_atlas.mat -out ',filename,'_regatlas'];
            system(cmd);
            cmd = [fsldir,'fslmaths ',filename,'_regatlas -mas ', ...
                atlasdir,'/atlas_t1_brain_mask ',filename,'_regatlas'];
            system(cmd);
            system('rm f_to_atlas.mat');
            filename = [filename,'_regatlas'];
            % Using FSL's flirt to register ther functional scan to the
            % atlas brain image
            fprintf('Done\n')
            % ___________________________________________________________ %
            
            
            
            % Smoothing of Functional Scan ______________________________ %
            % 
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Spatial smoothing ... '])
            cmd = [fsldir,'fslmaths ', filename, ' -s ', num2str(sm), ...
                ' ',filename,'_sm'];
            system(cmd);
            filename = [filename,'_sm'];
            % Using FSL's fslmaths tool to smooth the functional scans by
            % however much the user specified
            fprintf('Done\n')
            % ___________________________________________________________ %
            
            

            % Temporal Filtering of functional data _____________________ %
            % 
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Temporal filtering ... '])
            % Starting this step, we are no longer using FSL for
            % preprocessign. Rather, we're loading the nifti files into
            % MATLAB and performing some steps on them directly. This is
            % why the NIfTI toolbox provided with this function is
            % critical.
            nii = load_untouch_nii([filename,'.nii.gz']);
            f = nii.img;
            [X,Y,Z,T] = size(f);
            % Getting the dimensions of the functional scan
            f = reshape(f,[X*Y*Z,T]); 
            f = f';
            % Reshaping the scan into a 2D matrix of time over voxels
            tc = f;
            % Renaming the functional timeseries as tc
            SF = 1/TR';
            % Calculating the sampling frequency using the TR
            x = size(tc,1);
            % x is the number of timepoints
            tc_mirror = tc(end:-1:1,:);
            % This the mirror of the timeseries
            tc_extended = [tc_mirror; tc; tc_mirror];
            % Extend the timecourse in a periodical way
            tc_fft = fft(tc_extended);
            % Getting the fast fourier transform
            if fil(1) > 0
                high_cutoff = round(size(tc_extended,1)*fil(1)/SF);
                tc_fft(1:high_cutoff,:) = 0;
            end
            % Conducting the high pass filtering
            if fil(2) > 0
                low_cutoff = round(size(tc_extended,1)*fil(2)/SF);
                tc_fft(low_cutoff:end,:) = 0;
            end
            % Conducting the low pass filtering
            tc_filtered = ifft(tc_fft);
            % Inverse fast fourier transform
            temp = 2*real(tc_filtered);
            tc_filtered = temp(x+1:2*x,:);
            % Getting the final filtered functional scan
            f = tc_filtered';
            f = reshape(f,[X,Y,Z,T]);
            nii.img = f;
            save_untouch_nii(nii,[filename,'_fil.nii.gz']);
            cmd = [fsldir,'fslmaths ',filename,'_fil -mas ', atlasdir, ...
                '/atlas_t1_brain_mask ',filename,'_fil'];
            system(cmd);
            % Just to cut down on file size, I'm multiplying the functional
            % scan by the brain mask again
            filename = [filename,'_fil'];
            fprintf('Done\n')
            clear f1 tc_filtered temp tc_fft low_cutoff high_cutoff
            clear tc_extended tc_mirror x
            % ___________________________________________________________ %
            
            
            
            % Global signal regression __________________________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Global signal regression ... '])
            if gsr == 1
                nii = load_untouch_nii([filename,'.nii.gz']);
                f = nii.img;
                % Loading the functional scan into Matlab
                f_gsr = f;
                % Predefining the global signal regressed functional scan
                f_mean = zeros(1,size(f,4));
                % Predefining what will me the mean of the functional scan
                for t = 1:size(f,4)
                    temp = f(:,:,:,t);
                    temp_sum = sum(temp(:));
                    temp_length = length(find(temp(:)));
                    temp_mean = temp_sum/temp_length;
                    f_mean(t) = temp_mean;
                end
                % Taking the mean of the functional timeseries. This is the
                % functional timecourse.
                f_mean_zsc = zscore(f_mean);
                % Z-scoring the mean
                for x = 1:size(f,1)
                    for y = 1:size(f,2)
                        for z = 1:size(f,3)
                            v = double(squeeze(f(x,y,z,:)));
                            beta = (f_mean_zsc * f_mean_zsc') \ ...
                                (f_mean_zsc * v);
                            f_gsr(x,y,z,:) = v' - f_mean_zsc * beta;
                        end
                    end
                end
                % Regressing the global signal from the functional
                % timeseries
                f = f_gsr; 
                filename = [filename,'_gsr'];
                nii.img = f;
                save_untouch_nii(nii,[filename,'.nii.gz']);
                % Saving the global signal regressed functional scan
                fprintf('Done\n')
            else
                fprintf('Skipped\n')
            end
            % Conducting global signal regression on the functional data if
            % the user asked for it
            % ___________________________________________________________ %
            
            
            
            % White Matter and CSF signal regression ____________________ %
            % 
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'WM/CSF signal regression ... '])
            if wmcsfr == 1
                if gsr == 0
                    nii = load_untouch_nii([filename,'.nii.gz']);
                    f = nii.img;
                    % Loading the functional scan
                end
                wm_mask = load_nii('t1_reorient_regatlas_brain_wm.nii.gz');
                wm_mask = wm_mask.img;
                % Loading white matter mask
                csf_mask = ...
                    load_nii('t1_reorient_regatlas_brain_csf.nii.gz');
                csf_mask = csf_mask.img;
                % Loading the CSF mask
                wmcsf_mask = double(logical(wm_mask + csf_mask));
                % Creating a mask of the white matter and CSF together and
                % making sure that it is a binary matrix
                f_wmcsf = double(f);
                % Predefining the matrix that will hold just the timeseries
                % from the white matter and CSF 
                for t = 1:size(f,4)
                    f_wmcsf(:,:,:,t) = f_wmcsf(:,:,:,t) .* wmcsf_mask;
                end
                % Removing all the gray matter voxels from f_wmcsf
                f_wmcsf_r = f_wmcsf;
                % Predefining the matrix that will be the regressed signal
                % from the white matter and CSF
                f_wmcsf_mean = zeros(1,size(f,4));
                % Predefining what will me the mean of the functional scan
                for t = 1:size(f,4)
                    temp = f_wmcsf(:,:,:,t);
                    temp_sum = sum(temp(:));
                    temp_length = length(find(temp(:)));
                    temp_mean = temp_sum/temp_length;
                    f_wmcsf_mean(t) = temp_mean;
                end
                % Taking the mean of the wmcsf signal. This is the
                % functional timecourse of the wmcsf.
                f_wmcsf_mean_zsc = zscore(f_wmcsf_mean);
                % Z-scoring the mean
                for x = 1:size(f,1)
                    for y = 1:size(f,2)
                        for z = 1:size(f,3)
                            v = double(squeeze(f(x,y,z,:)));
                            beta = (f_wmcsf_mean_zsc*f_wmcsf_mean_zsc') ...
                                \ (f_wmcsf_mean_zsc * v);
                            f_wmcsf_r(x,y,z,:) = ...
                                v' - f_wmcsf_mean_zsc * beta;
                        end
                    end
                end
                % Regressing the white matter and CSF signal from the
                % functional timeseries 
                f = f_wmcsf_r;
                nii.img = f;
                save_untouch_nii(nii,[filename,'_wmcsfr.nii.gz']);
                % Saving the white matter/CSF regressef functiona scan
                filename = [filename,'_wmcsfr'];
                fprintf('Done\n')
            else
                fprintf('Skipped\n')
            end
            % Conducting white matter and CSF signal regression on data in
            % case the user specified.
            % ___________________________________________________________ %
            

            
            % Z-scoring all voxels in functional scan ___________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Z-scoring ... '])
            if gsr ~= 1
                nii = load_untouch_nii([filename,'.nii.gz']);
                f = nii.img;
                % Loading the functional scan in case it wasn't already
                % loaded in the previous step
            end
            mask = load_nii([atlasdir,'/atlas_t1_brain_mask.nii.gz']);
            mask = mask.img;
            % Loading the T1 mask
            [x,y,z,~] = size(f);
            % Assigning variables to the dimensions of the scan
            for j = 1:x
                for k = 1:y
                    for l = 1:z
                        if mask(j,k,l) == 1
                            voxel = f(j,k,l,:);
                            % Working on one brain voxel at a time
                            voxel_zsc = zscore(voxel);
                            % Z-scoring the voxel signal
                            f(j,k,l,:) = voxel_zsc;
                            % Adding the zscored value to the zscored
                            % functional scan matrix
                        end
                    end
                end
            end
            nii.img = f;
            save_untouch_nii(nii,[filename,'_zsc.nii.gz']);
            filename = [filename,'_zsc'];
            clear nii x y z j k l
            fprintf('Done\n')
            % Saving the newly z-scored functional scan
            % ___________________________________________________________ %
            
            
            
            % Slicing functional images into 3D matrices ________________ %
            %
            T = datetime('now');
            time = whatsthetime(T);
            fprintf([time,sub,scan,'Slicing into 3D matrices ... '])
            if slice == 1
                [~,~,z,~] = size(f);
                % Getting dimensions of functional scan
                num_boxes = z;
                % Predefining the number of boxes the output image is going
                % to have.
                if rem(sqrt(z),1) >= eps
                    % If the brain in the z dimension is not a perfect
                    % square
                    while rem(sqrt(num_boxes),1) >= eps
                        % Well then while that's the case
                        num_boxes = num_boxes + 1;
                        % Adding to the number of boxes
                    end
                end
                slices = cell(num_boxes,1);
                % Predefining cell array that will hold timeseries for all
                % slices
                for j = 1:z
                    % Going through all slices
                    temp = squeeze(f(:,:,j,:));
                    % This is the timeseries for a slice
                    temp = rot90(temp,1);
                    % Rotating the image to make anterior face upwards
                    slices{j} = temp;
                    % Adding timeseries for each slice in each cell
                end
                [x,y,t] = size(slices{1});
                % Updating these values because of the rotation
                empty_space = zeros(x,y,t);
                % Empty space to go in empty boxes
                for j = z+1:z+(num_boxes-z)
                    % Going through the rest of the cells in slices
                    slices{j} = empty_space;
                    % Adding empty space into cell
                end
                f_sliced = zeros(x*sqrt(num_boxes),y*sqrt(num_boxes),t);
                % Predefining the output image shape
                count = 0;
                % Will be counting through the boxes
                for row = 1:sqrt(num_boxes)
                    % For every row of brain slices
                    for column = 1:sqrt(num_boxes)
                        % For every column of brain slices
                        count = count + 1;
                        % Updating the count
                        f_sliced((row*x-x+1):(row*x),(column*y-y+1): ...
                            (column*y),:) = slices{count};
                        % Organizing brain slices in rows and columns
                    end
                end
                func = f_sliced;
                filename = [filename,'_3D'];
                save([filename,'.mat'],'func','-v7.3');
                fprintf('Done\n')
            else
                fprintf('Skipped\n')
            end
            % The abve steps are slicing the 4D functional scans into 3D
            % functional scans by removing one of the spatial dimensions.
            % ___________________________________________________________ %
            
        T = datetime('now');
        time = whatsthetime(T);
        fprintf([time,sub,scan,'Finished preprocessing\n'])
        
        if nargin > 12
            string = [sub,scan,'Functional scan preprocessed'];
            p.pushNote([],string,'');
            % Keeping the user updated with progress for each scan
        end
        
        end
        % Ending everything that we were doing to each functional scan
        
end
                            
warning('on','all');

end



function time = whatsthetime(T)
if length(num2str(hour(T))) == 1
    hr = ['0',num2str(hour(T))];
else 
    hr = num2str(hour(T));
end
if length(num2str(minute(T))) == 1
    min = ['0',num2str(minute(T))];
else
    min = num2str(minute(T));
end
time = [hr,':',min,' - '];
end


% % Note: It may be the case that the machine used to run this function runs
% % out of memory during the filtering step. If that is the case consistently
% % and a machine with larger RAM is not available, then the following code
% % can replace the temporal filtering code above. It cuts the functional
% % scan in four parts, which hopefully should get over the memory problem.
% % However, this methodology is far from ideal and takes longer. Another
% % alternative is using MATLAB's built-in filters.
%             
%         % Temporal Filtering of functional data _____________________ %
%         % 
%         T = datetime('now');
%         time = whatsthetime(T);
%         fprintf([time,sub,scan,'Temporal filtering ... '])
%         nii = load_untouch_nii([filename,'.nii.gz']);
%         f = nii.img;
%         [X,Y,Z,T] = size(f);
%         f = reshape(f,[X*Y*Z,T]); f = f';
%         f1 = f(1:floor(T/4),:);
%         f2 = f(floor(T/4)+1:floor(T/2),:);
%         f3 = f(floor(T/2)+1:(floor(T/2)+floor(T/4)),:);
%         f4 = f((floor(T/2)+floor(T/4))+1:end,:);
%         save('timecourse2','f2');
%         save('timecourse3','f3');
%         save('timecourse4','f4');
%         clear f2 f3 f4
%         % To save memory in case of limited RAM, I'm going to split the
%         % functional scan into four parts and conduct filtering on each
%         % of the parts separately. This isn't ideal but very helpful.
%         tc = f;
%         clear f2 f3 f4
%         SF = 1/TR';
%         x = size(tc,1);
%         tc_mirror = tc(end:-1:1,:);
%         tc_extended = [tc_mirror; tc; tc_mirror];
%         tc_fft = fft(tc_extended);
%         high_cutoff = round(size(tc_extended,1)*(hp/SF));
%         tc_fft(1:high_cutoff,:) = 0;
%         low_cutoff = round(size(tc_extended,1)*(lp/SF));
%         tc_fft(low_cutoff:end,:) = 0;
%         tc_filtered = ifft(tc_fft);
%         temp = 2*real(tc_filtered);
%         tc_filtered = temp(x+1:2*x,:);
%         f = tc_filtered;
%         save('timecourse1','f1');
%         clear f1 tc_filtered temp tc_fft low_cutoff high_cutoff
%         clear tc_extended tc_mirror x
%         % Filtering the first quarter of the functional scan
%         load('timecourse2');
%         tc = f2;
%         clear f2
%         SF = 1/TR';
%         x = size(tc,1);
%         tc_mirror = tc(end:-1:1,:);
%         tc_extended = [tc_mirror; tc; tc_mirror];
%         tc_fft = fft(tc_extended);
%         high_cutoff = round(size(tc_extended,1)*(hp/SF));
%         tc_fft(1:high_cutoff,:) = 0;
%         low_cutoff = round(size(tc_extended,1)*(lp/SF));
%         tc_fft(low_cutoff:end,:) = 0;
%         tc_filtered = ifft(tc_fft);
%         temp = 2*real(tc_filtered);
%         tc_filtered = temp(x+1:2*x,:);
%         f2 = tc_filtered;
%         save('timecourse2','f2');
%         clear f2 tc_filtered temp tc_fft low_cutoff high_cutoff
%         clear tc_extended tc_mirror x    
%         % Filtering the second quarter of the functional scan
%         load('timecourse3');
%         tc = f3;
%         clear f3
%         SF = 1/TR';
%         x = size(tc,1);
%         tc_mirror = tc(end:-1:1,:);
%         tc_extended = [tc_mirror; tc; tc_mirror];
%         tc_fft = fft(tc_extended);
%         high_cutoff = round(size(tc_extended,1)*(hp/SF));
%         tc_fft(1:high_cutoff,:) = 0;
%         low_cutoff = round(size(tc_extended,1)*(lp/SF));
%         tc_fft(low_cutoff:end,:) = 0;
%         tc_filtered = ifft(tc_fft);
%         temp = 2*real(tc_filtered);
%         tc_filtered = temp(x+1:2*x,:);
%         f3 = tc_filtered;
%         save('timecourse3','f3');
%         clear f3 tc_filtered temp tc_fft low_cutoff high_cutoff
%         clear tc_extended tc_mirror x   
%         % Filtering the third quarter of the functional scan
%         load('timecourse4');
%         tc = f4;
%         clear f4
%         SF = 1/TR';
%         x = size(tc,1);
%         tc_mirror = tc(end:-1:1,:);
%         tc_extended = [tc_mirror; tc; tc_mirror];
%         tc_fft = fft(tc_extended);
%         high_cutoff = round(size(tc_extended,1)*(hp/SF));
%         tc_fft(1:high_cutoff,:) = 0;
%         low_cutoff = round(size(tc_extended,1)*(lp/SF));
%         tc_fft(low_cutoff:end,:) = 0;
%         tc_filtered = ifft(tc_fft);
%         temp = 2*real(tc_filtered);
%         tc_filtered = temp(x+1:2*x,:);
%         f4 = tc_filtered;
%         clear tc_filtered temp tc_fft low_cutoff high_cutoff
%         clear tc_extended tc_mirror x
%         % Filtering the fourth quarter of the functional scan
%         load('timecourse1');
%         load('timecourse2');
%         load('timecourse3');
%         delete('timecourse1.mat');
%         delete('timecourse2.mat');
%         delete('timecourse3.mat');
%         delete('timecourse4.mat');
%         f = [f1; f2; f3; f4]';
%         f = reshape(f,[X,Y,Z,T]);
%         nii.img = f;
%         save_untouch_nii(nii,[filename,'_fil.nii.gz']);
%         clear f1 f2 f3 f4
%         cmd = [fsldir,'fslmaths ',filename,'_fil -mas ', atlasdir, ...
%             '/atlas_t1_brain_mask ',filename,'_fil'];
%         system(cmd);
%         % Finally, putting the whole functional scan back together and
%         % saving it as a new file.
%         filename = [filename,'_fil'];
%         fprintf('Done\n')
%         % ___________________________________________________________ %

