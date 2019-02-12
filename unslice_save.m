function output = unslice_save(input_image, dimensions, ...
    output_image_filename, voxel_size, origin)

% _______________________________________________________________________ %
%
% During the preprocessing pipeline, the anz_preprocess script slices the
% brain so that each timepoint can be observed in one plane.
%
% This is helpful for visualization within MATLAB etc. but not great for
% visualizing in FSLview or publishing in abstracts/papers. 
%
% Hence, I wrote this script to undo that slicing and give you a 3D image
% of the brain (or in the case of functional scans, a 4D image).
%
% You can choose to output that unsliced image or you can also ask the
% function to save the outputted image as NIfTI file for you to view using
% FSLview or a similar visualization software.
%
% Inputs 
%
% input_image - this is the preprocessed, sliced image - can be either
% anatomical (2D matrix) or functional (3D matrix)
%
% dimensions - this is a row vector with the desired dimensions of your
% output image. if you don't know this, then open the original standard
% space atlas that you used during preprocessing using FSLview or a similar
% software and look at its dimensions. for example, if you used a 2mm MNI
% atlas, then your dimensions will be [91 109 91 t] with t being the number
% of timepoints in case you are unslicing a timeseries.
%
% output_image_filename - this is the filename that you would like your
% NIfTI file to have. This is an optional input, and if you do not enter
% anything here, then the script will not save a NIfTI file with your
% output image. Make sure that this is entered as a string that ends with
% '.nii'.
%
% voxel_size - this is 1x3 row vector that has the dimensions of a voxel in
% the image. to get this value, open up the standard space atlas that you
% used during preprocessing using FSLview or a similar software and look at
% the dimensions (in FSLview, you do that by clicking the 'i' button).
%
% origin - this is a 1x3 row vector that has the origin of the image. to
% get this value, open up the standard space atlas that you used during
% preprocessing using FSLview or a similar software and look at the
% origin of the image (in FSLview, these are the three values next to the
% (X,Y,Z) coordinates in the bottom left of the window).
%
% 12/12/16 - Anzar Abbas 
%
% _______________________________________________________________________ %
 
[a, b, ~] = size(input_image);
% Getting the dimensions of our image first trying functional then
% anatomical input image

x = dimensions(2);
y = dimensions(1);
% Getting the dimensions of the image that we'll need right now

rows = a/x;
columns = b/y;
% Finding out how many rows and columns our input image has 

boxes = cell(x*y,1);
% Creating a cell array that will hold all the boxes 

row_count = 0;
box_count = 0;
% We'll need these counts in the following loop

while row_count < rows 
    % For every row of brain slices 
    
    col_count = 0;
    % We'll need this count in the following loop
    
    row_count = row_count + 1;
    % Updating the row count 
    
    while col_count < columns 
        % For every column of brain slices until we reach the end of the
        % row
        
        col_count = col_count + 1;
        % Updating the column count 
        
        box_count = box_count + 1;
        % Updating the box count 
                
        temp = input_image((row_count*x-x+1):(row_count*x),...
            (col_count*y-y+1):(col_count*y),:);
        % This is the timeseries for one box from the input 
        
        temp = rot90(temp,3);
        % We rotated the images 90 degrees when we sliced them. I am
        % un-rotating them here.
        
        boxes(box_count,1) = {temp};
        % Pulling from the inputted matrix and putting into the cell array
        
    end
    
end

% Now we have each of the boxes from the inputted image into our cell

z = dimensions(3);
% This is the number of slices that we have 

output = zeros(dimensions);
% Predefining our ouput image 

for i = 1:z
    % Going through all the slices 
    
    if length(dimensions) == 4
        % In the case of a timeseries 
        
        output(:,:,i,:) = boxes{i,1};
        % Adding one box at a time to the output image matrix 
        
    end
    
    if length(dimensions) == 3
        % In the case of an anatomical image 
        
        output(:,:,i) = boxes{i,1};
        % Adding one box at a time to the output image matrix 
        
    end
    
end

if nargin > 3
    
    nii = make_nii(output, voxel_size, origin);
    save_nii(nii,output_image_filename);
    
end
