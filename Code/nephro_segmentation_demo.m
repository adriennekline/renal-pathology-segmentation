%==========================================================================
%   Segmentation and quantification demo code for Justin chun lipid droplet staining
%
%   Implemented by Adrienne Kline
%   University of Calgary
%   Nov. 2020
%   all rights reserved 
%--------------------------------------------------------------------------
%% import nuceli images and color images
clear all 
close all
clc
%cd('/Users/Adrienne Kline/Documents/MATLAB/Chunlab/Lipid_analysis');
%% INPUT LIPID FILE
% lipid = imread('21288_spectra_548.tif'); % load in corresponding 548 channel for lipid quantification 
% figure(1)
% imshow(lipid)
%% INPUT STRUCTURAL IMAGE
%struct = imread('19436_glom2_all_struct.tif');
allFiles = dir('SF12-11795 a_691.tif');
struct_name = allFiles.name;
struct_name2 = strsplit(struct_name,{'.','_'});
struct = imread(struct_name);
save_name = strcat(struct_name2{1,1},struct_name2{1,2});

figure(2)
imshow(struct)
[xi,yi] = getpts % double click to select center of the glom and close image

x = size(struct,1);
y = size(struct, 2);

x_position = (xi/x)*200;
y_position = (yi/y)*200;

prompt = {'Enter width of glom relative to image size (decimal format) using most conservative estimate:'};
dlgtitle = 'Glom diameter:';
dims = [2 50];
definput = {'0.48'};
answer = inputdlg(prompt,dlgtitle,dims,definput)
glom_radius_proportion = str2double(answer)/2;

%% create structuring elements based on size of glom:
answer = str2double(answer);
element = 0.48 - answer;
if element > 0
    modification = 1 - element;
    modification = round(modification);
elseif element < 0
    modification = 1 + element;
    modification = round(modification);
else element = 0
    modification = 1;
end

%% Process Structural Image 
struct_gray = rgb2gray(struct); % convert from rgb to grayscale 
figure (3)
imshow(struct_gray)

sigma = 7; % sigma for the Guassian filter 

struct_gray_Gfilt = imgaussfilt(struct_gray,sigma); % filter image with guassian filter 
figure (4)
imshow(struct_gray_Gfilt)

%% Isolate background based on limited exposure
structbw_image_blackbits = struct_gray_Gfilt <= 2; % find 'black' background parts of image
figure (5)
imshow(structbw_image_blackbits)

%% Chan Vese segmentation
% use a guassian filtered image to remove excess noise and boundary
% fluctuations
% set to a starting medium mask
% iterate 200 times
% segment out the glomerulus
glom_chan = chenvese(struct_gray_Gfilt,'medium',195,0.2,'chan',x_position, y_position, glom_radius_proportion);

seg_normalsize = imresize(glom_chan,[size(struct_gray, 1) size(struct_gray, 1)]); % convert image back to original size 
figure (7)
imshow(seg_normalsize)

% develop a structuring element with 90 radius and 4 spokes - as a circle
% shape 
se = strel('disk',90*modification, 4);
% se = se*modification;

% Close holes in mask for glom based on the segmentation using the structuring
% element
seg_areafilt_close = imclose(seg_normalsize, se);
figure (8)
imshow(seg_areafilt_close)

% erode mask for glom based on the same size structuring element - se 
seg_areafilt_close_erode = imerode(seg_areafilt_close,se);
figure (9)
imshow(seg_areafilt_close_erode)

% dilate the mask for the glom using the same size structuring element
seg_areafilt_close_erode_dilate = imdilate(seg_areafilt_close_erode,se);
figure (10)
imshow(seg_areafilt_close_erode_dilate)

% Remove background black segments and 'black' parts from mask for
% refinement
seg_areafilt_close_erode_dilate_blackremove=double(seg_areafilt_close_erode_dilate) - double(structbw_image_blackbits);
figure (11)
imshow(seg_areafilt_close_erode_dilate_blackremove)

% convert from type double back to a logical type 
seg_areafilt_close_erode_dilate_blackremove_logical = logical(seg_areafilt_close_erode_dilate_blackremove==1);
figure (12)
imshow(seg_areafilt_close_erode_dilate_blackremove_logical)

% Filter the mask based on area keeping the largest element
seg_areafilt_close_erode_dilate_blackremove_areafilt = bwareafilt(seg_areafilt_close_erode_dilate_blackremove_logical,1);
figure (13)
imshow(seg_areafilt_close_erode_dilate_blackremove_areafilt)

%showcase outline of glom segmentation
figure (15)
imshow(struct);
hold on

yellow_mask = seg_areafilt_close_erode_dilate_blackremove_areafilt;
yellow_boundaries = bwboundaries(yellow_mask);
visboundaries(yellow_boundaries, 'Color', 'y')

%% Segmenation of tubules and interstitium
% this is based on what is not glom, or background and what is left in the
%structural image

% generate the inverse image of what is not black/background (this is a
% logical image)
bw_blackbits_inv = ~structbw_image_blackbits;

% Mask our dilated glom mask on the structural gray filtered image 
tubules = bsxfun(@times, struct_gray_Gfilt, cast(~seg_areafilt_close_erode_dilate,'like', struct_gray_Gfilt));
figure (17)
imshow(tubules)

% Perform Chane Vese segmentation on the whole image
% gets rid of glom
% gets rid of background/black
tubules_and_interstitium_gray = bsxfun(@times, tubules, cast(~structbw_image_blackbits,'like', tubules));
figure (18)
imshow(tubules_and_interstitium_gray)

tubules_and_interstitium_bw = im2bw(tubules_and_interstitium_gray,0);
figure (19)
imshow(tubules_and_interstitium_bw)

% Threshold to identify brighter tubules
bw_tubulesonly = im2bw(tubules_and_interstitium_gray,0.075);
figure (20)
imshow(bw_tubulesonly)

struc_element = strel('disk',10*modification, 4);

tubules_mask_close = imclose(bw_tubulesonly, struc_element);
figure (21)
imshow(tubules_mask_close)

% Cast mask for tubules on struct image
figure (22)
imshow(struct);
hold on
cyan_mask = tubules_mask_close;
cyan_boundaries = bwboundaries(cyan_mask);
visboundaries(cyan_boundaries, 'Color', 'c')
%tubules_extracted = bsxfun(@times, struct, cast(bw,'like', struct));

% segment out intersitium based on what is left that is not glom and not
% tubule!
% remove tubule mask from tubules and interstitium mask
interstit_only_mask_double = tubules_and_interstitium_bw - tubules_mask_close;
interstit_only_mask_logical = logical(interstit_only_mask_double==1);
figure (23)
imshow(interstit_only_mask_logical)

figure (24)
imshow(struct);
hold on
green_mask = interstit_only_mask_logical;
green_boundaries = bwboundaries(green_mask);
visboundaries(green_boundaries, 'Color', 'g')


%% Lipid Channel Processing 
gray_lipid = rgb2gray(lipid); % convert from rgb to binary 
struct_bw_lipid = im2bw(lipid,0);
figure (25)
imshow(struct_bw_lipid)

lipid_glom = bitand(struct_bw_lipid,seg_areafilt_close_erode_dilate_blackremove_areafilt);
figure (26)
imshow(lipid_glom)
title('Lipids in Glom')
hold on
yellow_mask = seg_areafilt_close_erode_dilate_blackremove_logical;
yellow_boundaries = bwboundaries(yellow_mask);
visboundaries(yellow_boundaries, 'Color', 'y')
hold off

lipid_tubules = bitand(struct_bw_lipid,tubules_mask_close);
figure (27)
imshow(lipid_tubules)
title('Lipids in tubules')
hold on
cyan_mask = tubules_mask_close;
cyan_boundaries = bwboundaries(cyan_mask);
visboundaries(cyan_boundaries, 'Color', 'c')
hold off

lipid_interstit = bitand(struct_bw_lipid,interstit_only_mask_logical);
figure (28)
imshow(lipid_interstit)
title('Lipids in interstitium')
hold on
green_mask = interstit_only_mask_logical;
green_boundaries = bwboundaries(green_mask);
visboundaries(green_boundaries, 'Color', 'g')
hold off

area_of_glom = nnz(seg_areafilt_close_erode_dilate_blackremove_logical);
area_of_tubules = nnz(tubules_mask_close);
area_of_interstitium = nnz(interstit_only_mask_logical);

area_of_lipid_glom =nnz(lipid_glom);
area_of_lipid_tubules = nnz(lipid_tubules);
area_of_lipid_interstit = nnz(lipid_interstit);

prop_lipid_in_glom = area_of_lipid_glom/area_of_glom;
prop_lipid_in_tubules = area_of_lipid_tubules/area_of_tubules;
prop_lipid_in_interstitium = area_of_lipid_interstit/area_of_interstitium;

Lipid_data = [prop_lipid_in_glom prop_lipid_in_tubules prop_lipid_in_interstitium];

% data=ones(10,4);     %Sample 2-dimensional data
data_cells = num2cell(Lipid_data);     %Convert data to cell array
col_header={'Glomerulus','Tubules','Interstitium'};     %Row cell array (for column labels)
% row_header(1:10,1)={'Time'};     %Column cell array (for row labels)
output_matrix=[col_header; data_cells];     %Join cell arrays
xlswrite(save_name,output_matrix); 