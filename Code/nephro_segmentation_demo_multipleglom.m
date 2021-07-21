%==========================================================================
%%   Segmentation and quantification demo code for Justin chun lipid droplet staining
%
%   Implemented by Adrienne Kline
%   University of Calgary
%   Copyright (c) 2020
%--------------------------------------------------------------------------
%% import nuceli images and color images
clear all 
close all
clc
%cd('/Users/Adrienne Kline/Documents/MATLAB/Chunlab/Lipid_analysis');
%% INPUT LIPID FILE
lipid = imread('sample13011_a_lipid.tif'); % load in corresponding 548 (lipid) channel for lipid quantification 
figure(1)
imshow(lipid)
%% INPUT STRUCTURAL IMAGE
%struct = imread('19436_glom2_all_struct.tif');
structureFile = dir('SF1219002_e_struct.tif'); % load in corresponding 691 (structure) channel for lipid quantification 
struct_name = structureFile.name;
struct_name2 = strsplit(struct_name,{'.','_'});
struct = imread(struct_name);
save_name = strcat(struct_name2{1,1},struct_name2{1,2},struct_name2{1,3});

figure(2)
imshow(struct)
title('Select all glomeruli centers')
[xi,yi] = getpts % double click to select center of the glom and close image

x = size(struct,1);
y = size(struct, 2);

for a = 1:size(xi,1)
    x_position(a) = (xi(a)/x)*200;
    y_position(a) = (yi(a)/y)*200;
end

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
    %modification = round(modification);
elseif element < 0
    modification = 1 + element;
   % modification = round(modification);
else element = 0
    modification = 1;
end

%% Process Structural Image 
struct_gray = rgb2gray(struct); % convert from rgb to grayscale 
figure (3)
imshow(struct_gray)

% Outline a tubule 
imshow(struct)
title('Sample a portion of a tubule')
tubule_mask = drawfreehand;
tubule_mask.FaceAlpha = 1;
tubule_mask.FaceSelectable = false;
tubule_mask_freehand = createMask(tubule_mask,struct_gray);

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
for e = 1:size(x_position,2)
    glom_chan = chenvese(struct_gray_Gfilt,'medium',185,0.2,'chan',x_position(e), y_position(e), glom_radius_proportion);

    seg_normalsize = imresize(glom_chan,[size(struct_gray, 1) size(struct_gray, 2)]); % convert image back to original size 
    figure (7)
    imshow(seg_normalsize)

    % develop a structuring element with 90 radius and 4 spokes - as a circle
    % shape 
    se = strel('disk', round((size(struct_gray,1)*0.0439453125)*modification), 4);

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
    
    gloms{e} = seg_areafilt_close_erode_dilate_blackremove_areafilt;
end

all_glom = zeros(size(struct,1),size(struct,1));
for o = 1:size(gloms,2)
    all_glom = all_glom + gloms{o};
end

all_glom_logical = logical(all_glom);
figure
imshow(all_glom_logical)

Glom_struct = bsxfun(@times, struct, cast(all_glom_logical,'like', struct));
figure
figure (6)
imshow(Glom_struct)


%showcase outline of glom segmentation
figure (15)
imshow(struct);
hold on

yellow_mask = all_glom_logical;
yellow_boundaries = bwboundaries(yellow_mask);
visboundaries(yellow_boundaries, 'Color', 'y')

%% Segmenation of tubules and interstitium
% this is based on what is not glom, or background and what is left in the
%structural image

% generate the inverse image of what is not black/background (this is a
% logical image)
bw_blackbits_inv = ~structbw_image_blackbits;

% Mask our dilated glom mask on the structural gray filtered image 
tubules = bsxfun(@times, struct_gray_Gfilt, cast(~all_glom_logical,'like', struct_gray_Gfilt));
%tubules = bsxfun(@times, struct_gray_Gfilt, cast(~seg_areafilt_close_erode_dilate,'like', struct_gray_Gfilt));
figure (17)
imshow(tubules)

% Perform Chane Vese segmentation on the whole image
% gets rid of glom
% gets rid of background/black
tubules_and_interstitium_gray = bsxfun(@times, tubules, cast(~structbw_image_blackbits,'like', tubules));
figure (18)
imshow(tubules_and_interstitium_gray)

tubules_and_interstitium_bw = imbinarize(tubules_and_interstitium_gray,0);
figure (19)
imshow(tubules_and_interstitium_bw)

% tubules_interstitium_gray_Gfilt = imgaussfilt(tubules_and_interstitium_gray,sigma); % filter image with guassian filter 
% figure (20)
% imshow(tubules_interstitium_gray_Gfilt)

tubule_sample = bsxfun(@times, struct_gray, cast(tubule_mask_freehand,'like', struct_gray));
figure (21)
imshow(tubule_sample)

[label, number_of_tubules] = bwlabel(tubule_mask_freehand);
[r,c] = find(label == 1);
for e = 1:size(r,1)
    tubule_gray(e) = (tubule_sample(r(e,1),c(e,1),1));
end  

clear e

mean_tubule_gray = mean(tubule_gray)/255;
std_tubule_gray = std(double(tubule_gray));
rescale_std_tubule_gray = (std_tubule_gray/255)/2;
mean_tubule_gray_adjust = mean_tubule_gray - rescale_std_tubule_gray;
% level = graythresh(tubules_interstitium_gray_Gfilt);
% level_adjust = level + 0.01;

% Threshold to identify brighter tubules
%bw_tubulesonly = im2bw(tubules_and_interstitium_gray,0.075);
bw_tubulesonly = imbinarize(tubules_and_interstitium_gray,mean_tubule_gray_adjust);
figure (21)
imshow(bw_tubulesonly)

tubule_radi = (answer/30)*(size(struct_gray,1)); 

bw_tubulesonly_areafilt = bwareafilt(bw_tubulesonly,[pi*(tubule_radi.^2) size(struct_gray,1)*size(struct_gray,2)]);
figure (22)
imshow(bw_tubulesonly_areafilt)

struc_element = strel('disk', round((size(struct_gray,1)*0.0068359375)*modification), 4);

tubules_mask_close = imclose(bw_tubulesonly_areafilt, struc_element);
figure (22)
imshow(tubules_mask_close)

% Cast mask for tubules on struct image
figure (23)
imshow(struct);
hold on
cyan_mask = tubules_mask_close;
cyan_boundaries = bwboundaries(cyan_mask);
visboundaries(cyan_boundaries, 'Color', 'c')

tubules_struct = bsxfun(@times, struct, cast(tubules_mask_close,'like', struct));
figure
imshow(tubules_struct) 
%tubules_extracted = bsxfun(@times, struct, cast(bw,'like', struct));

% segment out intersitium based on what is left that is not glom and not
% tubule!
% remove tubule mask from tubules and interstitium mask
interstit_only_mask_double = tubules_and_interstitium_bw - tubules_mask_close;
interstit_only_mask_logical = logical(interstit_only_mask_double==1);
figure (24)
imshow(interstit_only_mask_logical)

figure (25)
imshow(struct);
hold on
green_mask = interstit_only_mask_logical;
green_boundaries = bwboundaries(green_mask);
visboundaries(green_boundaries, 'Color', 'g')

interstit_struct = bsxfun(@times, struct, cast(interstit_only_mask_logical,'like', struct));
figure
imshow(interstit_struct) 


%% Lipid Channel Processing 
gray_lipid = rgb2gray(lipid); % convert from rgb to binary 
struct_bw_lipid = im2bw(lipid,0);
figure (26)
imshow(struct_bw_lipid)

lipid_glom = bitand(struct_bw_lipid,all_glom_logical);
figure (27)
imshow(lipid_glom)
title('Lipids in Glom')
hold on
yellow_mask = all_glom_logical;
yellow_boundaries = bwboundaries(yellow_mask);
visboundaries(yellow_boundaries, 'Color', 'y')
hold off

lipid_tubules = bitand(struct_bw_lipid,tubules_mask_close);
figure (28)
imshow(lipid_tubules)
title('Lipids in tubules')
hold on
cyan_mask = tubules_mask_close;
cyan_boundaries = bwboundaries(cyan_mask);
visboundaries(cyan_boundaries, 'Color', 'c')
hold off

lipid_interstit = bitand(struct_bw_lipid,interstit_only_mask_logical);
figure (29)
imshow(lipid_interstit)
title('Lipids in interstitium')
hold on
green_mask = interstit_only_mask_logical;
green_boundaries = bwboundaries(green_mask);
visboundaries(green_boundaries, 'Color', 'g')
hold off

area_of_glom = nnz(all_glom_logical);
area_of_tubules = nnz(tubules_mask_close);
area_of_interstitium = nnz(interstit_only_mask_logical);

area_of_lipid_glom =nnz(lipid_glom);
area_of_lipid_tubules = nnz(lipid_tubules);
area_of_lipid_interstit = nnz(lipid_interstit);

prop_lipid_in_glom = area_of_lipid_glom/area_of_glom;
prop_lipid_in_tubules = area_of_lipid_tubules/area_of_tubules;
prop_lipid_in_interstitium = area_of_lipid_interstit/area_of_interstitium;

Lipid_data = [prop_lipid_in_glom prop_lipid_in_tubules prop_lipid_in_interstitium];

%% Output to Excel
data_cells = num2cell(Lipid_data);     %Convert data to cell array
col_header={'Glomerulus','Tubules','Interstitium'};     %Row cell array (for column labels)
output_cell=[col_header; data_cells];%Join cell arrays
save_name = strcat(save_name, '.xlsx');
writecell(output_cell, save_name); 



