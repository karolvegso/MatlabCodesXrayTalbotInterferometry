%**************************************************************************
% this program opens M=5 images for forground and background
% measured by X-ray talbot interferometry 
% and calculate differential phase (dph) image 
% Author: Karol Vegso
% Affiliation: Institute of Physics, Slovak Academy of Sciences
%**************************************************************************
clear all
close all
%**************************************************************************
tic
% path to folder with forground images
path_to_fg_images='d:\Matlab\Xray_Talbot_Interferometry\XTI_dph_abs_vis_image_calculator_Momose_approach_CPU\';
% path to background images
path_to_bg_images='d:\Matlab\Xray_Talbot_Interferometry\XTI_dph_abs_vis_image_calculator_Momose_approach_CPU\';
% path to output folder to save final differentila phase image
path_to_output_folder='d:\Matlab\Xray_Talbot_Interferometry\XTI_dph_abs_vis_image_calculator_Momose_approach_CPU\';
% number of images in fringe scanning technique 
% the same number for forground and background
M=5;
% open forground images and store them in 3D matrix
% root image name for forground
root_image_name_fg='fg';
% root image name for background
root_image_name_bg='bg';
% number of digits in image numbering
number_digits=6;
% size of image
% size of image in horizontal direction, number of columns
image_size_cols=1536;
% size of image in vertical direction, number of rows
image_size_rows=512;
% image size
image_size=image_size_cols*image_size_rows;
% precision
% here, it is unsigned 16 bit integer
precision='uint16';
% order of reading bytes
order_read_bytes='ieee-le';
% create image buffer to store single forground image
image_buffer_fg=zeros(image_size, 1);
% create image buffer for final M forground images
% this is 3D set or stack of M images
image_store_fg=zeros(image_size_rows, image_size_cols, M);
% process M forground images
for index_0=1:M
    % create image_number as string
    image_number=num2str(index_0);
    % add number digits
    image_number=pad(image_number, number_digits, 'left');
    % replace empty spaces with zeros
    image_number=replace(image_number,' ', '0');
    % full path to forground image
    path_to_fg_image=strcat(path_to_fg_images, root_image_name_fg, image_number,'.raw');
    % create fileID
    fileID=fopen(path_to_fg_image);
    % read image
    image_buffer_fg=fread(fileID, image_size, precision, order_read_bytes);
    current_image=reshape(image_buffer_fg, [image_size_cols,image_size_rows]);
    current_image=transpose(current_image(:,:));
    % show current image
    imshow(current_image,[]);
    fclose(fileID);
    image_store_fg(:,:,index_0)=current_image(:,:);
    pause(1);
end
%**************************************************************************
% create image buffer to store single background image
image_buffer_bg=zeros(image_size, 1);
% create image buffer for final M background images
% this is 3D set or stack of M images
image_store_bg=zeros(image_size_rows, image_size_cols, M);
% process M background images
for index_0=1:M
    % create image_number as string
    image_number=num2str(index_0);
    % add number digits
    image_number=pad(image_number, number_digits, 'left');
    % replace empty spaces with zeros
    image_number=replace(image_number,' ', '0');
    % full path to background image
    path_to_bg_image=strcat(path_to_bg_images, root_image_name_bg, image_number,'.raw');
    % create fileID
    fileID=fopen(path_to_bg_image);
    % read image
    image_buffer_bg=fread(fileID, image_size, precision, order_read_bytes);
    current_image=reshape(image_buffer_bg, [image_size_cols,image_size_rows]);
    current_image=transpose(current_image(:,:));
    % show current image
    imshow(current_image,[]);
    fclose(fileID);
    image_store_bg(:,:,index_0)=current_image(:,:);
    pause(1);
end
%**************************************************************************
% calculate differential phase image
% define phase step
phase_step=(2*pi)/M;
% complex numbers image for foreground
complex_numbers_image_fg=zeros(image_size_rows, image_size_cols);
% complex numbers image for background
complex_numbers_image_bg=zeros(image_size_rows, image_size_cols);
% phase image for foreground
phase_image_fg=zeros(image_size_rows, image_size_cols);
% phase image for background
phase_image_bg=zeros(image_size_rows, image_size_cols);
% amplitude image for foreground
amp_image_fg=zeros(image_size_rows, image_size_cols);
% amplitude image for background
amp_image_bg=zeros(image_size_rows, image_size_cols);
% offset image for foreground
offset_image_fg=zeros(image_size_rows, image_size_cols);
% offset image for background
offset_image_bg=zeros(image_size_rows, image_size_cols);
% initiliaze differential phase image
dph_image=zeros(image_size_rows, image_size_cols);
% initiliaze absorption image
abs_image=zeros(image_size_rows, image_size_cols);
% initiliaze visibility image
vis_image=zeros(image_size_rows, image_size_cols);
%**************************************************************************
% start to measure time for calculation fo differential phase image
% calculate complex numbers image for forground
for index_0=1:M
    complex_numbers_image_fg(:,:)=complex_numbers_image_fg(:,:)+...
        image_store_fg(:,:,index_0)*exp(i*(index_0-1)*phase_step);
    offset_image_fg(:,:)=offset_image_fg(:,:)+image_store_fg(:,:,index_0);
end
phase_image_fg=angle(complex_numbers_image_fg(:,:));
amp_image_fg=2*abs(complex_numbers_image_fg(:,:))/M;
offset_image_fg=offset_image_fg(:,:)/M;
% calculate complex numbers image for background
for index_0=1:M
    complex_numbers_image_bg(:,:)=complex_numbers_image_bg(:,:)+...
        image_store_bg(:,:,index_0)*exp(i*(index_0-1)*phase_step);
    offset_image_bg(:,:)=offset_image_bg(:,:)+image_store_bg(:,:,index_0);
end
phase_image_bg=angle(complex_numbers_image_bg(:,:));
amp_image_bg=2*abs(complex_numbers_image_bg(:,:))/M;
offset_image_bg=offset_image_bg(:,:)/M;
% calculate final differntial phase image
dph_image=phase_image_fg(:,:)-phase_image_bg(:,:);
% calculate final absorption image
abs_image=offset_image_fg(:,:)./offset_image_bg(:,:);
% calculation final visibility image
vis_image=(amp_image_fg(:,:)./offset_image_fg(:,:))./(amp_image_bg(:,:)./offset_image_bg(:,:));
% stop to measure time for calculation fo differential phase image
% unwrap dph image
%dph_image=unwrap(dph_image(:,:));
% show final differential phase image
imshow(dph_image(:,:), []);
% show final absorption image
imshow(abs_image(:,:), []);
% show final visibility image
imshow(vis_image(:,:), []);
% reshape differential phase image
dph_image_reshaped=transpose(dph_image(:,:));
dph_image_reshaped=reshape(dph_image_reshaped(:,:), [image_size,1]);
% reshape absorption image
abs_image_reshaped=transpose(abs_image(:,:));
abs_image_reshaped=reshape(abs_image_reshaped(:,:), [image_size,1]);
% reshape visibility image
vis_image_reshaped=transpose(vis_image(:,:));
vis_image_reshaped=reshape(vis_image_reshaped(:,:), [image_size,1]);
%**************************************************************************
% define paths to save final images
% path to final differential phase image or dph image
path_to_final_dph_image=strcat(path_to_output_folder, 'dph.raw');
% path to final absorption image or abs image
path_to_final_abs_image=strcat(path_to_output_folder, 'abs.raw');
% path to final visibility image or vis image
path_to_final_vis_image=strcat(path_to_output_folder, 'vis.raw');
%**************************************************************************
% create fileID
fileID=fopen(path_to_final_dph_image, 'w');
% save differential phase image as 64bit real image with little bit endian
% ordering
fwrite(fileID,dph_image_reshaped,'double', 'ieee-le.l64');
fclose(fileID);
%**************************************************************************
% create fileID
fileID=fopen(path_to_final_abs_image, 'w');
% save absorption image as 64bit real image with little bit endian
% ordering
fwrite(fileID,abs_image_reshaped,'double', 'ieee-le.l64');
fclose(fileID);
%**************************************************************************
% create fileID
fileID=fopen(path_to_final_vis_image, 'w');
% save visilbility image as 64bit real image with little bit endian
% ordering
fwrite(fileID,vis_image_reshaped,'double', 'ieee-le.l64');
fclose(fileID);
%**************************************************************************
toc
%**************************************************************************
% end of program
%**************************************************************************