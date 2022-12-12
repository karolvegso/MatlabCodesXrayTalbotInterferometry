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
% path to folder with forground images
path_to_fg_images='d:\Matlab\Xray_Talbot_Interferometry\XTI_dph_image_calculator_sinus_fitting_CPU\';
% path to background images
path_to_bg_images='d:\Matlab\Xray_Talbot_Interferometry\XTI_dph_image_calculator_sinus_fitting_CPU\';
% path to output folder to save final differentila phase image
path_to_output_folder='d:\Matlab\Xray_Talbot_Interferometry\2017A_pp\final_image\';
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
% generate phase coordinates for sinusidal fitting
phase=zeros(1,M);
for index_0=1:M
    phase(1,index_0)=(index_0-1)*phase_step;
end
% phase image for forground
phase_image_fg=zeros(image_size_rows, image_size_cols);
% phase image for background
phase_image_bg=zeros(image_size_rows, image_size_cols);
% initiliaze differential phase image
dph_image=zeros(image_size_rows, image_size_cols);
%**************************************************************************
% perform sinusiodal fit on the foreground data
%**************************************************************************
% maximum number of iterations
k_max=1000;
% epsilon 1
epsilon_1=1.0e-8;
% epsilon 2
epsilon_2=1.0e-8;
% tau
tau=1.0e-3;
%**************************************************************************
% perform sinusiodal fitting pixel by pixel on the foreground data
%**************************************************************************
tic
for index_0=1:image_size_rows
    for index_1=1:image_size_cols
        %tic
        y_data(1,:)=image_store_fg(index_0,index_1,:);
        % initial parameters
        % amplitude
        x01=(max(y_data(1,:))-min(y_data(1,:)))/2;
        % phase shift
        x02=pi/2;
        % offset
        x03=mean(y_data(1,:));
        % initial vector
        x0=[x01 x02 x03];
        fit_result=levmar_sinus(phase, y_data, x0, k_max, epsilon_1, epsilon_2, tau);
        phase_image_fg(index_0,index_1)=fit_result(1,2);
        %toc
    end
end
toc
%**************************************************************************
% perform sinusiodal fit on the background data
%**************************************************************************
% maximum number of iterations
k_max=1000;
% epsilon 1
epsilon_1=1.0e-8;
% epsilon 2
epsilon_2=1.0e-8;
% tau
tau=1.0e-3;
%**************************************************************************
% perform sinusiodal fitting pixel by pixel on the background data
%**************************************************************************
tic
for index_0=1:image_size_rows
    for index_1=1:image_size_cols
        %tic
        y_data(1,:)=image_store_bg(index_0,index_1,:);
        % initial parameters
        % amplitude
        x01=(max(y_data(1,:))-min(y_data(1,:)))/2;
        % phase shift
        x02=pi/2;
        % offset
        x03=mean(y_data(1,:));
        % initial vector
        x0=[x01 x02 x03];
        fit_result=levmar_sinus(phase, y_data, x0, k_max, epsilon_1, epsilon_2, tau);
        phase_image_bg(index_0,index_1)=fit_result(1,2);
        %toc
    end
end
toc
%**************************************************************************
% calculate final differential phase image
dph_image=phase_image_fg(:,:)-phase_image_bg(:,:);
% unwrap dph image
% dph_image=unwrap(dph_image(:,:));
% show final differential phase image
imshow(dph_image(:,:), []);
% reshape differential phase image
dph_image_reshaped=transpose(dph_image(:,:));
dph_image_reshaped=reshape(dph_image_reshaped(:,:), [image_size,1]);
% path to final differential phase image or dph image
path_to_final_dph_image=strcat(path_to_output_folder, 'dph.raw');
% create fileID
fileID=fopen(path_to_final_dph_image, 'w');
% save differential phase image as 64bit real image with little bit endian
% ordering
fwrite(fileID,dph_image_reshaped,'double', 'ieee-le.l64');
fclose(fileID);
%**************************************************************************
% end of program
%**************************************************************************
%**************************************************************************
% beginning of function for sinusidal fitting
%**************************************************************************
function [x_new] = levmar_sinus(t_data, y_data, x0, k_max, epsilon_1, epsilon_2, tau)
    %**************************************************************************
    % beginning of algorithm
    %**************************************************************************
    % initial iteration variable
    k=0;
    ni=2;
    M=length(t_data(1,:));
    % define Jacobian matrix
    J=zeros(M,3);
    % fill Jacobian matrix
    J(:,1)=(-1.0)*sin(t_data(1,:)+x0(1,2));
    J(:,2)=(-1.0)*x0(1,1)*cos(t_data(1,:)+x0(1,2));
    J(:,3)=-1;
    % calculate A matrix
    A(:,:)=transpose(J(:,:))*J(:,:);
    % caclulate function f
    f=zeros(1,M);
    f=y_data(1,:)-x0(1,1)*sin(t_data(1,:)+x0(1,2))-x0(1,3);
    % calculate g
    g=transpose(J(:,:))*transpose(f(1,:));
    % calculate g norm
    g_norm=sqrt(sum(g(:,1).*g(:,1)));
    % boolean variable
    found_bool=(g_norm <= epsilon_1);
    % define mi
    mi=tau*max(diag(A(:,:)));
    % initialize x vector
    x=x0(1,:);
    while ((~found_bool) & (k < k_max))
        % increase iteration by one
        k=k+1;
        B=A+mi*eye(3);
        h_lm=(-1.0)*transpose(g)*inv(B);
        h_lm_norm=sqrt(sum(h_lm(1,:).*h_lm(1,:)));
        x_norm=sqrt(sum(x(1,:).*x(1,:)));
        if (h_lm_norm <= epsilon_2*(x_norm+epsilon_2))
            found_bool=1;
        else
            x_new=x(1,:)+h_lm(1,:);
            % calculate F(x)
            % caclulate function f
            f=zeros(1,M);
            f=y_data(1,:)-x(1,1)*sin(t_data(1,:)+x(1,2))-x(1,3);
            F_x=0.5*sum(f(1,:).*f(1,:));
            % calculate F(x_new)
            f=zeros(1,M);
            f=y_data(1,:)-x_new(1,1)*sin(t_data(1,:)+x_new(1,2))-x_new(1,3);
            F_x_new=0.5*sum(f(1,:).*f(1,:));
            % ro denominator
            ro_denominator=0.5*h_lm(1,:)*(mi*transpose(h_lm(1,:))-g(:,1));
            % calculate ro - gain ratio
            ro=(F_x-F_x_new)/ro_denominator;
            if (ro > 0)
                x=x_new(1,:);
                % define Jacobian matrix
                J=zeros(M,3);
                % fill Jacobian matrix
                J(:,1)=(-1.0)*sin(t_data(1,:)+x(1,2));
                J(:,2)=(-1.0)*x(1,1)*cos(t_data(1,:)+x(1,2));
                J(:,3)=-1;
                % calculate A matrix
                A(:,:)=transpose(J(:,:))*J(:,:);
                % caclulate function f
                f=zeros(1,M);
                f=y_data(1,:)-x(1,1)*sin(t_data(1,:)+x(1,2))-x(1,3);
                % calculate g
                g=transpose(J(:,:))*transpose(f(1,:));
                % calculate g norm
                g_norm=sqrt(sum(g(:,1).*g(:,1)));
                % boolean variable
                found_bool=(g_norm <= epsilon_1);
                % calculate mi
                mi=mi*max([(1/3) (1-(2*ro-1)^3)]);
                % define ni
                ni=2;
            else
                % calculate mi
                mi=mi*ni;
                % calculate ni
                ni=2*ni;
            end
        end    
    end
    %**************************************************************************
    % end of algorithm 
    %**************************************************************************
end
%**************************************************************************
% End of function for sinusiodal fit
%**************************************************************************