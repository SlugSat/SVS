clc;
clear;
close all;

%Read the image and converts it to an RGB matrix
P=imread('Uniformity_Experiment.jpg');

Gray=rgb2gray(P); %Converts image to gray scale to find relative intensity
[row column]=size(Gray);
myMeanFunction=@(block_struct) mean2(block_struct.data); %Assigns mean2 as a function for image processing
blockMean=blockproc(Gray,[ceil(row/10) ,ceil(column/10)], myMeanFunction); %Uses mean2 for 1/10 of the image

%Finds the Maximum and Minimum Image Intensity on each block.The
%Non-uniformity is printed
max = max(max(blockMean));
min=min(min(blockMean));
fprintf('The Non-Uniformity is %.1f%%\n',((max-min)./(max+min)).*100);

%Display the relative image intensity
figure(1);
surf(blockMean);
view(2);
colorbar
xlabel('X-Axis Block');
ylabel('Y-Axis Block');
zlabel('Average Image Intensity');


