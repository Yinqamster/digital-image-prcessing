%DIP16 Assignment 2
%Edge Detection
%In this assignment, you should build your own edge detection and edge linking 
%function to detect the edges of a image.
%Please Note you cannot use the build-in matlab edge and bwtraceboundary function
%We supply four test images, and you can use others to show your results for edge
%detection, but you just need do edge linking for rubberband_cap.png.
clc; clear all;
% Load the test image
% image='rubberband_cap.png';
image='single_key.png';
% image='keys_set.png';
% image='moon.jpg';
% image='test5.jpg';
imgTest = im2double(imread(image));

imgTestGray = rgb2gray(imgTest);
figure; clf;
imshow(imgTestGray);

%now call your function my_edge, you can use matlab edge function to see
%the last result as a reference first
% img_edge_before = edge(imgTestGray);
% figure;clf;
% imshow(img_edge_before);
% img_edge = my_edge(filter_gray_image);

img_edge = my_edge(imgTestGray);

if (strcmp(image,'rubberband_cap.png'))
    figure;clf;
    % imshow(img_edge);
    background = im2bw(imgTest, 1);
    imshow(background);
    %using imtool, you select a object boundary to trace, and choose
    %an appropriate edge point as the start point 
    imtool(img_edge);
    %now call your function my_edgelinking, you can use matlab bwtraceboundary 
    %function to see the last result as a reference first. please trace as many 
    %different object boundaries as you can, and choose different start edge points.
    % Bxpc = bwtraceboundary(img_edge, [197, 329], 'N');
    % Bxpc = bwtraceboundary(img_edge, [126, 232], 'N');
    Bxpc1 = my_edgelinking(img_edge, 125, 232);
    Bxpc2 = my_edgelinking(img_edge, 300, 210);
    Bxpc3 = my_edgelinking(img_edge, 300, 210);
    Bxpc4 = my_edgelinking(img_edge, 196, 77);
    Bxpc5 = my_edgelinking(img_edge, 152, 87);
    Bxpc6 = my_edgelinking(img_edge, 93, 296 );
    hold on
    % figure;clf;
    plot(Bxpc1(:,2), Bxpc1(:,1), 'w', 'LineWidth', 1);
    plot(Bxpc2(:,2), Bxpc2(:,1), 'w', 'LineWidth', 1);
    plot(Bxpc3(:,2), Bxpc3(:,1), 'w', 'LineWidth', 1);
    plot(Bxpc4(:,2), Bxpc4(:,1), 'w', 'LineWidth', 1);
    plot(Bxpc5(:,2), Bxpc5(:,1), 'w', 'LineWidth', 1);
    plot(Bxpc6(:,2), Bxpc6(:,1), 'w', 'LineWidth', 1);
end