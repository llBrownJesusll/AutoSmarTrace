% TestPointProcess.m
% This script loads the required files and runs the visualization.

clear; clc; close all;

% --- SETUP ---
% 1. Load the image to be analyzed
try
    input_image = imread('E:\Projects and Work\SFU\Forde Lab\AutoSmarTrace\40 Lp\Noisy300dots.jpg');
    % Convert to grayscale if it's a color image
    if size(input_image, 3) == 3
        input_image = rgb2gray(input_image);
    end
catch
    error('Could not load "Noisy300dots.jpg". Make sure it is in the same folder.');
end

% 2. Load the pre-trained neural network
try
    net_data = load('NetworkGen2.mat');
    net = net_data.net;
catch
    error('Could not load "NetworkGen2.mat". Make sure it is in the same folder.');
end


% --- EXECUTION ---
% Call the modified function to process the image and visualize the steps.
% The 'chains' variable will contain the final traced data.
disp('Starting chain detection process...');
disp('Press any key in the Command Window to advance to the next visualization step.');
chains = PointProcess_visual(input_image, net);

disp('Chain detection complete.');