clear;
clc;

imgPath = "E:\Projects and Work\SFU\Forde Lab\AutoSmarTrace\40 Lp\Noisy300dots.jpg";

A = imread(imgPath);
if size(A, 3) == 3
    Agray = rgb2gray(A);
else
    Agray = A;
end

se = strel('disk', 12);

bg = imopen(Agray, se);

tophat = Agray - bg;

im1 = medfilt2(tophat, [3 3]);

figure('Name', 'White tophat demo', 'NumberTitle', 'off');
subplot(2,2,1), imshow(Agray), title('Original Grayscale');
subplot(2,2,2), imshow(bg), title('Background (opening)');
subplot(2,2,3), imshow(tophat), title('Top-Hat (A - background)');
subplot(2,2,4), imshow(im1), title('Median filter - remove salt & pepper noise');



     
