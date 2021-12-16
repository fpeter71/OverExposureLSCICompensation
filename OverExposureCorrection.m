function [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation, steps)
%OverExposureCorrection Correction of laser speckle contrast calculation of
%oversatured image sequences
%   [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation, steps) sweeps all 
%   images found in the file list, calculates their speckle contrast of NxN
%   neighborhood and average for the sequence. K_raw is the standard
%   calulation, K_corrected is corrected for overexposure.
%   R_saturationratio [0..1] is the ratio of saturated pixels in the sliding
%   NxN window averaged for the sequence. The number of iterations (steps)
%   can be 1 or 2. One iteration quicker, less noisy, may fail at very low
%   contrast value and high saturating pixel ratio (~30%). Two iterations better,
%   with more numerical noise at high saturation ratio.
%
%   The image files should be of full path, single channel of type double,
%   uint8 or uint16. Saturation is calculated by the given I_saturation (e.g.
%   8-bit images 255, double normalized 1.0, etc.)
%
%   Example
%   --------
% cls
% clear 
% % collection of files with certain extension from a folder
% path = 'path to image sequence\';
% d = dir([path '*.tiff']); % change extension as needed
% for i = 1:length(d)
%     name = d(i).name;
%     filenames{i} = [path name];
% end
% 
% % call the correction
% N = 11;
% IterationSteps = 2; % 1 or 2
% I_saturation = 65500; % 255 for uint8, 2^16-1 for uint16, 1.0 for double. Adjust to actual camera maximum.
% [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation, IterationSteps);
%
% % show results
% figure
% subplot 131
% imagesc(1./K_raw.^2,[1 15]);
% title('Raw contrast map');
% colorbar
% 
% subplot 132
% imagesc(100 * R_saturationratio);
% title('Saturation ratio [%]');
% colorbar
% 
% subplot 133
% imagesc(1./K_corrected.^2,[1 15]);
% title('Corrected contrast map');
% colorbar
% colormap(parula)

%   Copyright 2021 SZTAKI, INSTITUTE FOR COMPUTER SCIENCE AND CONTROL
%   Peter Foldesy, Mate Siket, Adam Nagy, Imre Janoki, foldesy@sztaki.hu

disp('Correction of overexposure in laser speckle contrast imaging');

% threshold values
I_saturation_0 = 1.0 * I_saturation;
I_saturation_1 = 0.8 * I_saturation;
I_saturation_2 = 0.6 * I_saturation;

% contrast maps
K_0 = 0;
K_1 = 0;
K_2 = 0;

% saturation ratios
R_0 = 0;
R_1 = 0;
R_2 = 0;

% calculating contrasts
disp('Calculating contrast maps');
ImageNumber = length(filenames);
for i = 1:ImageNumber
    % fist interation
    imagein = imread(filenames{round(i)});
    imagein = double(imagein);
    
    % artificial saturation and saturation ratio
    R_0 = R_0 + conv2(imagein >= I_saturation_0, ones(N)/N/N,'same');
    R_1 = R_1 + conv2(imagein >= I_saturation_1, ones(N)/N/N,'same');
    if steps > 1
        R_2 = R_2 + conv2(imagein >= I_saturation_2, ones(N)/N/N,'same');
    end
    
    % contrast calculation with the given thresholds
    imagein( imagein >= I_saturation_0 ) = I_saturation_0;
    m = conv2(imagein, ones(N)/N/N,'same');
    st = stdfilt(imagein,ones(N));
    K_0 = K_0 + st./m;
    
    imagein( imagein >= I_saturation_1 ) = I_saturation_1;
    m = conv2(imagein, ones(N)/N/N,'same');
    st = stdfilt(imagein,ones(N));
    K_1 = K_1 + st./m;
   
    if steps > 1
        imagein( imagein >= I_saturation_2 ) = I_saturation_2;
        m = conv2(imagein, ones(N)/N/N,'same');
        st = stdfilt(imagein,ones(N));
        K_2 = K_2 + st./m;
    end
end

% normalization for the sequence
R_0 = R_0 / ImageNumber;
R_1 = R_1 / ImageNumber;
if steps > 1
    R_2 = R_2 / ImageNumber;
end

K_0 = K_0 / ImageNumber;
K_1 = K_1 / ImageNumber;
if steps > 1
    K_2 = K_2 / ImageNumber;
end

disp('Correction iterations');
% working with kappa^2
K_0 = K_0.^2;
K_1 = K_1.^2;
if steps > 1
    K_2 = K_2.^2;
end

% iterations
K_A = K_0 - R_0 .* (K_1 - K_0)./(R_1 - R_0);
if steps > 1
    K_B = K_1 - R_1 .* (K_2 - K_1)./(R_2 - R_1);
    K_C = K_A - R_0 .* (K_B - K_A)./(R_1 - R_0);
    % at near saturated neighborhood the K_B-K_A may change sign due to dK/dR
    % inflection point
    K_C = abs(K_C);   
else
    K_A = K_A .*(1+R_0);
end

% correction
if steps > 1
    K_D = K_C + R_0./(1 + K_C).*K_C.^2;
else
    K_D = K_A + R_0./(1 + K_A).*K_A.^2;
end
        
% merging with the orginal contrast map
K_D(isnan(K_D)) = K_0(isnan(K_D)); % elimianting possible inf, and nan errors
K_D(isinf(K_D)) = K_0(isnan(K_D));

% optional step to inpainting diverging values
K_D = regionfill(K_D, K_D < 0.01);

tobereplaced = double(R_0 > 0);
K_corrected = tobereplaced.*K_D + (1-tobereplaced).*K_0;

% return to kappa, contrast instead of kappa^2
K_corrected = sqrt(K_corrected);
K_raw = sqrt(K_0);
R_saturationratio = R_0;

disp('Ready');


