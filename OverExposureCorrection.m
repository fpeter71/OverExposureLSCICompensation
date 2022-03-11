function [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation)
%OverExposureCorrection Correction of laser speckle contrast calculation of
%oversatured image sequences
%   [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation) sweeps all 
%   images found in the file list, calculates their speckle contrast of NxN
%   neighborhood and average for the sequence. K_raw is the standard
%   calulation, K_corrected is corrected for overexposure.
%   R_saturationratio [0..1] is the ratio of saturated pixels in the sliding
%   NxN window averaged for the sequence. 
%   The correction factors of the linear extrapolation are set for low
%   contrast range.
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
% I_saturation = 65500; % 255 for uint8, 2^16-1 for uint16, 1.0 for double. Adjust to actual camera maximum.
% [K_raw, K_corrected, R_saturationratio] = OverExposureCorrection(filenames, N, I_saturation);
%
% % show results
% figure
% subplot 131
% imagesc(1./K_raw.^2,[1 15]);
% title('Raw 1/contrast^2 map');
% colorbar
% 
% subplot 132
% imagesc(100 * R_saturationratio);
% title('Saturation ratio [%]');
% colorbar
% 
% subplot 133
% imagesc(1./K_corrected.^2,[1 15]);
% title('Corrected 1/contrast^2  map');
% colorbar
% colormap(parula)

%   Copyright 2022 SZTAKI, INSTITUTE FOR COMPUTER SCIENCE AND CONTROL
%   Peter Foldesy, Mate Siket, Adam Nagy, Imre Janoki, foldesy@sztaki.hu

disp('Correction of overexposure in laser speckle contrast imaging');

% threshold values
I_saturation_0 = 1.0 * I_saturation;

% contrast maps
K_0 = 0;

% saturation ratios
R_0 = 0;

% calculating contrasts
disp('Calculating contrast maps');
ImageNumber = length(filenames);
for i = 1:ImageNumber
    % fist interation
    imagein = imread(filenames{round(i)});
    imagein = double(imagein);
    
    % artificial saturation and saturation ratio
    R_0 = R_0 + conv2(imagein >= I_saturation_0, ones(N)/N/N,'same');
    
    % contrast calculation with the given thresholds
    imagein( imagein >= I_saturation_0 ) = I_saturation_0;
    m = conv2(imagein, ones(N)/N/N,'same');
    st = stdfilt(imagein,ones(N));
    K_0 = K_0 + st./m;
end

% normalization for the sequence
R_0 = R_0 / ImageNumber;
K_0 = K_0 / ImageNumber;


disp('Correction iterations');
K_C = K_0 ./ (1-R_0 + eps);

c1 = -0.8;
q1 = -0.85;
q2 = 0.25;
K_C = K_C .* (1+c1*R_0)./(1+q1*R_0+q2*R_0.^2);
        
% merging with the orginal contrast map
K_C(isnan(K_C)) = K_0(isnan(K_C)); % elimianting possible inf, and nan errors
K_C(isinf(K_C)) = K_0(isnan(K_C));

% optional step to inpainting diverging values
K_C = regionfill(K_C, K_C < 0.01);

tobereplaced = double(R_0 > 0);

% return to kappa, contrast instead of kappa^2
K_corrected = tobereplaced.*K_C + (1-tobereplaced).*K_0;
K_raw = K_0;
R_saturationratio = R_0;

disp('Ready');

