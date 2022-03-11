clear

% collection of files with certain extension from a folder
path = 'folder of images\';
d = dir([path '*.tiff']); % change extension as needed
for i = 1:length(d)
    name = d(i).name;
    filenames{i} = [path name];
end
% or just a single image as
% filenames = {'sample.tiff'};

N = 7; % used neighborhood for NxN spatial contrast
[K_raw, K_corrected, R_saturationratio] = ...
    OverExposureCorrection(filenames, N, 255);

% show results
subplot 131
imagesc(1./K_raw.^2,[0 3]);
title('Raw 1/contrast^2 map');
colorbar

subplot 132
imagesc(100 * R_saturationratio);
title('Saturation ratio [%]');
colorbar

subplot 133
imagesc(1./K_corrected.^2,[0 3]);
title('Corrected 1/contrast^2 map');
colorbar
