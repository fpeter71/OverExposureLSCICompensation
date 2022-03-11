%{
%%
clear
% collection of files with certain extension from a folder
path = 'C:\Users\folde\OneDrive - SZTAKI\Matlab\stereo\IllumiKez\'; % IllumiSweep_Papir IllumiSweep_3_20ms IllumiSweep_1
d = dir([path '*.tiff']); % change extension as needed
for i = 1:1 %length(d)
    name = d(i).name;
    filenames{i} = [path name];
end
% call the correction
%}
load sampleimages

N = 7;
[K_raw, K_corrected, R_saturationratio] = ...
    OverExposureCorrection({'sample.tiff'}, N, 255, 1);
% show results

subplot 131
imagesc(1./K_raw.^2,[0 3]);
title('Raw contrast map');
colorbar

subplot 132
imagesc(100 * R_saturationratio);
title('Saturation ratio [%]');
colorbar

subplot 133
imagesc(1./K_corrected.^2,[0 3]);
title('Corrected contrast map');
colorbar
%colormap(jet)
