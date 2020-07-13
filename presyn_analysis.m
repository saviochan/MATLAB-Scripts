%%
clear; clc; close all;

%%
start_baseline_range = (150:1300);
access_peak1_range = (1310:1340);
access_peak2_range = (2300:2400);
long_baseline_range = (2645:7720);
epsp1peak_range = (7930:8350);
epsp2peak_range = (8440:8700);

%%
current_dir = pwd;
[FileName,PathName] = uigetfile('*.abf','SELECT THE .abf FILE');
cd(current_dir);
[D,si,h] = abfload(FileName);


%%
% % Get all .abf files as a structure
% abf_files = dir('*.abf');
% num_files = length(abf_files);
% 
% if num_files == 0, 
%     error('No .abf files found.'); 
% end


%%
% Reshape the file into a 2-dimensional matrix
[a,b,c]=size(D);
d = reshape(D,a,c); % Get rid of headstage channel dimension

% % Change here to set range of traces.
% d=D(:,1:18);

%%
mean_start_baseline = mean(d(start_baseline_range),2);
%mean_access_peak2 = mean(d(access_peak2_range),2);
mean_long_baseline = mean(d(long_baseline_range),2);
%mean_epsp1peak = mean(d(epsp1peak_range),2);
%mean_epsp2peak = mean(d(epsp2peak_range),2);


%%

numsweeps = size(d,2); % Determine the number of sweeps
ephys_measures = zeros(numsweeps,4); % Pre-allocate memory

for n = 1:numsweeps
    [acess_peak2loc, access_peak2mag] = peakfinder(d(access_peak2_range, n))
    [epsp1loc, epsp1mag] = peakfinder(d(epsp1peak_range, n), 10, 0, -1, true, false)
    [epsp2loc, epsp2mag] = peakfinder(d(epsp2peak_range, n), 10, 0, -1, true, false)
    
    acess_mag = max(access_peak2mag) - mean(d(start_baseline_range, n))
    epsp1_mag = min(epsp1mag) - mean(d(long_baseline_range, n))
    epsp2_mag = min(epsp2mag) - mean(d(long_baseline_range, n))
    ppr = (epsp2_mag/epsp1_mag)
    
    ephys_measures(n,:) = [acess_mag, epsp1_mag, epsp2_mag, ppr]
    
end

%% Save data output
FileName(1:8) % Display Filename
xlswrite(strcat(FileName,'_ephys-measures.xls'), ephys_measures)
  