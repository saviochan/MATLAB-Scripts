% Script finds firing frequency of VC'ed cell attached recordings
%
% Written by Arin Pamukcu
% Last editted August 22, 2016

clear all; close all;


%filename = '18912002Z.abf';
filename = '19124008.abf';

% Define time for finding resting potential
starttime = 4800000 %3000000;
endtime = 6000000 %5000000;
time = starttime:endtime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Acquire data into <data> variable, plot
d=abfload(filename);
dnew=d(time);
plot(d);

% Find peaks
%[pks,locs] = findpeaks(dnew,'MinPeakDistance',100,'MinPeakHeight',-10);
[pks,locs] = findpeaks(dnew,'MinPeakDistance',100,'MinPeakHeight',5);
figure; plot(time,dnew,time(locs),pks,'or')
xlabel('Time (ms)');
ylabel('Current (pA)');
% axis tight

% Rough frequency
timeinsec = length(time) / 10000; % convert time to seconds
freq1 = length(pks) / timeinsec;

% Find mean interval between peaks
cycles = diff(locs);
meanCycle = mean(cycles);
StdevCycle = std(cycles);
freq2 = 1/(meanCycle*0.0001);

% Add rate at each peak
FracCycles = 0.0001.*cycles;
instantrate = 1./FracCycles;

figure; plot(instantrate);


CVcycles = StdevCycle / meanCycle;

% Values to copy and paste on excel or gSheets
copy2excel = {'peaks per time','mean peak interval';freq1,freq2};