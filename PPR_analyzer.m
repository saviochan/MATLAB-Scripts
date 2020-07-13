% PPR_analyzer determines the decay times of stimulus evoked IPSCs over a 
% user specified range of either 5% - 95% or 10% - 90% of baseline from 
% PClamp .abf files.
%
% This program analyzes data collected by the pClamp protocol:
% 00_E-S_Electrical_Protocol_160721_2p
% and with parameters: HP = -80 mV; 2 pulses; 5ms duration; 20 HZ 
%
% Created by Dr. Harry S. Xenias, PhD
% Northwestern University 
% Dept. of Physiology
% Chicago, IL 60611
% 
%
%
% Last modified: February 16, 2017
% Updated: 

%% **************************************************
%  **************************************************
% %% CHECK HARD CODED PARAMETER LISTS !!!
% 
% % LINE 43: for # of pulses
% % LINE 53: for FileName & trace range
% % LINE 108: for index selection
%
%  **************************************************
%  **************************************************

%% Clear command screen, close figures, and clear data.
clc; close all; clear;


%% Specify the sampling rate, stimulation frequency, duration, and ISI
fprintf('PLEASE VERIFTY THAT THESE DEFAULT VARIABLE SETTINGS\n');
fprintf('ARE CORRECT AND CHANGE IN THE SCRIPT AS NEEDED:\n');

% ! ! ! CHANGE THESE PARAMETERS HERE AS NEEDED ! ! !
sampling_rate = 10^4 % This is the *defalt* sampling rate used in most protocols.
stim_freq = 20 % In units of Hz.
ISI = 1/stim_freq * sampling_rate % This is the inter-stimulus interval in sampling points.
stim_duration = 5 * 10^-3 % In units of ms.

% See "Isolate individual IPSCs" section below where the parameters
% start_point, pulsenums, and cycles are used:

%%

%start_point = 3910; % For the 2.5 s protocol
start_point = 16801 % For the 5 s protocol%
pulsenums = [2];
%pulsenums = [1 2 5 10]
cycles = length(pulsenums)

fprintf('\nHit any key to continue...\n');

%pause;

%% Select .abf file, plot individual traces and mean of traces. 
[d,si,h,FileName] = selectabffile; % Call the function selectabffile.
% % Use the next two lines to hardcode load a data file and restrict the
% % range:
%F = abfload('15d01004.abf');
F = abfload(FileName);
%d = F(:,1,1:15);
d = F(:,1,1:end);
[numtraces, M] = ftraceplot(d); % Use function ftraceplot to plot traces.
numtraces

%% Select the starting point of the first stimulus (optional GUI)

%prompt = 'What is starting point of the first stimulus (in units of sampling points? ';
%start_point = input(prompt)

fprintf('\nSelect the starting point of the first stimulus before the first amplitude response.\n');
%x = floor(fptselect(M));
%start_point = x(1);

%% Select and calculate the baseline
fprintf('\nSelect the offset segment to be subtracted.\n');

% x1 = floor(fptselect(M));
% baseline_start = x1(1);
% x2  = floor(fptselect(M,baseline_start));
% baseline_end = x2(1);
baseline_start = 2200;
baseline_end = 4700;

close;

D = squeeze(d(:,1,:)); % Reduce this 3D matrix "page" to a 2D matrix.
Dzeroed = bsxfun(@minus,D,mean(d(baseline_start:baseline_end,:))); % Subtract baselines of traces from each full trace.
figure, plot(Dzeroed);
title('Offset Subtraction of Traces');

%% Select end points for transients to plot access resistance over time

% Calculate the mean of the subtracted offset traces
Mzeroed = mean(Dzeroed,2); 

% Select the segement of the *negative transient* to track access
fprintf('\nSelect the data segment containing the negative peak transient.\n');
% x1 = floor(fptselect(Mzeroed));
% startpt = x1(1);
% x2 = floor(fptselect(Mzeroed,startpt));
% endpt = x2(1);
startpt = 1150;
endpt = 1550;

close;

% Truncate the data to only the selected segment above
R = Dzeroed(startpt:endpt,:);

% Get the mins over the selected segment, which will be the peak access
Rmins = min(R,[],1);
 
Rmean = mean(Rmins);

% Find the indices of the values that lie about the mean within 20% of the max-min diff
index = find((Rmins <= Rmean - 0.2*(min(Rmins) - max(Rmins)) & (Rmins >= Rmean + 0.2*(min(Rmins) - max(Rmins)))));
index = (1:numtraces);
%index = (1:40);
S = ['\nThere were', ' ', num2str(length(index)),' traces within 20 percent of the mean amplitude.\n'];
fprintf(S);
R20values = Rmins(index);

% Note the above is within 20% of the max-min difference from the mean
% transient amplitude
figure, plot(Dzeroed(:,index));
title('Traces with access within 20% of mean transient amplitude'); 

figure
line([1 numtraces], [Rmean Rmean],'Color','g','LineWidth',2.5)
hold on;

% The traces within the 20% criteria range
N = NaN(1,numtraces);
N(1,index) = R20values;
Y = NaN(1,numtraces);
Y(1,index) = Rmins(index);
scatter(1:numtraces,Y,'k');

% The traces outside the 20% criteria range
a = NaN(1,numtraces);
i = find(isnan(N));
a(i) = Rmins(i);
hold on,
scatter(1:numtraces,a,'r');

axis([1 length(Rmins) 2*Rmean 0]);
title('Peak Transient Amplitudes vs. Trace Number.');
xlabel('Peak Transient Amplitude');
ylabel('Trace Number');

% Plot mean trace of only sweeps within 20% of mean transient amplitude
figure, plot(mean(Dzeroed(:,index),2));
title('Average trace with 20 percent transient error criterion.') 

   %% MANUAL OPTION here to change the index of the selected traces
    % Change here the index after seeing the transients if so desired:
    %I = 1:numtraces;
    
%% Isolate individual IPSCs
% Reiteratively go through each indexed sweep wtih the 20 percent error
% criterion to then isolate individual IPSCs

M = mod(index, cycles);
maxdecay = 0; % For determining the maximum decay for xlim settings (below)
maxylim = 0; % For determing the maximum ylim for later axis setting (below)

for cycle = 1:cycles
    cycle
    %if cycle == 4 % For 4 pulses
    if cycle == cycles % Recall the last set has mod = 0;
    %if M == 0 % This is NONSENSE
        I = index(find(M == 0));
        mod(index,cycle);
    %elseif cycle ~= 4
    %elseif M ~= 0 % This is NONSENSE
    elseif cycle ~= cycles
        I = index(find(M == (cycle)));
        mod(index,cycle);
    end;
   
    Z = Dzeroed(:,I);
    Zmean = mean(Z,2);
    %cycle
    pulses = pulsenums(cycle);
    traces = length(I);
    for p = 1:pulses;
        interval_start = start_point + (p-1)*ISI + 2*p*stim_duration;
        interval_end = start_point + (p * ISI) + 2*p*(p * ISI);
        %interval_end = start_point + 200;
        interval_start = 4900 + (p-1)*(ISI + stim_duration);
        if p ~=pulses
            interval_end = interval_start + 200;
        elseif p == pulses
            interval_end = length(d);
        end;
        
        %interval_end = 4900 + p*(ISI + stim_duration);
        
        IPSCs = Z(interval_start:interval_end,:); % Restrict data range to evoked IPSCs and traces of a given cycle
      
        pulsepeaks = min(IPSCs);
        
        %locations = find(IPSCs == pulsepeaks);
        ROI = interval_start:interval_end;
        ROIs= bsxfun(@times, ones(length(ROI), length(I)), ROI');
        
        l = bsxfun(@eq, IPSCs, pulsepeaks);
        locations = ROIs(l) + start_point;
        
        A{cycle,p} = pulsepeaks;
        
        L{cycle,p} = locations;
        
        
        
    end;
    
    Amp = [];
    D = [];

    for trace = 1:length(I)
        amp_span = Dzeroed(interval_start:interval_end,I(trace));
        min_amp = unique(min(amp_span));
        min_loc = unique(min(find(amp_span == min_amp)));
        amp_decay = Dzeroed(interval_start:end,I(trace));
        
        % Move along forward from the peak to find where the 95% value lays
        for n = min_loc:length(amp_span)
            if amp_decay(n) >= (0.95)*min_amp % Remember dealing with negative values
                first_pt = n
                disp('hello') 
                break;
            end;
        end
        
        % Move along forward from the peak to find where the 5% value lays
        for n = min_loc:length(amp_decay)
            if amp_decay(n) >= (0.05)*min_amp
                second_pt = n
                break;
            end;
        end 
        
        amp_decay = second_pt - first_pt
        D = [D amp_decay];
        Amp = [Amp min_amp];
    end
    
    Amps_Decays{cycle} = [Amp' D'];
    
meantrace = mean(Dzeroed(:,I),2);

meanIPSC = meantrace(interval_start:interval_end,:)

meanpeak = min(meanIPSC)

meanlocation = find(meantrace == meanpeak)

% Be sure to avoid multipe values of the mean amplitude!
unique_meanlocation=min(meanlocation); 

% Move along forward from the peak to find where the 95% value lays
for n = unique_meanlocation:length(meantrace)
    if meantrace(n) >= (0.95)*meanpeak % Remember dealing with negative values
        first_pt = n
        break;
    end;
end

% Move along forward from the peak to find where the 5% value lays
for n = unique_meanlocation:length(meantrace)
    if meantrace(n) >= (0.05)*meanpeak
       second_pt = n
        break;
    end;
end

% For subplot formating of axis (see below)
time_decay = second_pt - first_pt
    if time_decay > maxdecay
        maxdecay = time_decay;
    end
    %meanDecays(cycle) = time_decay;
% Plot the mean last ampltude of the given cycle
%figure, 

figure(p+20)
plot(meantrace), 

yl = ylim;
line([interval_start interval_start], [yl(1) yl(2)]); 

figure(100), subplot(cycles,2,2*(cycle-1)+1)
plot(meantrace), 

subplot(cycles,2,2*(cycle-1)+2)

offset = first_pt - interval_start; 
plot(meantrace(interval_start:second_pt + 500))
hold on
plot((offset + 1:offset + time_decay + 1), meantrace(first_pt:second_pt));
meantrace(first_pt) % Leave without semicolon to report mean 95% value.

% For subplot formatting of axis (see below)
if max(abs(ylim))  maxylim
    maxylim = max(abs(ylim));
end;

% Make also a copy of each subplot as its own figure
figure(200 + cycle)
offset = first_pt - interval_start; 
plot(meantrace(interval_start:second_pt + 500))
hold on
plot((offset + 1:offset + time_decay + 1), meantrace(first_pt:second_pt));


%decay_time = upperbound - decay_range(1);
decay_time = time_decay;
%n = length(decay_time)*(1/sampling_rate)*1000;
n = decay_time*(1/sampling_rate)*1000;
N = ['The decay time is', ' ', num2str(n), ' ms.'];
disp(N);

% Assign meanpeak values to a vector:
meanPeaks(cycle) = meanpeak;

% Assign decay values to a vector:
meanDecays(cycle) = n;
%meanAD(cycle) = [meanPeaks' meanDecays'];
end;

meanAD= [meanPeaks' meanDecays'];

for cycle = 1:cycles
    
    % Make also a copy of each subplot as its own figure
    figure(200 + cycle)
    xlim([0 (2*ISI + maxdecay)]);
    ylim([-maxylim 10]);
    
    % Neatly reformat x and y limits of the subplot
    figure(100)
    subplot(cycles,2,2*(cycle-1)+2)
    xlim([0 (2*ISI + maxdecay)]);
    ylim([-maxylim 10]);
end


X = 1:pulsenums(end)
Y = NaN(1,pulsenums(end))

Y(pulsenums(1:end)) = meanDecays(1:end) % More general case
figure, scatter(X,Y)

%% Save Data Files to a .MAT file
    
% IndividualAmplitudeData = struct('IndividualAmplitudes',A,'IndividualLocations',L);
% Individual_AmpDecay_Data = struct('IndividualDecays',Amps_Decays);
% Mean_AmpDecay_Data = struct('MeanAmplitudesDecayData', meanAD);
% 
% FileData = struct('IndividualAmplitudeData', IndividualAmplitudeData, ...
%     'Individual_AmpDecay_Data', Individual_AmpDecay_Data, ...
%     'Mean_AmpDecay_Data', Mean_AmpDecay_Data);
% name = strcat(FileName(1:end-4),'_Data','.mat');
% 
% save(name, 'FileData');

IndividualAmplitudeData = struct('A',A,'L',L);
IndividualAmpDecayData = struct('Amps_Decays',Amps_Decays);
MeanAmpDecayData = struct('meanAD', meanAD);

FileData = struct('index', index, 'numtraces', numtraces,'IndividualAmplitudeData', IndividualAmplitudeData, ...
    'Individual_AmpDecayData', IndividualAmpDecayData, ...
    'MeanAmpDecayData', MeanAmpDecayData);
name = strcat(FileName(1:end-4),'_Data','.mat');

%save(name, 'FileData');
% Show in the command window the meanAD data for copy and paste
format shortG
abs(meanAD(:,1)) % Mean absolute data
meanAD(:,2) % Mean decay data

%For the case of two pulses only for amplitudes and PPR:
%[mean(A{1,1}); mean(A{1,2}); mean(A{1,2})/mean(A{1,1})]
[-mean(A{1,1}); -mean(A{1,2}); mean(A{1,2})/mean(A{1,1});meanAD(:,2)]

copy2excel = {'1st Amplitude','2nd Amplitude','PPR', 'Decay Time'; -mean(A{1,1}), -mean(A{1,2}), mean(A{1,2})/mean(A{1,1}), meanAD(:,2)}
data = [-mean(A{1,1}), -mean(A{1,2}), mean(A{1,2})/mean(A{1,1}), meanAD(:,2)];
xlswrite(strcat(FileName,'_PPR-Data.xls'), data);
%%  ! ! ! HOW TO ACCESS SAVED DATA ! ! !
% 
% F = load('15811018_Data.mat')
% 
% F = 
% 
%    FileData: [1x1 struct]
%    F.FileData
%
%    ans = 
%
%                      index: [1x40 double]
%                  numtraces: 40
%    IndividualAmplitudeData: [4x10 struct]
%    Individual_AmpDecayData: [1x4 struct]
%           MeanAmpDecayData: [1x1 struct]
% 
% % E.G. TO ACCESS THE 2ND CYCLE AND 2 PULSES DATA: 
% F.FileData.IndividualAmplitudeData(2,2).A  % USES PARENTHESES () !!
%
% % TO ACCESS THE AMP_DECAY DATA:
% F.FileData.Individual_AmpDecayData(2).Amps_Decays
% 
% ans =
% 
%    1.0e+03 *
% 
%    -0.2696    4.0390
%    -0.2426    2.8290
%    -0.2069    2.5040
%    -0.1581    2.5100
%    -0.1642    3.0470
%    -0.2509    3.2350
%    -0.1475    3.7790
%    -0.1217    2.4820
%







