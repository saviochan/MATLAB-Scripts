clear all; close all;

% file='16n16012.abf';
% % Enter in usable sweep
% %cnew=21;
% 
% % Load the file
% [d,si]=abfload(file);

current_dir = pwd;
[FileName,PathName] = uigetfile('*.abf','SELECT THE .abf FILE');
cd(current_dir);
[d,si,h] = abfload(FileName);

% % Change here to set range of traces.
% D=d;
% d=D(:,:,1:40);

%%
% Reshape the file into a matrix
[a,b,c]=size(d);
dnew=reshape(d,a,c);

dnewstart=1;
%c=22; % Use this to limit to a range of traces
dnewend=c;
% Plotting a figure
figure; plot(dnew(:,dnewstart:dnewend));
figure, plot(dnew(:,c));

% Defining time points to detect spikes in a 500ms current step
startpulse=1159;
endpulse=6000;

% Find baseline of each sweep
baseline=mean(dnew(1:startpulse,:));
baselinegood=baseline(dnewstart:dnewend);
RMP=mean(baselinegood);

% Look at sweepmodes of certain sweeps
sweepmode=mean(dnew(startpulse:endpulse,:));
%figure; plot(sweepmode);

% Making a NaN matrix 
allpeaks=NaN(100,dnewend);
allpeaklocs=NaN(100,dnewend);
alltroughs=NaN(100,dnewend);
alltlocs=NaN(100,dnewend);

% Putting peak values into the NaN matrix 
for j=dnewstart:dnewend;
    [peaks,locs]=findpeaks(dnew(startpulse:endpulse,j),'minpeakheight',sweepmode(j)+15,'minpeakdistance',50,'minpeakprominence',20);
    if peaks;
        allpeaks(1:length(peaks),j)=peaks;
        allpeaklocs(1:length(locs),j)=locs;
    end
end

% Plot peaks found by findpeaks
p=c;
figure;
plot(dnew(:,p)); hold on;
plot(allpeaklocs(:,p)+startpulse-1,allpeaks(:,p),'r*');
hold off;

% Action potentials per sweep
numpeaks=isfinite(allpeaks);
sumnumpeaks=sum(numpeaks);
sumnumpeakstrans=transpose(sumnumpeaks);
[maxpeaksvalue,maxpeakssweep]=max(sumnumpeaks);
maxpeakscurrent=-150+25*(maxpeakssweep-1);

% Interspike interval
ISI=diff(allpeaklocs);
ISImineach=min(ISI);
ISImineachtrans=transpose(ISImineach);
ISImaxeach=max(ISI);
ISImaxeaechtrans=transpose(ISImaxeach);
ISIminall=min(ISImineach)/10;

% Rheobase
arrayfirstAP=find(sumnumpeaks>0);
firstAPsweep=arrayfirstAP(1);
rheobase=-150+25*(firstAPsweep-1);

% Find threshold and spike height for first peak
firstAPloc=allpeaklocs(1,firstAPsweep);
firstAPlatency=firstAPloc/10;
firstAP=allpeaks(1,firstAPsweep);
intfirstAP=interp1(firstAPloc+startpulse-500:firstAPloc+startpulse-1,dnew(firstAPloc+startpulse-500:firstAPloc+startpulse-1,firstAPsweep),firstAPloc+startpulse-500:0.25:firstAPloc+startpulse-1,'spline');

slope=diff(intfirstAP);
[maxslopey,maxslopex]=max(slope);
maxdvdt=10*maxslopey;
slopechange=diff(slope);

thresholdarray=find(slopechange>0.1);
threshold=thresholdarray(1);
thresholdvalue=intfirstAP(threshold);

figure; plot(dnew(:,firstAPsweep)); hold on;
plot(firstAPloc+startpulse-500:firstAPloc+startpulse-1,dnew(firstAPloc+startpulse-500:firstAPloc+startpulse-1,firstAPsweep),'o',firstAPloc+startpulse-500:0.25:firstAPloc+startpulse-1,intfirstAP,':.');

%plot(firstAPloc+startpulse-500:0.25:firstAPloc+startpulse-1,intfirstAP,':.'); hold on;
plot(threshold/4+firstAPloc+startpulse-500.25,thresholdvalue,'c*');

spikeheight=firstAP-thresholdvalue;

% Find fAHP of first spike
[fAHPpeaktroughy,fAHPpeaktroughx]=min(dnew(firstAPloc+startpulse-1:firstAPloc+startpulse+100,firstAPsweep));
plot(fAHPpeaktroughx+firstAPloc+startpulse-2,fAHPpeaktroughy,'r*');
fAHPfirstpeak=abs(fAHPpeaktroughy-thresholdvalue);

% Find resistance of the cell
baseline=dnew(1250,1:(firstAPsweep-1));
baselinetrans=transpose(baseline);
figure; plot(dnew(:,1:(firstAPsweep-1))); hold on;
plot(1250,baseline,'r*'); hold off;

currentinjection=-150:25:(rheobase-25);

figure; plot(currentinjection,baseline,'mo');
line=lsline;
p=polyfit(get(line,'xdata'),get(line,'ydata'),1);
resistance=p(1)*1000;


    % Find threshold and spike height for ten sweeps
    tensweeps=find(sumnumpeaks>=9 & sumnumpeaks<=11);
    if tensweeps
        tensweepsnumberAP=sumnumpeaks(tensweeps);
        tensweep=tensweeps(1);
        tensweepnumberAP=sumnumpeaks(tensweep);
        %inttensweep=interp1(1:10000,dnew(1:10000,tensweep),1:0.25:10000,'spline');
    
        tensweepsweep=dnew(:,tensweep);

        figure; plot(tensweepsweep); hold on;
        %plot(1:0.25:10000,inttensweep,'o');
        %inttensweepslope=diff(inttensweep);
        %plot(1:0.25:9999,inttensweepslope(1:end-3),':.');
        plot(allpeaklocs(:,tensweep)+startpulse-1,allpeaks(:,tensweep),'r*');
    
        inttensweepslopechange=diff(dnew(:,tensweep));
        plot(inttensweepslopechange);
    
        thresholdsweeparray=NaN(tensweepnumberAP,1);
        for k=1:tensweepnumberAP;
            m=allpeaklocs(k,tensweep)+startpulse-1;

            thresholdsweeparraytemp=find(inttensweepslopechange(m-50:m)>0.10);
            thresholdsweeparray(k)=thresholdsweeparraytemp(1)+m-49;
        end

        thresholdsweeparrayvalue=tensweepsweep(thresholdsweeparray);
        plot(thresholdsweeparray,thresholdsweeparrayvalue,'c*');

        tensweeppeaks=allpeaks(:,tensweep);
        tensweeppeaks=tensweeppeaks(isfinite(tensweeppeaks));
        tensweepsspikeheight=tensweeppeaks-thresholdsweeparrayvalue;

        tensweepISI=ISI(:,tensweep);

        % Find fAHP of ten spike sweep
        fAHPpeaktroughy=NaN(tensweepnumberAP,1);
        fAHPpeaktroughx=NaN(tensweepnumberAP,1);
        for n=1:tensweepnumberAP;
            q=allpeaklocs(n,tensweep)+startpulse-1;
            [fAHPpeaktroughy(n),fAHPpeaktroughx(n)]=min(dnew(q:q+100,tensweep));
            fAHPpeaktroughx(n)=fAHPpeaktroughx(n)+q-1;
        end

        plot(fAHPpeaktroughx,fAHPpeaktroughy,'r*');

        fAHPtensweep=abs(thresholdsweeparrayvalue-fAHPpeaktroughy);
    end
    
       
% Width at Half Max of First Peak
[peak1,loc1]=findpeaks(dnew(startpulse:endpulse,firstAPsweep),'minpeakdistance',50,'minpeakheight',sweepmode(firstAPsweep)+15,'minpeakprominence',20);
dnewshift=dnew(:,firstAPsweep)+abs(thresholdvalue);
%plot(dnewshift);

[peak,loc,width,prom]=findpeaks(dnewshift(startpulse:endpulse),'minpeakdistance',50,'minpeakprominence',20);
figure; findpeaks(dnewshift(startpulse:endpulse),'minpeakdistance',50,'minpeakprominence',20,'Annotate','extents','WidthReference','halfheight');

FWHM=width(1)/10;

% Properties of -150mV Current Step
% Find RMP of 1st sweep
rmpmed1=mean(dnew(1:1100,1));


% Trough of sweeps 1
[troughy,troughx]=min(dnew(startpulse:2500,1));
troughvalue=troughy;
trough=rmpmed1-troughvalue;
% Steady state of sweep 1
steadystatevalue=mean(dnew(5000:5600,1));
steadystate=rmpmed1-steadystatevalue;
% Sag
sag=steadystatevalue-troughvalue;
sagratio=steadystate/trough;
% Plot the points found
figure;
plot(dnew(:,1)); hold on; 
plot(1000,rmpmed1,'k*');
plot(troughx+startpulse-1,troughy,'r*');
plot(5800,steadystatevalue,'g*'); hold off;

% Copy to excel
if tensweeps;
    copy2excel=[rmpmed1,trough,troughvalue,steadystate,steadystatevalue,sag,sagratio,RMP,resistance,rheobase,thresholdvalue,spikeheight,maxdvdt,FWHM,fAHPfirstpeak,maxpeaksvalue,maxpeakscurrent,ISIminall,firstAPlatency,tensweep,tensweepnumberAP];
    else
    copy2excel=[rmpmed1,trough,troughvalue,steadystate,steadystatevalue,sag,sagratio,RMP,resistance,rheobase,thresholdvalue,spikeheight,maxdvdt,FWHM,fAHPfirstpeak,maxpeaksvalue,maxpeakscurrent,ISIminall,firstAPlatency];
end

%% Harry's added code: display Spiking values and other parameters

FileName(1:8) % Display Filename
sumnumpeaks'
output = [rmpmed1,trough,troughvalue,steadystate,steadystatevalue,sag,sagratio,RMP,resistance,rheobase,thresholdvalue,spikeheight,maxdvdt,FWHM,fAHPfirstpeak,maxpeaksvalue,maxpeakscurrent,ISIminall,firstAPlatency,tensweep,tensweepnumberAP]'

%% Save data output
xlswrite(strcat(FileName,'_sumnumpeaks.xls'), sumnumpeaks')
xlswrite(strcat(FileName,'_ephys-parameters.xls'), output)
