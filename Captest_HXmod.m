clc, clear all; close all;

%% Calculate capacitance
% Instructions: Put script in same folder as files to be analyzed. Set as
% your matlab directory.
% Written by Vivian Hernández version 3/20/2015. Email me if something
% doesn't make sense or if you find mistakes: vivianmhg@gmail.com
%
% Modified by HSX 3/23/2020
%

capstart=5.156e4;
basestart=5.9e4;
baseend=6.1e4;

files = dir('*.abf');
alldat=[];
tau1s=[];
tau2s=[];
weightedtaus=[];
fastpercent=[];
slowpercent=[];
ttrapz=[];
Cmall=[];
Rtotall=[];
deltaI_all=[];
ra_all=[];
goodnessrsquare=[];

fprintf('Type Group name using underscores, i.e. GS_Control_dSPN_Capacitance\n');

GroupName  = input('Type Group name here: ','s');

for i=1:length(files)
    
    file=files(i).name;
    [d,si]=abfload(file); 
    
    [a,b,c]=size(d);
    dn=reshape(d,[a,c,b]);
    d2=mean(dn,2);
    d2=d2*-1;                   %flip so easier to work with
    
    d2=interp1(1:length(d2),d2,1:.1:10000);  %interpolate so there are now 10 points for every original data point
    %sampling interval is 100 microseconds, sampleperms=100because interpolated, 
    %baseline=mode(d2(50:1050));
    
    basepulse=mode(d2(basestart:baseend));  % Baseline after cap trans decays
    d2=d2-basepulse;                        % Shift data relative to base
    
%     figure;plot(d2)
%     title(file);
    
    [k, ind]=max(d2(1:basestart)); %Find peak of cap trans
    
    ideal90=.9*k; 
    ideal10=.10*k; 
    
    [min90, ind90]=min(abs(d2(ind:ind+4500)-ideal90));
    [min10, ind10]=min(abs(d2(ind:ind+4500)-ideal10));

    
    % Fit to biexponential curve
    curveforfitting=d2(ind+ind90:ind+ind10);
    x=1:length(curveforfitting);

    fo=fitoptions('METHOD','NonlinearLeastSquares'); 
    [curve goodness] = fit(x',curveforfitting', 'exp2',fo );

    tau1=(-1/curve.b);
    tau2=(-1/curve.d);

    totamp=curve.a+curve.c;
    weightedtau=((curve.a/totamp)*tau1 +(curve.c/totamp)*tau2);

    xfits=1:length(d2(ind90:ind10));
    yfits=curve.a*(exp(curve.b*x)) + curve.c*(exp(curve.d*x));
  
% % Plot fit
%     figure;plot(d2(ind+ind90:ind+ind10),'k-'); hold on;
%     plot(xfits,yfits,'r-');

    tau1s=[tau1s, tau1];
    tau2s=[tau2s,tau2];
    goodnessrsquare=[goodnessrsquare, goodness.rsquare];
    weightedtaus=[weightedtaus,weightedtau];
    fastpercent=[fastpercent, (curve.a/totamp)*100];
    slowpercent=[slowpercent, (curve.c/totamp)*100];
    
    %Find area under curve
    diffcurrforAC=d2(capstart:basestart); %Calculate area under this portion of the curve
    diffcurrforAC(diffcurrforAC<0)=0;
     
    baseline=mode(d2(50:1050)); %relative to basepulse because basepulse is at 0
    deltaI=abs(baseline); 
    deltaI_all=[deltaI_all, deltaI];

    trp=trapz(diffcurrforAC)*(10); %in A*10-12 * units where units=(100*s*10^-3) ***need pCoulombs (A*s*1-^-12)
   
    %Find capacitance
    correctionfactor=1000*(weightedtau/100)*(deltaI); %correction = deltaI (pA) * tau (ms) = 10^-15, *1000 = 10^-12
    Cm=((trp+correctionfactor)/10)/1000; %C=area under curve(pC)/applied voltage(mV) (divide by 1000 to get to 10^-12)
    ttrapz=[ttrapz, trp];
    Cmall=[Cmall, Cm]; %area = mV * samples *
    
    %Total membrane resistance
    Rtot=(10/deltaI)*1000;%Rm=deltaV/I  - in megaohms (mV/pA)
    Rtotall=[Rtotall,Rtot]; %pA / mV = ohm*9

    %Access resistance
    ra=1000*(weightedtau/100)/(Cm); %(sE-3/ FE-12)= 10^9 ohm  %s/ohm=F
    ra_all=[ra_all,ra];
    

    
end
Rmall=Rtotall-ra_all; %Membrane resistance
fprintf('[Cmall, Rmall, ra_all,goodnessrsquare]');
alldat=[Cmall',Rmall',ra_all',goodnessrsquare'] %All data in one
%files.name
%xlswrite('RC_Control_dSPN.xls', alldat)
DataSetName = strcat(GroupName, '.xls');
%xlswrite('RC_Control_dSPN.xls', alldat)
xlswrite(DataSetName, alldat)


%%
fprintf('File names (in order):\n');

for i=1:length(files)
    
    n=files(i).name;
    fprintf('%s\n', n(1:8))
end

