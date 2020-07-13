%%
clear all; %close all; clc

current_dir = pwd;
[FileName,PathName] = uigetfile('*.txt','SELECT THE .txt FILE');
cd(current_dir);

% %%
% %fileID = fopen('nums1.txt','r');
% fileID = fopen(FileName,'r');
% %A = fscanf(fileID,formatSpec)
% %A = fscanf(fileID, '%s', [200, 200])
% A = fscanf(fileID, 'w')
% fclose(fileID);
% 
% % infile = 'hawaii_coast.dat.txt';
% % data = load(infile, '-ascii');

%%
matrix = load(FileName, '-ascii');
%matrix = [1 1; 1 -1; -1 1; -1 -1; 1 0; -1 0; 0 1; 0 -1];
x = matrix(:,1); % column 1 of the data text file is assigned the variable x
y = matrix(:,2); % column 2 is assigned the variable y

%%
% regionprops(true(size(matrix)), matrix,  'WeightedCentroid')

%%
% matrix=matrix/sum(matrix(:));
% [m,n]=size(matrix);
% [I,J]=ndgrid(1:m,1:n);
%    centroid=[dot(I(:),matrix(:)),  dot(J(:),matrix(:))]

%%
xcentroid = sum(x)/sum(ones(1,length(x)))
ycentroid = sum(y)/sum(ones(1,length(y)))

%%
P =nan(length(matrix),2);

data_pt = 1;

while (data_pt <= length(matrix))
    X = (x(data_pt) - xcentroid);
    Y = (y(data_pt) - ycentroid);
    r = sqrt(X^2 + Y^2);
%     %theta = atan ((y(data_pt) - ycentroid)/(x(data_pt) - xcentroid));
%     if Y >= 0 && X >= 0
%         theta = acos(X/r)
%     elseif Y >= 0 && X < 0
%         theta = acos(X/r) + pi/2
%     elseif Y < 0 && X < 0
%         theta = acos(X/r) + pi
%     elseif Y < 0 && X >= 0
%         theta = acos(X/r) + 2*pi/3
%     end

% %theta = atan ((y(data_pt) - ycentroid)/(x(data_pt) - xcentroid));
%     if Y >= 0
%         theta = acos(X/r)
%     elseif Y < 0
%         theta = acos(X/r) + pi
%         
% %     elseif Y < 0 && X < 0
% %         theta = acos(X/r) + pi
% %     elseif Y < 0 && X >= 0
% %         theta = acos(X/r) + 3*pi/2
%     end

    if Y >= 0 
        theta = acos(X/r);
    elseif Y < 0

        if X <= 0
            %X = -X;
            theta = (acos(-X/r) + pi);
        elseif X > 0
            theta = ((pi/2 -acos(X/r)) + 3*pi/2);
            %theta = 0;
        end   
    end

    P(data_pt,:) = [r theta];
    data_pt = data_pt + 1;
    
end

%polarhistogram(theta,nbins)
figure
%nbins  = 25;
nbins = ceil(log10(length(x))*10);

%h = polarhistogram(P(:,2),nbins);
h = polarhistogram(P(:,2),'BinEdges', linspace(0,2*pi,nbins));

response = [num2str(nbins),' ', 'bins were used'];
display(response);

%%
figure, scatter(x,y, 'filled', 'k');
axis([0 1 0 1]);
hold on
scatter(xcentroid, ycentroid, 'filled', 'r');
figure, polarplot(P(:,2), P(:,1));

xlswrite(strcat(FileName(1:end-4),'_PolarVariables.xls'), P)


%%
display('BinEdges');
h.BinEdges';
display('BinValues');
h.Values';

%%
c = ksdensity([x,y], [x,y]);
figure, h = scatter(x, y, 50,c,'filled')
colormap(viridis), colorbar
hold on, scatter(xcentroid, ycentroid, 55, [128,128,128]/256, 'filled')
axis([0 1 0 1]);
title('Viridis');

%%
% [rho, theta] = cart2pol(x-xcentroid, y-ycentroid);
% [t,r] = meshgrid(theta,rho);
% [xx,yy] = pol2cart(t,r);
% P = peaks(xx,yy);
% %P = peaks(c,c);
% figure('color','white');
% polarplot3d(P);
% view([-18 76]);
% 

%% End here to prevent the random generator portion to run
% The below code was to randomly sample from the combined
% HuCD data and plot the polarhistograms 

return

%%
%The below code was to randomly sample from the combined
% HuCD data and plot the polarhistograms 
%FileName = 'Combined_HuCD_CartesianCoordinates.txt';

clear all; close all; clc

current_dir = pwd;
[FileName,PathName] = uigetfile('*.txt','SELECT THE .txt FILE');
cd(current_dir);

D = load(FileName, '-ascii');
%matrix = [1 1; 1 -1; -1 1; -1 -1; 1 0; -1 0; 0 1; 0 -1];

% FileName = '00_Npas1combined19305.txt';
% D = load(FileName, '-ascii');

data_length = length(D);
sample_size = ceil(0.5*data_length) % Select 50% of the data to sample
indx_sample = randperm(data_length,sample_size);

matrix = D(indx_sample,:);

x = matrix(:,1); % column 1 of the data text file is assigned the variable x
y = matrix(:,2); % column 2 is assigned the variable y

xcentroid = sum(x)/sum(ones(1,length(x)))
ycentroid = sum(y)/sum(ones(1,length(y)))

% Go then to section line #38 and run down the sections
