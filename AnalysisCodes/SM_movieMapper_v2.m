function varargout = SM_movieMapper_v2(varargin)
% SM_MOVIEMAPPER_V2 MATLAB code for SM_movieMapper_v2.fig
%      SM_MOVIEMAPPER_V2, by itself, creates a new SM_MOVIEMAPPER_V2 or raises the existing
%      singleton*.
%
%      H = SM_MOVIEMAPPER_V2 returns the handle to a new SM_MOVIEMAPPER_V2 or the handle to
%      the existing singleton*.
%
%      SM_MOVIEMAPPER_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_MOVIEMAPPER_V2.M with the given input arguments.
%
%      SM_MOVIEMAPPER_V2('Property','Value',...) creates a new SM_MOVIEMAPPER_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_movieMapper_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_movieMapper_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_movieMapper_v2

% Last Modified by GUIDE v2.5 27-May-2022 10:28:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SM_movieMapper_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @SM_movieMapper_v2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SM_movieMapper_v2 is made visible.
function SM_movieMapper_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_movieMapper_v2 (see VARARGIN)

% Choose default command line output for SM_movieMapper_v2
handles.output = hObject;
handles.ch1peaks=[];
handles.ch2peaks=[];
handles.M21x=[]; %transformation matrixes, explicitly, no matlab canned transformation routines were used
handles.M21y=[];
handles.M12x=[];
handles.M12y=[];

handles.PointCor=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_movieMapper_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_movieMapper_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)

%% 1. Allows User to Select Tif file---------------
DefaultFilePath='Z:\Valter Z\September 26, 2017\';
[fileName,pathName]=uigetfile('*.tif','Select ".tif" file',DefaultFilePath);
%% 1.-----------------------------------------------.

%% 2. gather Info about Tif image-----------.
FileTif=strcat(pathName,fileName);
handles.FileTif=FileTif;
InfoImage=imfinfo(FileTif);
handles.mImage=InfoImage(1).Width;
handles.nImage=InfoImage(1).Height;
handles.NumberImages=length(InfoImage);
%% 2.-------------------------------------.
%% 3. Reads and averages tif stack images and gives an average image.
firstFrame = imread(FileTif, 'index', 1);
Isum = uint16(zeros(size(firstFrame,1),size(firstFrame,2))); %create space for the sum image.
for frame = 1:5%handles.NumberImages %%%%%% Change later - so far is averaging the first five frames
    AvgImg = imread(FileTif, 'index', frame);
    Isum = Isum+AvgImg;
end
%AvgImg = Isum./handles.NumberImages;
AvgImg = Isum./5;
%% 3.-------------------------------------.
%% 4. split Image for different channel
howManyColors=2;
if howManyColors==2
    CH1AvgImg=AvgImg(1:handles.nImage,1:handles.mImage/2);
    CH2AvgImg=AvgImg(1:handles.nImage,handles.mImage/2+1: handles.mImage);
    %% 4.-------------------------------------.
    %% 5. Find background for different channel - do not like this disk business...
    backgroundC1 = imopen(CH1AvgImg,strel('disk',35));
    backgroundC2 = imopen(CH2AvgImg,strel('disk',35));
    %% 5.-------------------------------------.
    %% 6. make BG corrected image for different channel
    %CH1AvgImg=CH1AvgImg-backgroundC1;
   % CH2AvgImg=CH2AvgImg-backgroundC2;
    %% 6.-------------------------------------.
    %% 7. plot image for different channel
    axes(handles.beadImage_CH1);
    imshow(imadjust(CH1AvgImg)) %show image CH1
    axes(handles.beadImage_CH2);
    imshow(imadjust(CH2AvgImg)) %show image CH2
    %% 7.-------------------------------------.
end
%% 8. Initialize all different parameters.
handles.PointCor=uint16([size(CH1AvgImg,1)/2,size(CH1AvgImg,2)/2,size(CH2AvgImg,1)/2,size(CH2AvgImg,2)/2;
                     size(CH1AvgImg,1)/2,size(CH1AvgImg,2)/2,size(CH2AvgImg,1)/2,size(CH2AvgImg,2)/2;
                     size(CH1AvgImg,1)/2,size(CH1AvgImg,2)/2,size(CH2AvgImg,1)/2,size(CH2AvgImg,2)/2]);

handles.CH1Image(:,:,1)=CH1AvgImg;
handles.CH2Image(:,:,1)=CH2AvgImg;
handles.CH1Image(:,:,2)=uint16(zeros(size(CH1AvgImg))); %why is it keeping all those zeros...
handles.CH2Image(:,:,2)=uint16(zeros(size(CH2AvgImg)));
handles.CH1Image(:,:,3)=uint16(zeros(size(CH1AvgImg)));
handles.CH2Image(:,:,3)=uint16(zeros(size(CH2AvgImg)));
handles.CH1AvgImgDup=CH1AvgImg;%keep duplicate images for reset.
handles.CH2AvgImgDup=CH2AvgImg;%keep duplicate images for reset

axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
%synchronize the individual axis limits across several figures or subplots within a figure. 
zoom on;
%% 8.-------------------------------------.

handles.PointCor=[];
guidata(hObject, handles); %set handles to access from other functions.   

%% 9. Update Status:
statement=['Movie: ' fileName ', is loaded successfully.'];
set(handles.status,'String',statement);
drawnow;pause(0.1);
guidata(hObject, handles);%set handles to access from other functions.  

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
handles.CH1Image(:,:,1)=handles.CH1AvgImgDup;
handles.CH2Image(:,:,1)=handles.CH2AvgImgDup;
CH1AvgImg=handles.CH1AvgImgDup;
CH2AvgImg=handles.CH2AvgImgDup;
handles.CH1Image(:,:,2)=uint16(zeros(size(CH1AvgImg)));
handles.CH2Image(:,:,2)=uint16(zeros(size(CH2AvgImg)));
handles.CH1Image(:,:,3)=uint16(zeros(size(CH1AvgImg)));
handles.CH2Image(:,:,3)=uint16(zeros(size(CH2AvgImg)));
axes(handles.beadImage_CH1);
imshow(imadjust(CH1AvgImg)) %show image CH1
axes(handles.beadImage_CH2);
imshow(imadjust(CH2AvgImg)) %show image CH2

handles.PointCor=[];

%% 9. Update Status:
statement=['Reset everything.'];
set(handles.status,'String',statement);
drawnow;pause(0.1);

guidata(hObject, handles);%set handles to access from other functions.


% --- Executes on button press in peakfinder.
function peakfinder_Callback(hObject, eventdata, handles)

%% 1. Gather available Information.
radius=str2double(get(handles.peakRadiusCh1,'String'));
thr=str2double(get(handles.threshCH1,'String'));    


%% Find and Plot CH1 peaks
CH1AvgImg=handles.CH1Image(:,:,1);
CH1AvgImg=uint16(CH1AvgImg);
%% 2. estimate Threshold 
Int=reshape(CH1AvgImg,[size(CH1AvgImg,1)*size(CH1AvgImg,2) 1]);
n=100;
threshold=mean(Int);
while n>0.1
    foreG=Int(Int>threshold);
    backG=Int(Int<threshold);
    thfinal=(mean(foreG)+mean(backG))/2;
    n=abs(thfinal-threshold);
    threshold =thfinal;
end
threshold=threshold*thr;
p=FastPeakFind(CH1AvgImg,threshold/4,(fspecial('gaussian', 7,1)), radius+1);
chpeaks=reshape(p',2,numel(p)/2);
chpeaks=chpeaks';
chpeaks=[chpeaks(:,2) chpeaks(:,1)]; %reshaping of the peak list...
I=diag(CH1AvgImg(chpeaks(:,2),chpeaks(:,1)));
TF=ones(size(chpeaks,1),1);
for i=1:size(chpeaks,1)
if I(i)>threshold*1.5
    TF(i)=0;
end
end
chpeaks(logical(not(TF)),:) = [];
Distcap=radius+1;
TF=zeros(size(chpeaks,1),1);
for i=1:size(chpeaks,1)
    C=chpeaks(i,:);
    D = pdist2(chpeaks,C);
    D(i,:)=[];
    [val,~]= min(D);
    if val<Distcap
        TF(i)=1;
    end
end
chpeaks(logical(TF),:) = [];
%% 5. Plot  CH1 all.....
axes(handles.beadImage_CH1);
handles.ch1peaks=chpeaks;
imshow(imadjust(CH1AvgImg));
axis off
hold on
scatter(chpeaks(:,2), chpeaks(:,1),'MarkerEdgeColor', 'g')
hold off


%% Find and Plot CH2 peaks
CH2AvgImg=handles.CH2Image(:,:,1);
CH2AvgImg=uint16(CH2AvgImg);
% 2. estimate Threshold 
Int=reshape(CH2AvgImg,[size(CH2AvgImg,1)*size(CH2AvgImg,2) 1]);
n=100;
threshold=mean(Int);
while n>0.1
    foreG=Int(Int>threshold);
    backG=Int(Int<threshold);
    thfinal=(mean(foreG)+mean(backG))/2;
    n=abs(thfinal-threshold);
    threshold =thfinal;
end
threshold=threshold*thr;
clearvars p chpeaks
p=FastPeakFind(CH2AvgImg,threshold/4,(fspecial('gaussian', 7,1)), radius+1);
chpeaks=reshape(p',2,numel(p)/2);
chpeaks=chpeaks';
chpeaks=[chpeaks(:,2) chpeaks(:,1)]; %reshaping of the peak list...
I=diag(CH2AvgImg(chpeaks(:,2),chpeaks(:,1)));
TF=ones(size(chpeaks,1),1);
for i=1:size(chpeaks,1)
if I(i)>threshold*1.5
    TF(i)=0;
end
end
chpeaks(logical(not(TF)),:) = [];
Distcap=radius+1;
TF=zeros(size(chpeaks,1),1);
for i=1:size(chpeaks,1)
    C=chpeaks(i,:);
    D = pdist2(chpeaks,C);
    D(i,:)=[];
    [val,~]= min(D);
    if val<Distcap
        TF(i)=1;
    end
end
chpeaks(logical(TF),:) = [];
%% 5. Plot  CH1 all.....
axes(handles.beadImage_CH2);
handles.ch2peaks=chpeaks;
imshow(imadjust(CH2AvgImg));
axis off
hold on
scatter(chpeaks(:,2), chpeaks(:,1),'MarkerEdgeColor', 'r')
hold off
zoom on;
statement= strcat('We found, ', num2str(size(handles.ch1peaks,1)), ' peaks in CH1, and  ' , num2str(size(handles.ch2peaks,1)), ' peaks in CH2.');
set(handles.status,'String',statement);
drawnow;pause(0.1);
guidata(hObject,handles);%set handles to access from other functions.   


% --- Executes on button press in firstPoint.
function firstPoint_Callback(hObject, eventdata, handles)
statement= 'Select the first point by clicking on a point on CH1 and CH2';
set(handles.status,'String',statement);
drawnow;pause(0.1);
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom off;
whichPointAreYouMoving (hObject,1,handles,1)
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom on;

% --- Executes on button press in secondPoint.
function secondPoint_Callback(hObject, eventdata, handles)
statement= 'Select the second point by clicking on a point on CH1 and CH2';
set(handles.status,'String',statement);
drawnow;pause(0.1);
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom off;
whichPointAreYouMoving (hObject,1,handles,2)
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom on;

% --- Executes on button press in thirdPoint.
function thirdPoint_Callback(hObject, eventdata, handles)
statement= 'Select the third point by clicking on a point on CH1 and CH2';
set(handles.status,'String',statement);
drawnow;pause(0.1);
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom off;
whichPointAreYouMoving (hObject,1,handles,3)
axes(handles.beadImage_CH1);
axes(handles.beadImage_CH2);
linkaxes([handles.beadImage_CH1, handles.beadImage_CH2])
zoom on;


function whichPointAreYouMoving (hObject,~,handles,whichPoint)
%% 1. Gather available Information.
PointCor=handles.PointCor;
%% 1.----------------------------------

%% 2. Find which point and what is the coordinate of the clicked point.
if whichPoint==1
    set(handles.secondPoint,'value',0);
    set(handles.thirdPoint,'value',0);
    A=get(handles.firstPoint,'value');
elseif whichPoint==2
    set(handles.firstPoint,'value',0);
    set(handles.thirdPoint,'value',0);
    A=get(handles.secondPoint,'value');
elseif whichPoint==3
    set(handles.secondPoint,'value',0);
    set(handles.firstPoint,'value',0);
    A=get(handles.thirdPoint,'value');
end

if A==1 %if not in error
    set(handles.beadImage_CH1);
    [y1,x1]=ginput(1);
    pause();
    set(handles.beadImage_CH2);
    [y2,x2]=ginput(1);
    xCH1=round(x1);
    xCH2=round(x2);
    yCH1=round(y1);
    yCH2=round(y2);
    
    PointCorThis=[xCH1,yCH1,xCH2,yCH2];
    PointCor(whichPoint,:)=PointCorThis;
    handles.PointCor=PointCor;
    guidata(hObject, handles);
    
    CH1AvgImg=handles.CH1Image(:,:,1);
    CH1AvgImg=uint16(CH1AvgImg);
    CH2AvgImg=handles.CH2Image(:,:,1);
    CH2AvgImg=uint16(CH2AvgImg);
    
    %% 3. Re-plot Image for the full version:
    %% 3a. For CH1
    axes(handles.beadImage_CH1);
    imshow(imadjust(CH1AvgImg));
    axis off;
    hold on
    scatter(handles.ch1peaks(:,2), handles.ch1peaks(:,1),'MarkerEdgeColor', 'g');
    hold on
    for i=1:size(PointCor,1)
    scatter(PointCor(i,2), PointCor(i,1),'MarkerEdgeColor', 'y', 'LineWidth',2);
    end
    hold off
    
    %% 3a. -------------------------------
    %% 3b. For CH2
    axes(handles.beadImage_CH2);
    imshow(imadjust(CH2AvgImg));
    axis off;
    hold on
    scatter(handles.ch2peaks(:,2), handles.ch2peaks(:,1),'MarkerEdgeColor', 'r')
    hold on
    for i=1:size(PointCor,1)
    scatter(PointCor(i,4), PointCor(i,3),'MarkerEdgeColor', 'y', 'LineWidth',2);
    end
    hold off
    %% 3b. -------------------------------
    %% 3. -------------------------------
end





% --- Executes on button press in map.
function map_Callback(hObject, eventdata, handles)
%% 1. Gather available Information ----------------------------------------
PointCor=handles.PointCor;
PointCor=double(PointCor);
radius=str2double(get(handles.peakRadiusCh1,'String'));
% Gather Info of the images.
CH1AvgImg=handles.CH1Image(:,:,1);
CH1AvgImg=uint16(CH1AvgImg);
CH2AvgImg=handles.CH2Image(:,:,1);
CH2AvgImg=uint16(CH2AvgImg);
% All peaks in channel 1 and 2
ch1peaks=handles.ch1peaks;
ch2peaks=handles.ch2peaks;
%% 1.----------------------------------------------------------------------
%% 2. Create the initial Transformation with three selected point pairs---.
% Create matrix :[1 1 1; x1 x2 x3; y1 y2 y3]
ICh1=[1 1 1;PointCor(:,1)';PointCor(:,2)']; %initial channel 1...
ICh2=[1 1 1;PointCor(:,3)';PointCor(:,4)'];
% ICh2= T12* ICh1
T12=ICh2* inv(ICh1); %transformation matrix

% Now use this transformation to find all CH2 points, given CH1 points.
% create matrix with many Ch1 peaks:[1 1 1....; x1 x2 x3....; y1 y2 y3....]
ACh1=[ones(1,size(ch1peaks,1));ch1peaks(:,1)';ch1peaks(:,2)']; %advanced
% use this matrix to find peaks in Ch2: we call it predicted CH2...
pCh2=T12*ACh1; %predicted channel 2 points equivalent to channel 1 points

%NB : so far all of them
pCh2peak=[pCh2(2,:)' pCh2(3,:)']; %non-integers, predicted
% Checking this output...
% figure(2)
% imshow(imadjust(CH2AvgImg))
% hold on;
% plot(ch2peaks(:,2),ch2peaks(:,1),'ro')
% hold on;
% plot(pCh2peak(:,2),pCh2peak(:,1),'go')
%% 2. ---------------------------------------------------------------------
%% 3 This part is probably unnecessary...Makes it easier down the road.----
ch2peaks=sortrows(ch2peaks,2); %sort rows based on column 2
pCh2peak=sortrows(pCh2peak,2); %sorted locations of peaks in channel 2, based on locations of channel 1 peaks
%% 3. ---------------------------------------------------------------------

%% 4: Find a set of points that can be used for the second round of fitting--
% For this I first tried the strategy to find all possible predicted points
% that are within few(1 or 2) pixel away from the actual point and then use
% them as seed for the next round. But this did not work very well. I think
% the reason is, it was only sampleing points from the good central part of
% the image. That doesnot represent the whole image and hence create very
% bad mapping.  This part is commented down below.
%% :::::::::::::::::::::::::::::::::::::::::::::::: %%
%% Instead I tried a different strategy.
% devide the image in 16 blocks and then choose two points from each block
% for the next round. This provides a more broad range of sampeling.
m=0;
InfoT=[];
% i and j represents the image is devided in 4x4 16 sections.
for i=1:4
    for j=1:4
        m=m+1;
        k{m}=ch2peaks(ch2peaks(:,1)>size(CH2AvgImg,1)*(j-1)/4 & ch2peaks(:,1)<size(CH2AvgImg,1)*j/4 & ch2peaks(:,2)>size(CH2AvgImg,2)*(i-1)/4 & ch2peaks(:,2)<size(CH2AvgImg,2)*i/4,:);
        pk{m}=pCh2peak(pCh2peak(:,1)>size(CH2AvgImg,1)*(j-1)/4 & pCh2peak(:,1)<size(CH2AvgImg,1)*j/4 & pCh2peak(:,2)>size(CH2AvgImg,2)*(i-1)/4 & pCh2peak(:,2)<size(CH2AvgImg,2)*i/4,:);
        Info=[];
        c=pk{m};  %predicted peak location(s) array
        d=k{m}; %original peak locations array
        if size(c,1)>=2 %only works if there are at least 2 points for 1/16th of an image
            for g=1:length(c)
                C=c(g,:); %C is one predicted point
                D = pdist2(d,C);
                %pdist2(X,Y) returns a matrix D containing the Euclidean distances between each pair of observations
                [val,ind]= min(D); %finds value and internal # (among the points in 1/16 of the point that
                %is the closest to the predicted point
                %it is not coordinates of any particular point
                Info=[Info; val ind d(ind,:) C];
                %structure of Info distance index, actual positions
                %(integer), predicted positions (integear)
                %Info contains some coordinates. Are they original or
                %predicted? I guess the first two coordinates are original,
                %and last two (non-integer) are predicted.
            end
            %after the end of this cycle it should have found the point
            %that is closest to prediction
            a=sortrows(Info,1);
            b=a(1:2,:);
            InfoT=[InfoT;b];
        end
    end
end
% InfoT has all the info for the images...
% InfoT=[(valof minimum diff) (ind of that val) (corresponding CH2 actual peak x) (corresponding CH2 actual peak y)
%      (corresponding CH2 predicted peak x) (corresponding CH2 predicted peak y)

%% Checking
% figure(2)
% imshow(imadjust(CH2AvgImg))
% hold on;
% plot(InfoT(:,4),InfoT(:,3),'ro')
% hold on;
% plot(InfoT(:,6),InfoT(:,5),'go')
%% 4. ---------------------------------------------------------------------


%% 5. Make a list of the selected points.
% format InfoT properly.
F_pC2=[ones(1,size(InfoT,1));InfoT(:,5)';InfoT(:,6)']; % formated predicted CH2 peaks (columns 5 and 6 of Info)

F_AC1=inv(T12)*F_pC2; % Found corresponding formated actual CH1 peaks

%OK, so above we transformed channel 1 peaks to produce "predicted channel
%2 peaks" but the results here look integer, so inverse transform must be
%working well

% matching CH1 and CH2 pairs. All of these are actual points (not predicted)
%why, given that we just calculated F_AC1 from predicted CH2 peaks?
MCH1CH2pairs=[F_AC1(2,:)' F_AC1(3,:)' InfoT(:,3) InfoT(:,4)]; % 3, 4 original channel 2 points
%would be safer if we had here original channel one points in 1 and 2, but it is
%admittedly difficult to find channel one points who are the counterparts
%of selected channel 2 points.

% Thus, the first two items are channel 1 points converted from channel 2
% and
Xc1=MCH1CH2pairs(:,2); %twice transformed channel 1 points
Yc1=MCH1CH2pairs(:,1); %twice transformed channel 1 points
xc2=MCH1CH2pairs(:,4);
yc2=MCH1CH2pairs(:,3);
%These pairs are used to fit the rest of the points.

%% Checking...
%  figure(3)
%  imshow(imadjust(CH1AvgImg))
%  hold on;
%  plot(ch1peaks(:,2),ch1peaks(:,1),'ro')
%  hold on;
%  plot(Xc1,Yc1,'go')
%  figure(4)
%  imshow(imadjust(CH2AvgImg))
%  hold on;
%  plot(ch2peaks(:,2),ch2peaks(:,1),'ro')
%  hold on;
%  plot(xc2,yc2,'go')
%% 5. ---------------------------------------------------------------------
%% 6. Create the 3rd order 2d polynomial Transformation with the selected 32 points.
B1 =makepol(Xc1,Yc1,3); %must be a self-made function, all possible combinations of all possible orders of z andf y
%this is set of all possible products of x and y, x3, x2y, y3, y2x, etc
%produces a set of 10 coefficiens for each of 32 points
B2 =makepol(xc2,yc2,3); %polynomial to be used to make channel 1 points from channel 2 points
%B-s are 32x10 arrays. 32 is 16x2=number of points, 10 is number of
%coefficients in 2D 3-rd order polynomial
M21x=Xc1'/B2'; %% Xc1=M21x*B2 - OK, but how does it even work, dividing 32x1 vector by 32x10 matrix
M21y=Yc1'/B2'; %% Yc1=M21y*B2
M12x=xc2'/B1'; %% xc2=M12y*B1
M12y=yc2'/B1'; %% yc2=M12y*B1

%% Find Ch2 peaks using CH1 peaks.
%the reason for messing around with the full set of peaks is just for
%visualization purposes. Matrix is based on the 32 peaks
% next is based on all 3567 peaks
fullB1=makepol(ch1peaks(:,2),ch1peaks(:,1),3);
%this one is based on the full set of Ch 1 peaks, without preselection
%I do not think these peaks were actually matched eirher...
predictedx2=(M12x*fullB1')';  %predicted channel 2 peaks based on ALL available channel 1 peaks
predictedy2=(M12y*fullB1')';

%% Find predicted locations of Ch1 peaks using CH2 peaks.
fullB2=makepol(ch2peaks(:,2),ch2peaks(:,1),3);
predictedX1=(M21x*fullB2')';
predictedY1=(M21y*fullB2')';
%% Checking...
% figure(5)
% imshow(imadjust(CH1AvgImg))
% hold on;
% plot(ch1peaks(:,2),ch1peaks(:,1),'ro')
% hold on;
% plot(predictedX1,predictedY1,'go')
% figure(6)
% imshow(imadjust(CH2AvgImg))
% hold on;
% plot(ch2peaks(:,2),ch2peaks(:,1),'ro')
% hold on;
% plot(predictedx2,predictedy2,'go')
%% 6. ---------------------------------------------------------------------
%% 7. Now plot the circles of the predicted circles in CH1 or in CH2.

% Find and remove if any of them are too close to the edge...
x=find(predictedX1>(radius+1)& predictedX1<(size(CH1AvgImg,1)-radius-1));
%find indiceas and values of non-zero elements
predictedX1=uint16(predictedX1(x));
predictedY1=uint16(predictedY1(x));
y=find(predictedY1>(radius+1)& predictedY1<(size(CH1AvgImg,2)-radius-1));
predictedX1=uint16(predictedX1(y));
predictedY1=uint16(predictedY1(y));


%% 3. Plot Image for the full version:
%% 3a. For CH1
axes(handles.beadImage_CH1);
imshow(imadjust(CH1AvgImg)); %shows image without any circles
axis off;
hold on;
scatter(handles.ch1peaks(:,2), handles.ch1peaks(:,1),'MarkerEdgeColor', 'g');
hold on
scatter(predictedX1, predictedY1,'MarkerEdgeColor', 'b');
hold off

%% 3a. -------------------------------
%% 3b. For CH2
axes(handles.beadImage_CH2);
imshow(imadjust(CH2AvgImg));
axis off;
hold on;
scatter(handles.ch2peaks(:,2), handles.ch2peaks(:,1),'MarkerEdgeColor', 'r');
hold on
scatter(predictedx2, predictedy2,'MarkerEdgeColor', 'c');
hold off



%% 8. Preparing to save...
handles.M21x=M21x;
handles.M21y=M21y;
handles.M12x=M12x;
handles.M12y=M12y;
%The _mapped file contains 4 column of data: M21x, M21y, M12x, M12y. where, CH1=M21*CH2 and CH2=M12*CH1
%the first column should contain the coefficients needed to produce CH1 x
%from channel 2 x and y


statement= strcat('The mapping is successfully completed.');
set(handles.status,'String',statement);
drawnow;pause(0.1);

guidata(hObject, handles);

function B =makepol(x,y,n) %NB: first argument is X, or rather a list of x
B=[];
for m=1:length(x)
    if n==3
        A=[1 x(m) y(m) x(m)*y(m) x(m)^2 y(m)^2 x(m)^2*y(m) x(m)*y(m)^2 x(m)^3 y(m)^3];
    end
    if n==4
        A=[1 x(m) y(m) x(m)*y(m) x(m)^2 y(m)^2 x(m)^2*y(m) x(m)*y(m)^2 x(m)^3 y(m)^3 x(m)^3*y(m) x(m)^2*y(m)^2 x(m)*y(m)^3 x(m)^4 y(m)^4];
    end
    
    B=[B;A];
end


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
M21x=handles.M21x;
M21y=handles.M21y;
M12x=handles.M12x;
M12y=handles.M12y;
Fullmodel=[M21x' M21y' M12x' M12y'];
FileTif=handles.FileTif;
pathName=FileTif(1:max(strfind(FileTif,'\')));
fileName=FileTif(max(strfind(FileTif,'\'))+1:end);
%FileExt=fileName(min(strfind(fileName,'.')):end);
fileNameOnly=fileName(1:min(strfind(fileName,'.'))-1);
outfilename=strcat(pathName,fileNameOnly,'_mapped.txt');
save(outfilename,'Fullmodel','-ascii');
A='The _mapped file contains 4 column of data: M21x, M21y, M12x, M12y. where, CH1=M21*CH2 and CH2=M12*CH1';
outfilename2=strcat(pathName,fileNameOnly,'_OutputDes.txt');
fid = fopen(outfilename2,'wt');
fprintf(fid, A);
fclose(fid);
statement= strcat('The mapping parameters are saved successfully.');
set(handles.status,'String',statement);
drawnow;pause(0.1);





function status_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function allowedError_Callback(hObject, eventdata, handles)
% hObject    handle to allowedError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of allowedError as text
%        str2double(get(hObject,'String')) returns contents of allowedError as a double


% --- Executes during object creation, after setting all properties.
function allowedError_CreateFcn(hObject, eventdata, handles)
% hObject    handle to allowedError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peakRadiusCh1_Callback(hObject, eventdata, handles)
% hObject    handle to peakRadiusCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peakRadiusCh1 as text
%        str2double(get(hObject,'String')) returns contents of peakRadiusCh1 as a double


% --- Executes during object creation, after setting all properties.
function peakRadiusCh1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peakRadiusCh1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end























function threshCH1_Callback(hObject, eventdata, handles)
% hObject    handle to threshCH1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshCH1 as text
%        str2double(get(hObject,'String')) returns contents of threshCH1 as a double


% --- Executes during object creation, after setting all properties.
function threshCH1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshCH1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
