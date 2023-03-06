function varargout = SM_Analysis(varargin)
% SM_ANALYSIS MATLAB code for SM_Analysis.fig
%      SM_ANALYSIS, by itself, creates a new SM_ANALYSIS or raises the existing
%      singleton*.
%
%      H = SM_ANALYSIS returns the handle to a new SM_ANALYSIS or the handle to
%      the existing singleton*.
%
%      SM_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_ANALYSIS.M with the given input arguments.
%
%      SM_ANALYSIS('Property','Value',...) creates a new SM_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_Analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_Analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_Analysis

% Last Modified by GUIDE v2.5 16-Nov-2018 13:03:59

% Begin initialization code - DO NOT EDIT
%Shiba: Changed int. at line 222 and 239 to 0.6 from 1.Modified 030718

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SM_Analysis_OpeningFcn, ...
    'gui_OutputFcn',  @SM_Analysis_OutputFcn, ...
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


% --- Executes just before SM_Analysis is made visible.
function SM_Analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_Analysis (see VARARGIN)

% Choose default command line output for SM_Analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_Analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_Analysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function showCurrentFrame(FileTifname,frameNum, handles) %self-made function that reads one frame and shows background-corrected image
%% Reads the current frame of the tif stack images.
currentFrame = imread(FileTifname, 'index', frameNum);
%background = imopen(currentFrame,strel('disk',15));
background=medfilt2(currentFrame,[64,64],'symmetric');
currentFrame=currentFrame-background;

currentFrame= imadjust(im2uint16(currentFrame));
%% end: Reads the first frame of the tif stack images.

%% makes different images for different channel and plots them
howManyColors=2;
if howManyColors==2
    CH1AvgImg=currentFrame(1:handles.nImage,1:handles.mImage/2);
    CH2AvgImg=currentFrame(1:handles.nImage,handles.mImage/2+1: handles.mImage);
    axes(handles.CH1);
    imshow(imadjust(CH1AvgImg)) %show image CH1
    axes(handles.CH2);
    imshow(imadjust(CH2AvgImg)) %show image CH2
end
%% end: makes different images for different channel and plots them










function analyze_Callback(hObject, ~, handles)
%% 1. Gather Info...
%from October 8th 2017 contains the option of working on the group of files
[filelist path fi] = uigetfile('*.tif','Select files on which you wish to perform drift correction.','MultiSelect','on');

%cd(path);

if iscell(filelist) == 0
    filelist2{1} = filelist;
    filelist = filelist2;
    clear filelist2;
end

%common processing parameters for all movies
M_matrix=handles.M_matrix;
radius=str2double(get(handles.radius,'string'));
r_region=radius;
inFrameUB=str2double(get(handles.initialFrameUB,'string'));
ch1Th=str2double(get(handles.ch1Th,'string'));
ch2Th=str2double(get(handles.ch2Th,'string'));
CreateTrace=get(handles.CreateTrace,'value');
edgePad=32;


M21x=M_matrix(:,1); %matrixes of parameters of the polynomials based on beadmapping - same for all files
M21y=M_matrix(:,2);
M12x=M_matrix(:,3);
M12y=M_matrix(:,4);



xx = 1:length(filelist);



%% setup statusbar
axes(handles.statusBar);
xtick=[];
xticklabel=[];
xlim([0 length(filelist)]);
ylim([0 20]);
rectangle('Position',[0,0,length(filelist)+1,21],'FaceColor',[1 1 1])

set(handles.movNumber,'string',['0' '/' num2str(length(filelist))]);




%cycle over multiple file names begins...
for kk = 1:length(filelist) %do for many files
    
    FileTifname = strcat(path,filelist{kk});
    
    linkaxes([handles.CH1, handles.CH2])
    
    %% Update Status:
    statement=['Analyzing: ' filelist{kk} ];
    set(handles.statusbox,'String',statement);
    set(handles.movieName,'String',filelist{kk});
    drawnow;pause(0.1);
    guidata(hObject, handles);
    
    
    
    % ref_movie = input_movie;
    
    %FileTifname=handles.FileTifname;
    
    handles.FileTifname=FileTifname; %just in case - making sure that the handles that used to be in "Open" get re/written
    InfoImage=imfinfo(FileTifname);
    handles.mImage=InfoImage(1).Width;
    handles.nImage=InfoImage(1).Height;
    handles.NumberImages=length(InfoImage);
    % showCurrentFrame(FileTifname,1, handles) %self-made function that shows - show the background-corrected frame of every movie
    guidata(hObject, handles);
    NumberImages=handles.NumberImages;
    
    % make an average image of the first few frames of the whole movie.
    firstFrame= imread(FileTifname, 'index', 1);
    TotImg=uint16(zeros(size(firstFrame)));
    
    
    imageAvg=min(NumberImages,inFrameUB);
    for frame=1:imageAvg
        currentFrame = imread(FileTifname, 'index', frame);
        TotImg=TotImg+currentFrame;
    end %end of for
    TotImg=TotImg./imageAvg;
    
    Iblur1 = imgaussfilt(TotImg,1);
    Iblur2 = imgaussfilt(TotImg,10);
    Bsub_TotImg=Iblur1-Iblur2;
    Img=imadjust(Bsub_TotImg);
    outfilename=[FileTifname(1:max(strfind(FileTifname,'.'))-1) '_img.tif'];
    imwrite(mat2gray(Img),outfilename)
    
    
    %% Sujay Used this commented part to test the effect of background blurring on the number of peaks.
    %     figure (4)
    %     h(1)= subplot(1,2,1)
    %     imshow(imadjust(TotImg)) %show image CH1
    %     D1=double(TotImg);
    %     thres = (max([min(max(D1,[],1))  min(max(D1,[],2))]))*0.1;
    %     p=FastPeakFind(TotImg,thres,fspecial('disk',2), 50);
    %     chpeaks=reshape(p',2,numel(p)/2);
    %     chpeaks=chpeaks';
    %     Ch1base=[chpeaks(:,2) chpeaks(:,1)];
    %     hold on;
    %     plot(Ch1base(:,2),Ch1base(:,1),'go');
    %     hold off;
    %
    %
    %     h(2)=subplot(1,2,2)
    %     imshow(imadjust(Img))
    %     linkaxes(h)
    %     D1=double(Img);
    %     thres = (max([min(max(D1,[],1))  min(max(D1,[],2))]))*0.1
    %     p=FastPeakFind(Img,thres,fspecial('disk',2), 50);
    %     chpeaks=reshape(p',2,numel(p)/2);
    %     chpeaks=chpeaks';
    %     Ch1base=[chpeaks(:,2) chpeaks(:,1)];
    %     hold on;
    %     plot(Ch1base(:,2),Ch1base(:,1),'go');
    %     hold off;
    %% -------------------------------------
    
    
    
    %% Update Status:
    statement='Creating the .tif image.';
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
    %% ---------------------------
    
    %% Seperate image as two channels
    CH1TotImg=Img(:,1:handles.mImage/2);
    CH2TotImg=Img(:,handles.mImage/2+1:handles.mImage);
    %% ----------
    
    %% Find peaks independently in CH1.
    D1=double(CH1TotImg);
%     thres = (max([min(max(D1,[],1))  min(max(D1,[],2))]))*0.4
%     p=FastPeakFind(CH1TotImg,thres,fspecial('disk',2), 32);
%     chpeaks=reshape(p',2,numel(p)/2);
%     chpeaks=chpeaks';
%     allPeaksCh1=[chpeaks(:,2) chpeaks(:,1)];
    
    p= PeakFinderCentroid(CH1TotImg,ch1Th,edgePad,'y',radius);
    p = p';
    yc1 = p(:,2);
    xc1 = p(:,1);
    allPeaksCh1=[xc1 yc1];
    
    
    
    
    
    %% Find peaks independently in CH2.
    clear p chpeaks;
    D2=double(CH2TotImg);
    D2(D2==0) = NaN;
   % thres = (max([min(max(D2,[],1))  min(max(D2,[],2))]))*0.4;
   % p=FastPeakFind(CH2TotImg,thres,fspecial('disk',2), 32);
%     chpeaks=reshape(p',2,numel(p)/2);
%     chpeaks=chpeaks';
%     allPeaksCh2=[chpeaks(:,2) chpeaks(:,1)];
%    

   p= PeakFinderCentroid(CH2TotImg,ch2Th,edgePad,'y',radius);
    p = p';
    yc2 = p(:,2);
    xc2 = p(:,1);
    allPeaksCh2=[xc2 yc2];

    
    %% Plot individual peaks on CH1 and CH2. These are All the Peaks found in both Channels independently.
    axes(handles.CH1);
    imshow(imadjust(CH1TotImg)) %show image CH1
    hold on;
    plot(allPeaksCh1(:,2),allPeaksCh1(:,1),'go');
    hold off;
    
    axes(handles.CH2);
    imshow(imadjust(CH2TotImg)) %show image CH2
    hold on;
    plot(allPeaksCh2(:,2),allPeaksCh2(:,1),'ro');
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now we find the common peaks and the CH1 only and CH2 only peaks.
    %fixedChanel= Ch2
    %variable=Ch1
    % First, Map all peaks of CH2 to CH1. We call the mapped peaks as predicted
    % CH1: preCH1.
    
    B2=makepol(allPeaksCh2(:,2),allPeaksCh2(:,1),3); %arguments go in x, y order
    Xc1=(M21x'*B2')';   %% Xc1=M21x*B2
    Yc1=(M21y'*B2')';
    preCh1=[Yc1 Xc1];  %predicted peaks in Ch1 by the peaks in Ch2.
    %just changing order x-y since first index in images is row, i.e. y
    %predicted locations of peaks in channel 1 based on locations in channel 2
    
  %% We create a Info matrix to gather all informations about the peaks.
  %% Colulmn1: Distance between Predicted peak in CH1 and the closest real peak in CH1.
  %% Column2: indices of the CH2 peaks in order of allPeaksCh2
  %% Column3 & 4: Coordinates of predicted CH1 peaks Y (column3) and X(column4)
  %% Column5: indices of the CH1 peaks in order of allPeaksCh1
  %% Column6 & 7: Coordinates of real CH1 peaks Y (column6) and X(column7)
  
    Info=zeros(length(preCh1),7);
    for g=1: length(preCh1)
        C=preCh1(g,:); %point (x,y)
        C2ind=g;
        D = pdist2(allPeaksCh1,C);
        [val,C1ind]= min(D);
        Info(g,:)=[val C2ind C C1ind allPeaksCh1(C1ind,:)];
    end
    
    %% Common peaks are those whose centers falls within the radius of the predicted peak.
    CommonPeaks=Info(Info(:,1)<=radius,:);
    %% Some peaks are mapped to multiple Real peaks of CH1. Filter those peaks.
    CommonPeaks=sortrows(CommonPeaks,[5,1]);
    [~,ia,~] = unique(CommonPeaks(:,5));
    %% Keep Only the unique peaks
    NewCommonPeaks=CommonPeaks(ia,:);
    NewCommonPeaks=sortrows(NewCommonPeaks,5);
   %% NewCommon Peaks are all the common peaks that is found in both chanels. 
    %% Rest of the peaks do not have a pair in CH1
    Temp=ones(length(allPeaksCh2),1);
    for i=1:length(NewCommonPeaks)
        Temp(NewCommonPeaks(i,2),1)=0;
    end
    nonZero=find(Temp);
    morepeaks2=Info(nonZero,:); %% Only the mapped values are important not the nearest values in CH1.
    % more peaks that do not match with peaks in CH1
    % So Just find values in the mapped positions. and ignore the nearest
    % particles.2342
    
    morepeaksCH2= allPeaksCh2(morepeaks2(:,2),:);
    
    clear Temp;
    Temp=ones(length(allPeaksCh1),1);
    for i=1:length(NewCommonPeaks)
        Temp(NewCommonPeaks(i,5),1)=0;
    end
    nonZero=find(Temp);   
    morepeaksCH1= allPeaksCh1(nonZero,:); % this peaks dso not have a matching pair in CH2.
    % So Just find values in the mapped positions. and ignore the nearest particles.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% So Now the list of CH1 peaks would be: NewCommonPeaks CH1s, morepeaksCH1 and mapped morepeaksCH2
    CH1CommonPeaks=[NewCommonPeaks(:,6) NewCommonPeaks(:,7)];
    %morepeaksCH1=morepeaksCH1;
    mappedpeaksfromCH2=[morepeaks2(:,3) morepeaks2(:,4)];
    %AllCH1Peaks=[CH1CommonPeaks;morepeaksCH1;mappedpeaksfromCH2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% and list of CH2 peaks would be: NewCommonPeaks CH2s, mapped morepeaksCH1 (we need mapping for this) and morepeaksCH2
    CH2CommonPeaks=allPeaksCh2(NewCommonPeaks(:,2),:);
    B1=makepol(morepeaksCH1(:,2),morepeaksCH1(:,1),3);
    Xc1=(M12x'*B1')';   %% Xc1=M21x*B2
    Yc1=(M12y'*B1')';
    mappedpeaksfromCH1=[Yc1 Xc1];  %predicted peaks in Ch1 by the peaks in Ch2.
    % morepeaksCH2=morepeaksCH2;
    %AllCH2Peaks=[CH2CommonPeaks;mappedpeaksfromCH1;morepeaksCH2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    AnCriteria=get(handles.AnCriteria,'Value');
    if AnCriteria==1
        % These are the mapped peaks in CH1. So if we reverse map them, we should get proper peaks in CH2.
        %These are the peaks in the CH1 Channel
        AllCH1Peaks=CH1CommonPeaks;
        AllCH2Peaks=CH2CommonPeaks;
        
        %% Plot individual peaks on CH1 and CH2
        axes(handles.CH1);
        imshow(imadjust(CH1TotImg)) %show image CH1
        hold on;
        plot(allPeaksCh1(:,2),allPeaksCh1(:,1),'go');
        plot(CH1CommonPeaks(:,2),CH1CommonPeaks(:,1),'yo');
        %  plot(morepeaksCH1(:,2),morepeaksCH1(:,1),'co');
        % plot(mappedpeaksfromCH2(:,2),mappedpeaksfromCH2(:,1),'mo');
        hold off;
        
        axes(handles.CH2);
        imshow(imadjust(CH2TotImg)) %show image CH2
        hold on;
        plot(allPeaksCh2(:,2),allPeaksCh2(:,1),'ro');
        plot(CH2CommonPeaks(:,2),CH2CommonPeaks(:,1),'yo');
        %  plot(mappedpeaksfromCH1(:,2),mappedpeaksfromCH1(:,1),'co');
        % plot(morepeaksCH2(:,2),morepeaksCH2(:,1),'mo');
        hold off;
        
    elseif AnCriteria==2
        %% Only the CH1 peaks are analyzed...
        AllCH1Peaks=[CH1CommonPeaks;morepeaksCH1];
        AllCH2Peaks=[CH2CommonPeaks;mappedpeaksfromCH1];
        
        %% Plot individual peaks on CH1 and CH2
        axes(handles.CH1);
        imshow(imadjust(CH1TotImg)) %show image CH1
        hold on;
        plot(allPeaksCh1(:,2),allPeaksCh1(:,1),'go');
        plot(CH1CommonPeaks(:,2),CH1CommonPeaks(:,1),'yo');
        plot(morepeaksCH1(:,2),morepeaksCH1(:,1),'go');
        % plot(mappedpeaksfromCH2(:,2),mappedpeaksfromCH2(:,1),'mo');
        hold off;
        
        axes(handles.CH2);
        imshow(imadjust(CH2TotImg)) %show image CH2
        hold on;
        plot(allPeaksCh2(:,2),allPeaksCh2(:,1),'ro');
        plot(CH2CommonPeaks(:,2),CH2CommonPeaks(:,1),'yo');
        plot(mappedpeaksfromCH1(:,2),mappedpeaksfromCH1(:,1),'go');
        % plot(morepeaksCH2(:,2),morepeaksCH2(:,1),'mo');
        hold off;
    elseif AnCriteria==3
        %% Only the CH2 peaks are analyzed. CH2 peaks are first Mapped to CH1 and then those predicted CH1 peaks are considered the real peaks...
        AllCH1Peaks=[CH1CommonPeaks;mappedpeaksfromCH2];
        AllCH2Peaks=[CH2CommonPeaks;morepeaksCH2];
        
        %% Plot individual peaks on CH1 and CH2
        axes(handles.CH1);
        imshow(imadjust(CH1TotImg)) %show image CH1
        hold on;
        plot(allPeaksCh1(:,2),allPeaksCh1(:,1),'go');
        plot(CH1CommonPeaks(:,2),CH1CommonPeaks(:,1),'yo');
        %  plot(morepeaksCH1(:,2),morepeaksCH1(:,1),'co');
        plot(mappedpeaksfromCH2(:,2),mappedpeaksfromCH2(:,1),'ro');
        hold off;
        
        axes(handles.CH2);
        imshow(imadjust(CH2TotImg)) %show image CH2
        hold on;
        plot(allPeaksCh2(:,2),allPeaksCh2(:,1),'ro');
        plot(CH2CommonPeaks(:,2),CH2CommonPeaks(:,1),'yo');
        %  plot(mappedpeaksfromCH1(:,2),mappedpeaksfromCH1(:,1),'co');
        plot(morepeaksCH2(:,2),morepeaksCH2(:,1),'ro');
        hold off;
    elseif AnCriteria==4
        AllCH1Peaks=[CH1CommonPeaks;morepeaksCH1;mappedpeaksfromCH2];
        AllCH2Peaks=[CH2CommonPeaks;mappedpeaksfromCH1;morepeaksCH2];
        
        %% Plot individual peaks on CH1 and CH2
        axes(handles.CH1);
        imshow(imadjust(CH1TotImg)) %show image CH1
        hold on;
        plot(allPeaksCh1(:,2),allPeaksCh1(:,1),'go');
        plot(CH1CommonPeaks(:,2),CH1CommonPeaks(:,1),'yo');
        plot(morepeaksCH1(:,2),morepeaksCH1(:,1),'go');
        plot(mappedpeaksfromCH2(:,2),mappedpeaksfromCH2(:,1),'ro');
        hold off;
        
        axes(handles.CH2);
        imshow(imadjust(CH2TotImg)) %show image CH2
        hold on;
        plot(allPeaksCh2(:,2),allPeaksCh2(:,1),'ro');
        plot(CH2CommonPeaks(:,2),CH2CommonPeaks(:,1),'yo');
        plot(mappedpeaksfromCH1(:,2),mappedpeaksfromCH1(:,1),'go');
        plot(morepeaksCH2(:,2),morepeaksCH2(:,1),'ro');
        hold off;
    end
  
    AllPeaks=[AllCH1Peaks AllCH2Peaks];
    AllPeaks=round(AllPeaks);
    
    
    %% We got all points in CH1 and CH2 and ready to find integral intensity at those points:
    nummols=size(AllPeaks,1);
    
    %% Update Status:
    statement=['Peaks found:  ' num2str(nummols)];
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
    
    set(handles.commonPeak,'String',num2str(length(CH2CommonPeaks)));
    drawnow;pause(0.1);
    
    set(handles.CH1Peak,'String',num2str(length(mappedpeaksfromCH1)));
    drawnow;pause(0.1);
    set(handles.CH2peak,'String',num2str(length(morepeaksCH2)));
    drawnow;pause(0.1);
    
    set(handles.totalAnalyzedPeak,'String',num2str(nummols));
    drawnow;pause(0.1);
    guidata(hObject, handles);
 
    pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Peaks found. Now calculate the intensity and make traces.
    
    
    donor=zeros(nummols,NumberImages);
    acceptor=zeros(nummols,NumberImages);
    
    x1= AllPeaks(:,2); %x in the Ch1
    y1= AllPeaks(:,1); %y in the Ch1
    x2= AllPeaks(:,4); %x in the Ch2
    y2= AllPeaks(:,3); %y in the Ch2
    
    
  tic  
    aa=0; %flags for too large signal values
 if CreateTrace==0   
    
    for frame=1:NumberImages
        currentFrame = imread(FileTifname, 'index', frame);
        %% Split into frames.
        CH1AvgImg=currentFrame(1:handles.nImage,1:handles.mImage/2);
        CH2AvgImg=currentFrame(1:handles.nImage,handles.mImage/2+1: handles.mImage);
        
        % Create a background image by setting the molecules to NAN.
        Ch1B=double(CH1AvgImg);
        Ch2B=double(CH2AvgImg);
        
        for m=1:nummols
            Ch1B(y1(m)-3:y1(m)+3,x1(m)-3:x1(m)+3)=nan;
            Ch2B(y2(m)-3:y2(m)+3,x2(m)-3:x2(m)+3)=nan;
        end
        
        for m=1:nummols
            %% CH1            
            region = CH1AvgImg(y1(m)-r_region:y1(m)+r_region,x1(m)-r_region:x1(m)+r_region); %first index is row, aka y.
            bg=(Ch1B(y1(m)-2*r_region:y1(m)+2*r_region,x1(m)-2*r_region:x1(m)+2*r_region));            
            background = nanmedian(reshape(bg,1,size(bg,1)*size(bg,2)));
            particle = double(region(r_region-2:r_region+4,r_region-2:r_region+4));
            particle = particle-background;
            I = sum(sum(particle));
            donor(m,frame) = I;
            if I>65536
                aa=1;
            end
            donor(m,frame) = I;
            %% CH2
            region = CH2AvgImg(y2(m)-r_region:y2(m)+r_region,x2(m)-r_region:x2(m)+r_region);
            bg=(Ch2B(y2(m)-2*r_region:y2(m)+2*r_region,x2(m)-2*r_region:x2(m)+2*r_region));
            background = nanmedian(reshape(bg,1,size(bg,1)*size(bg,2)));
            particle = double(region(r_region-2:r_region+4,r_region-2:r_region+4));
            particle = particle-background;
            I = sum(sum(particle));
            acceptor(m,frame) = I;
            if I>65536
                aa=1;
            end
            acceptor(m,frame) = I;
        end % end of cycle over number of molecules
        
        if mod(frame,10)==0  
            %% Update Status:
            statement=['Done with frame:  ' num2str(frame)];
            set(handles.statusbox,'String',statement);
            drawnow;pause(0.1);
            guidata(hObject, handles);
        end
    end %end of cycle over number of frames
    toc
    %
    if aa==1 %if some signals do not fit into 65536 (beads)
        donor=int16(donor/10); %matters only for beads. Prevents problems with bright bead signals not fitting into
        %16-bit integer
        acceptor=int16(acceptor/10);
    else
        donor=int16(donor);
        acceptor=int16(acceptor);
    end
    %
    peaktable = zeros(size(AllPeaks,1)*2,3);
    traces = zeros(size(AllPeaks,1)*2,size(donor,2));
    AllPeaks = round(AllPeaks);
    %one element is 4-element string
    %first two elements are coordinates of peaksToRecord1, third and fourth are coordinates of peaksToRecord2
    
    % Organize peaktable and traces in alternating rows
    %size_m_pairs=size(donor,2) %this is the length of the donor trace. Should be equal to number of frames
    for i = 1:size(AllPeaks,1)*2
        peaktable(i,1) = i; %pks file has just three columns, the first is number 1,2,3,,,,
        if mod(i,2) ~= 0 %mod(a,2) returns the remainder after division of a by 2, not zero - here are ODD numbers
            peaktable(i,2:3) = cat(2,AllPeaks((i+1)/2,2),AllPeaks((i+1)/2,1));
            traces(i,:) = donor((i+1)/2,:); %first index - peak, second index: frame
        else % even numbers
            peaktable(i,2:3) = cat(2,AllPeaks(i/2,4),AllPeaks(i/2,3));
            traces(i,:) = acceptor(i/2,:);
        end
    end
    traces = cat(1,1:size(traces,2),traces);
    
    %% Update Status:
    statement=['Writing the Pks file.'];
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
    
    outfilename=[FileTifname(1:max(strfind(FileTifname,'.'))-1) '.pks'];
    save(outfilename, 'peaktable', '-ascii');
    pause(1)
    
    %% Update Status:
    statement='Writing the .traces file.';
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
    
    len = size(traces,2); %size of second dimension of traces
    Ntraces = size(traces,1)-1; %size of the first dimension of traces minus 1
    outfilename=[FileTifname(1:max(strfind(FileTifname,'.'))-1) '.traces'];
    fid = fopen(outfilename,'w','ieee-le');
    %'l' or 'ieee-le' Little-endian BYTE ordering
    fwrite(fid,len,'int32');
    fwrite(fid,Ntraces,'int16');
    fwrite(fid,traces,'int16');
    fclose(fid);

     end 
    
    pause(1)
    
    %% Update Status:
    statement=['Done with movie: ' filelist{kk}];
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
     
    
    %% setup statusbar
    axes(handles.statusBar);
    xtick=[];
    xticklabel=[];
    xlim([0 length(filelist)])
    axis off
    for j=1:1000
        line([(xx(kk)-1)+(j-1)*0.001 (xx(kk)-1)+(j-1)*0.001],[0 20],'Color','green');
    end
    pause(1);
    drawnow
    
    set(handles.movNumber,'string',[num2str(kk) '/' num2str(length(filelist))]);
    
 
end %end of cycle over multiple files to be processed.







function B =makepol(x,y,n)
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






function true_threshold=FindThreshold(Img, rad, min_thr)
%uses background-corrected image, background close to zero
Alb=0.05;
Aub=0.10;
area_size=(2*rad+1)^2;
% im_size=length(Img); - this works only for square images
[width,height,z1] = size(Img);
im_size=width*height; %NB im_size is now redefined as image area!!!
imageMax=max(max(Img));
imageMin=1;
%temp_ints=reshape(double(Img),[1,im_size^2]);
temp_ints=reshape(double(Img),[1,im_size]);
%[X,h]=histcounts(temp_ints,'Normalization','probability','BinWidth',1);
%histogram of all possible intensities in the image
%there must be large percentage of zeros... almost 34% in our tests
temp_ints_sort=sort(temp_ints,'descend');

temp_thresh_initial=temp_ints_sort(50)*.1;
%this is the same type of threshold as used in beadmapper...
% one tenth of a fiftiest-highest intensity
fraction_area=(Alb+Aub)/2;
counter=0; %infinite loop protection; one gets out of while after 50 iterations
while (fraction_area>Alb && fraction_area<Aub )%looking for a rather low threshold, peaks covering 10% of total area
    [cent,~]=FastPeakFind(Img,temp_thresh_initial);
    %does fast peak finding with current threshold
    %if input threshold too low, will find many peaks and the area covered
    %with peaks will be high.
    num_points=length(cent)/2; %length of larger dimension - number of peaks
    area_points=area_size*num_points; %area covered by peaks; in pixels
    %fraction_area=area_points/(im_size^2); %area in pixels as fraction of image area covered by peaks
    fraction_area=area_points/(im_size);
    if fraction_area>Aub %if initial guess for threshold was high
        imageMin=temp_thresh_initial;
        temp_thresh_initial=(imageMax+imageMin)/2;
        fraction_area=(Alb+Aub)/2;
        counter=counter+1;
    elseif fraction_area<Alb %if initial guess for threshold is too small, and too many peaks are found, try again with threshold increased by 20%
        imageMax=temp_thresh_initial;
        temp_thresh_initial=(imageMin+imageMax)/2;
        fraction_area=(Alb+Aub)/2;
        counter=counter+1;
    elseif fraction_area>Alb && fraction_area<Aub
        fraction_area=(Alb+Aub);
        true_threshold =temp_thresh_initial;
    end
    %by the end of the movie there are very few peaks left, so the program
    %would decrease threshold too much trying to fill 10% of the area with
    %peaks... Need to prevent it from doing that.
    if counter>50 %infinite loop protection
        fraction_area=(Alb+Aub);
        true_threshold =temp_thresh_initial;
    end
end %end of while
if temp_thresh_initial<min_thr %if despite all efforts it aimed at too small threshold
    true_threshold=min_thr;
end %end of if

% counter=0; %infinite loop protection
% while (fraction_area>Alb && fraction_area<Aub )%looking for a rather low threshold, peaks covering 10% of total area
%     [cent,~]=FastPeakFind(Img,temp_thresh_initial);
%     %does fast peak finding with current threshold
%     %if input threshold too low, will find many peaks and the area covered
%     %with peaks will be high.
%     num_points=length(cent)/2; %length of larger dimension - number of peaks
%     area_points=area_size*num_points; %area covered by peaks; in pixels
%     fraction_area=area_points/(im_size^2); %area in pixels as fraction of image area covered by peaks
%     if fraction_area>Aub %if initial guess for threshold was high
%        temp_thresh_initial=temp_thresh_initial*1.05;
%        fraction_area=(Alb+Aub)/2;
%        counter=counter+1;
%     elseif fraction_area<Alb %if initial guess for threshold is too small, and too many peaks are found, try again with threshold increased by 20%
%        temp_thresh_initial=temp_thresh_initial/1.05;
%         fraction_area=(Alb+Aub)/2;
%         counter=counter+1;
%     elseif fraction_area>Alb && fraction_area<Aub
%         fraction_area=(Alb+Aub);
%         true_threshold =temp_thresh_initial;
%     end
%     if counter>50
%          fraction_area=(Alb+Aub);
%         true_threshold =temp_thresh_initial;
%     end
%         %by the end of the movie there are very few peaks left, so the program
%     %would decrease threshold too much trying to fill 10% of the area with
%     %peaks... Need to prevent it from doing that.
% end %end of while
% if temp_thresh_initial<min_thr %if despite all efforts it aimed at too small threshold
%     true_threshold=min_thr;
% end %end of if







function loadMap_Callback(hObject, ~, handles)
DefaultFilePath='Z:\Valter Z\September 26, 2017';
[fileName,pathName]=uigetfile('*.txt','Select the mapping file',DefaultFilePath);
name=strcat(pathName,fileName);
handles.beadFileName=name;
A = importdata(name,' ');
handles.M_matrix=A;
showName=strcat(pathName(1:5),'...\',fileName);
set(handles.showName,'String',showName);
guidata(hObject, handles); %set handles to access from other functions.



function initialFrameLB_Callback(hObject, eventdata, handles)
% hObject    handle to initialFrameLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initialFrameLB as text
%        str2double(get(hObject,'String')) returns contents of initialFrameLB as a double


% --- Executes during object creation, after setting all properties.
function initialFrameLB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initialFrameLB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initialFrameUB_Callback(hObject, eventdata, handles)
% hObject    handle to initialFrameUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initialFrameUB as text
%        str2double(get(hObject,'String')) returns contents of initialFrameUB as a double


% --- Executes during object creation, after setting all properties.
function initialFrameUB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initialFrameUB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function radius_Callback(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of radius as text
%        str2double(get(hObject,'String')) returns contents of radius as a double


% --- Executes during object creation, after setting all properties.
function radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function statusbox_Callback(hObject, eventdata, handles)
% hObject    handle to statusbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statusbox as text
%        str2double(get(hObject,'String')) returns contents of statusbox as a double


% --- Executes during object creation, after setting all properties.
function statusbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function movieName_Callback(hObject, eventdata, handles)
% hObject    handle to movieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of movieName as text
%        str2double(get(hObject,'String')) returns contents of movieName as a double


% --- Executes during object creation, after setting all properties.
function movieName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to movieName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in AnCriteria.
function AnCriteria_Callback(hObject, eventdata, handles)
% hObject    handle to AnCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AnCriteria contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AnCriteria


% --- Executes during object creation, after setting all properties.
function AnCriteria_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AnCriteria (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch1Th_Callback(hObject, eventdata, handles)
% hObject    handle to ch1Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch1Th as text
%        str2double(get(hObject,'String')) returns contents of ch1Th as a double


% --- Executes during object creation, after setting all properties.
function ch1Th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch1Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ch2Th_Callback(hObject, eventdata, handles)
% hObject    handle to ch2Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ch2Th as text
%        str2double(get(hObject,'String')) returns contents of ch2Th as a double


% --- Executes during object creation, after setting all properties.
function ch2Th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ch2Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CreateTrace.
function CreateTrace_Callback(hObject, eventdata, handles)
% hObject    handle to CreateTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CreateTrace
