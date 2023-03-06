function varargout = SM_DriftCorrection(varargin)
% SM_DRIFTCORRECTION MATLAB code for SM_DriftCorrection.fig
%      SM_DRIFTCORRECTION, by itself, creates a new SM_DRIFTCORRECTION or raises the existing
%      singleton*.
%
%      H = SM_DRIFTCORRECTION returns the handle to a new SM_DRIFTCORRECTION or the handle to
%      the existing singleton*.
%
%      SM_DRIFTCORRECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_DRIFTCORRECTION.M with the given input arguments.
%
%      SM_DRIFTCORRECTION('Property','Value',...) creates a new SM_DRIFTCORRECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_DriftCorrection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_DriftCorrection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_DriftCorrection

% Last Modified by GUIDE v2.5 14-Nov-2017 10:46:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SM_DriftCorrection_OpeningFcn, ...
    'gui_OutputFcn',  @SM_DriftCorrection_OutputFcn, ...
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


% --- Executes just before SM_DriftCorrection is made visible.
function SM_DriftCorrection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_DriftCorrection (see VARARGIN)

% Choose default command line output for SM_DriftCorrection
handles.output = hObject;

axes(handles.axes1);
xlabel('pixel movement(x)')
    ylabel('pixel movement(y)')
    title('Only consider integer steps.')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_DriftCorrection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_DriftCorrection_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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



% --------------------------------------------------------------------
function newFile_ClickedCallback(hObject, eventdata, handles)
load_Callback(hObject, eventdata, handles)



% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)

%--Drift Correction Settings----
smoothframes = str2double(get(handles.smoothFrame,'string')); % Number of frames over which to perform rolling average of drift coordinates
max_drift = str2double(get(handles.maxDrift,'string'));
maxRelShift = str2double(get(handles.maxRelShift,'string'));
delOriginal=get(handles.delOriginal,'value');

[filelist path fi] = uigetfile('*.tif','Select files on which you wish to perform drift correction.','MultiSelect','on');

%cd(path);

if iscell(filelist) == 0
    filelist2{1} = filelist;
    filelist = filelist2;
    clear filelist2;
end




%% Update Status:
statement='Movie loaded successfully.';
set(handles.statusbox,'String',statement);
drawnow;pause(0.1);
guidata(hObject, handles);


xx = 1:length(filelist);

%% setup statusbar
axes(handles.statusBar);
xtick=[];
xticklabel=[];
xlim([0 length(filelist)]);
ylim([0 20]);
rectangle('Position',[0,0,length(filelist)+1,21],'FaceColor',[1 1 1])

set(handles.movNumber,'string',['0' '/' num2str(length(filelist))]);





for kk = 1:length(filelist) %do for many files
    
    input_movie = strcat(path,filelist{kk});
    
    ref_movie = input_movie;
    
    %--------------------------------
    
    %--Calculate Drift---------------
    ref_info = imfinfo(ref_movie,'tif');
    frames = size(ref_info,1);
    xdim = ref_info(1,1).Width;
    ydim = ref_info(1,1).Height;
     %using only a part of the region to be FTransformed. Center
        %(quarter of image, area-wise) of the green channel
        min_pixel = 1; %% 1
        max_pixel = ydim;
    
    im = zeros(ydim,xdim,smoothframes);
    for p = 1:smoothframes
        im(:,:,p) = imread(ref_movie,'tif',p);
    end
    % start_frame = imread(ref_movie,'tif',1);
    start_frame = mean(im,3);
    start_frame2 = start_frame - imopen(start_frame,strel('disk',4));
    start_fft=fft2(start_frame2(min_pixel:max_pixel,min_pixel:max_pixel));
    ref_fft=start_fft;
    
    % imshow(start_frame2/max(max(start_frame2)));
    % pause;
    
    driftxy = [0 0];
    
    
    %% Update Status:
    statement='Working on frame: 1';
    set(handles.statusbox,'String',statement);
    set(handles.movieName,'String',filelist{kk});
    drawnow;pause(0.1);
    guidata(hObject, handles);
    lastgoodbackup(1,1)=0;
    lastgoodbackup(1,2)=0;
    
    for i = 2:frames-smoothframes+1 %for every file, go through frames
        if mod(i,10)==0
            
            %% Update Status:
            statement=strcat('Working on frame:', num2str(i+1));
            set(handles.statusbox,'String',statement);
            drawnow;pause(0.1);
            guidata(hObject, handles);
            
        end
        if i == 2 %first frame is already created. Second frame requires special treatment in case of smoothing?
            im = zeros(ydim,xdim,smoothframes); %% Sujay changed it from zeros(xdim,ydim,smoothframes);
            for q = 2:smoothframes+1
                im(:,:,q) = imread(ref_movie,'tif',q);
            end
        else
            imnew = imread(ref_movie,'tif',i+smoothframes-1);
            im = cat(3,im(:,:,2:end),imnew);
        end
        test_frame = mean(im,3);
        %     test_frame = imread(ref_movie,'tif',i);
        test_frame2 = test_frame - imopen(test_frame,strel('disk',4)); %4 is radius. 
        %The objects smaller than 9 pixels in diameter get removed by
        %imopen
        %so this is high frequency random noise filter rather than
        %background correction
        %then the next subtraction creats an image with only sharp peaks
        
        %     output = dftregistration(fft2(test_frame(128:384,128:384)),fft2(start_frame(128:384,128:384)),100);
       test_fft=fft2(test_frame2(min_pixel:max_pixel,min_pixel:max_pixel));
       %this is current frame
        output = dftregistration(test_fft,start_fft,100); %absolute shifts with respect to the first frame
        %compares current frame to start frame, therefore the offset is cumulative absolute shift  
        absoffset(1,1) = output(1,4); % These are absolute shifts
        absoffset(1,2) = output(1,3); %X-Y shift is done here.  
        reloffset(1,1) = absoffset(1,1)-driftxy(end,1);
        reloffset(1,2) = absoffset(1,2)-driftxy(end,2);
        if or(abs(absoffset(1,1)) > max_drift,abs(absoffset(1,2))>max_drift)
           % if the absolute shift is suspiciously large 
            if or(abs(reloffset(1,1)) > 1,abs(reloffset(1,2))>1) %more likely crazy, if relative shift is crazy too
                output = dftregistration(test_fft,ref_fft,100); %underestimated relative shift
                reloffset(1,1) = output(1,4); % These are underestimated relative shifts 
                reloffset(1,2) = output(1,3); %X-Y shift is done here.
                %the third variable returned by dft is rows i.e. Y
                if or(abs(reloffset(1,1)) > 1,abs(reloffset(1,2))>1) %most definitely crazy
%                     %if underestimated shift between neighbor frames is too
%                     %large
                    lastgoodbackup(1,1)=driftxy(end,1); %just keep the previous shifts and hope for the best
                    lastgoodbackup(1,2)=driftxy(end,2);
                else
                    lastgoodbackup(1,1)=driftxy(end,1)+reloffset(1,1); % correct a little
                    lastgoodbackup(1,2)=driftxy(end,2)+reloffset(1,2);
                    ref_fft=test_fft; %resets the reference so we compare with this frame next time
                end
            else %maybe the "honest" absolute shift is just too large
                lastgoodbackup(1,1)=reloffset(1,1)+driftxy(end,1); %add relative shift to the last absolute shift
                lastgoodbackup(1,2)=reloffset(1,2)+driftxy(end,2);
                ref_fft=test_fft; %resets the reference so we compare with this frame next time
           end
       else %meaning if absolute shifts are small enough (likely no crazyness)
            %but what if accidentally absolute shift stays within +/-max_drift but
            %relative shift is too large?
            if or(abs(reloffset(1,1)) > maxRelShift,abs(reloffset(1,2))>maxRelShift) %likely crazy, if relative shift is crazy 
                output = dftregistration(test_fft,ref_fft,100); %underestimated relative shift
                reloffset(1,1) = output(1,4); % These are underestimated relative shifts 
                reloffset(1,2) = output(1,3); %X-Y shift is done here.
                if or(abs(reloffset(1,1)) > maxRelShift,abs(reloffset(1,2))>maxRelShift) %most definitely crazy
%                     if underestimated shift between neighbor frames is too
%                     large
                    lastgoodbackup(1,1)=driftxy(end,1); %just keep the previous shifts and hope for the best
                    lastgoodbackup(1,2)=driftxy(end,2);
                else
                    lastgoodbackup(1,1)=driftxy(end,1)+reloffset(1,1); % correct a little
                    lastgoodbackup(1,2)=driftxy(end,2)+reloffset(1,2);
                    ref_fft=test_fft; %resets the reference so we compare with this frame next time
               end
            else %everything perfect
                lastgoodbackup(1,1)=absoffset(1,1); %absolute shifts are likely reasonable, backups them
                lastgoodbackup(1,2)=absoffset(1,2);
                ref_fft=test_fft; %resets the reference since otherwise it stays the first frame
           end
        end
        driftxy = cat(1,driftxy,lastgoodbackup); %absolute shifts
     end
    
    framestoadd(1:smoothframes-1,1) = driftxy(end,1); %last frames in case of smoothing have to be handled separately
    framestoadd(1:smoothframes-1,2) = driftxy(end,2);
    
    driftxy = cat(1,driftxy,framestoadd);
       
    % driftxy = fliplr(driftxy); %Test whether flipping x & y improves results
    driftx = driftxy(:,1);
    driftx_smooth = smooth(driftx,smoothframes);
    drifty = driftxy(:,2);
    drifty_smooth = smooth(drifty,smoothframes);
    
    c = linspace(1,5,length(driftx));
    axes(handles.axes1);
    scatter(driftx(1:1:end),drifty(1:1:end),38,c,'filled');
    xlabel('pixel movement(x)')
    ylabel('pixel movement(y)')
    title('Only consider integer steps.')
    hold on
    plot(driftx_smooth(1:end),drifty_smooth(1:end),'r-->');
    % axis square;
    % axis equal;
    axis normal;
    hold off
    %--------------------------------
    
    %driftxy=round(driftxy); 
    driftx=round(driftx);
    drifty=round(drifty);
    %need to do it here, so first shifts could accumulate, but the end result would be integer
    % pause;
    %Perform drift correction on input movie and save corrected movie as output
    
    
    %% Update Status:
    statement='Performing drift correction on and writing frame: 1';
    set(handles.statusbox,'String',statement);
    drawnow;pause(0.1);
    guidata(hObject, handles);
    
    
    for k = 1:frames
        if mod(k,10)==0
            %% Update Status:
            statement=strcat('writing frame:',num2str(k+1));
            set(handles.statusbox,'String',statement);
            drawnow;pause(0.1);
            guidata(hObject, handles);
        end
        I = imread(input_movie,'tif',k);
        tform = fitgeotrans([0, 0; 1, 1; -1, 1], [-driftx(k), -drifty(k); 1-driftx(k), 1-drifty(k); -1-driftx(k), 1-drifty(k)],'affine');
        ra = imref2d(size(I));
        I2 = imwarp(I,tform,'OutputView',ra);
        %         I2 = I;
        I2 = im2uint16(I2);
        I3 = im2double(I2);
        %         figure(4)
        %         imshow(I3/max(max(I3)));
        %         pause(0.001);
        if (frames*xdim*ydim*2)<2^32 %old version, no bif tiff
            if k == 1
                 imwrite(I2,strcat(input_movie(1:end-4),'_driftcorrected.tif'),'tif');
            else
                 imwrite(I2,strcat(input_movie(1:end-4),'_driftcorrected.tif'),'tif','WriteMode','append');
             end
        else %new version, big tiff
            clear options;
            options.big = true; % Use BigTIFF format
            options.append = true;
            saveastiff(I2, strcat(input_movie(1:end-4),'_driftcorrected.tif'), options); % 4GB Big TIFF file
        end
    end
    
   
    
    if delOriginal==1
        delete (input_movie)
    end
 
    
    
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
    
    
    
    
end
  %% Update Status:
            statement='DONE';
            set(handles.statusbox,'String',statement);
            drawnow;pause(0.1);
            guidata(hObject, handles);




function maxDrift_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function maxDrift_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smoothFrame_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function smoothFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delOriginal.
function delOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to delOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of delOriginal



%% Dependent Functions
function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    output=[error,diffphase];
        
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc); 
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    md2 = fix(m/2); 
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
% Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
  
    % Compute crosscorrelation and locate the peak 
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid 
        [max1,loc1] = max(CC);   
        [max2,loc2] = max(max1); 
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end  

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return
%
function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return
%
function res = saveastiff(data, path, options)
% options.color
%   : true or FALSE
%   : If this is true, third dimension should be 3 and the data is saved as a color image.
% options.compress
%   : 'no', 'lzw', 'jpeg' or 'adobe'.
%     Compression type.
%       'no'    : Uncompressed(Default)
%       'lzw'   : lossless LZW
%       'jpeg'  : lossy JPEG (When using JPEG compression, ImageWidth,
%                 ImageLength, and RowsPerStrip must be multiples of 16.)
%       'adobe' : lossless Adobe-style
% options.message
%   : TRUE or false.
%     If this is false, all messages are skipped. 
% options.append
%   : true or FALSE
%     If path is exist, the data is appended to an existing file.
%     If path is not exist, this options is ignored.
% options.overwrite
%   : true or FALSE
%     Overwrite to an existing file.
% options.big 
%   : true or FALSE, 
%     Use 64 bit addressing and allows for files > 4GB
% 
% Defalut value of 'options' is
%     options.color     = false;
%     options.compress  = 'no';
%     options.message   = true;
%     options.append    = false;
%     options.overwrite = false;
%     options.big       = false;
% 
% res : Return value. It is 0 when the function is finished with no error.
%       If an error is occured in the function, it will have a positive
%       number (error code).
%
% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

tStart = tic;
errcode = 0;
try
%% Init options parameter    
if nargin < 3 % Use default options
    options.color = false;
    options.compress = 'no';
    options.message = true;
    options.append = false;
    options.overwrite = false;
end
if ~isfield(options, 'message'),   options.message   = true; end
if ~isfield(options, 'append'),    options.append    = false; end
if ~isfield(options, 'compress'),  options.compress  = 'no';  end
if ~isfield(options, 'color'),     options.color     = false; end
if ~isfield(options, 'overwrite'), options.overwrite = false; end
if  isfield(options, 'big') == 0,  options.big       = false; end

if isempty(data), errcode = 1; assert(false); end
if (options.color == false && ndims(data) > 3) || ...
   (options.color == true && ndims(data) > 4)
    % Maximum dimension of a grayscale image is 3 of [height, width, frame]
    % Maximum dimension of a color image is 4 of [height, width, color, frame]
    errcode = 2; assert(false);
end

%% Get image informations
% http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
if ~options.color
    if ndims(data) >= 4, errcode = 2; assert(false); end;
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.Photometric = Tiff.Photometric.MinIsWhite;
%     tagstruct.Photometric = Tiff.Photometric.Mask;
%     tagstruct.Photometric = Tiff.Photometric.Separated;
else
    if ndims(data) >= 5, errcode = 2; assert(false); end;
    [height, width, cc, depth] = size(data); % cc: color channels. 3: rgb, 4: rgb with alpha channel
    if cc ~= 3 && cc ~= 4, errcode = 3; assert(false); end;
    tagstruct.Photometric = Tiff.Photometric.RGB;
%     tagstruct.Photometric = Tiff.Photometric.CIELab;
%     tagstruct.Photometric = Tiff.Photometric.ICCLab;
%     tagstruct.Photometric = Tiff.Photometric.ITULab;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.Palette;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.YCbCr;
end
tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % (RGB RGB,RGB RGB,RGB RGB), http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html
% (Unsupported)tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Separate; % (RRR RRR, GGG GGG, BBB BBB), % http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html

%% Complex number
% http://www.awaresystems.be/imaging/tiff/tifftags/samplesperpixel.html
if ~options.color && isreal(data) % Grayscale image with real numbers
    tagstruct.SamplesPerPixel = 1;
    data = reshape(data, height, width, 1, depth);
elseif ~options.color && ~isreal(data) % Grayscale image with complex numbers
    tagstruct.SamplesPerPixel = 2;
    data = reshape([real(data) imag(data)], height, width, 2, depth);
elseif options.color && isreal(data) % Color image with real numbers
    tagstruct.SamplesPerPixel = cc;
    if cc == 4
        tagstruct.ExtraSamples = Tiff.ExtraSamples.AssociatedAlpha; % The forth channel is alpha channel
    end
    data = reshape(data, height, width, cc, depth);
elseif options.color && ~isreal(data) % Color image with complex numbers
    tagstruct.SamplesPerPixel = cc * 2;
    if cc == 3
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 3); % 3(real)+3(imag) = 6 = 3(rgb) + 3(Extra)
    else
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 5); % 4(real)+4(imag) = 8 = 3(rgb) + 5(Extra)
    end
    data = reshape([real(data) imag(data)], height, width, cc*2, depth);
end

%% Image compression
% http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
switch lower(options.compress)
    case 'no'
        tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
        tagstruct.Compression = Tiff.Compression.LZW;
    case 'jpeg'
        tagstruct.Compression = Tiff.Compression.JPEG;
    case 'adobe'
        tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise
        % Use tag nubmer in http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
        tagstruct.Compression = options.compress;
end

%% Sample format
% http://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
switch class(data)
    % Unsupported Matlab data type: char, logical, cell, struct, function_handle, class.
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        if options.color
            errcode = 4; assert(false);
        end
    case {'single', 'double', 'uint64', 'int64'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        % (Unsupported)Void, ComplexInt, ComplexIEEEFP
        errcode = 5; assert(false);
end

%% Bits per sample
% http://www.awaresystems.be/imaging/tiff/tifftags/bitspersample.html
switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tagstruct.BitsPerSample = 64;
    otherwise
        errcode = 5; assert(false);
end

%% Rows per strip
maxstripsize = 8*1024;
tagstruct.RowsPerStrip = ceil(maxstripsize/(width*(tagstruct.BitsPerSample/8)*size(data,3))); % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html
if tagstruct.Compression == Tiff.Compression.JPEG
    tagstruct.RowsPerStrip = max(16,round(tagstruct.RowsPerStrip/16)*16);
end

%% Overwrite check
if exist(path, 'file') && ~options.append
    if ~options.overwrite
        errcode = 6; assert(false);
    end
end

%% Save path configuration
path_parent = pwd;
[pathstr, fname, fext] = fileparts(path);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
    cd(pathstr);
end

%% Write image data to a file
file_opening_error_count = 0;
while ~exist('tfile', 'var')
    try
        if ~options.append % Make a new file
            s=whos('data');
            if s.bytes > 2^32-1 || options.big
                tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
            else
                tfile = Tiff([fname, fext], 'w');
            end
        else
            if ~exist([fname, fext], 'file') % Make a new file
                s=whos('data');
                if s.bytes > 2^32-1 || options.big
                    tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
                else
                    tfile = Tiff([fname, fext], 'w');
                end
            else % Append to an existing file
                tfile = Tiff([fname, fext], 'r+');
                while ~tfile.lastDirectory(); % Append a new image to the last directory of an exiting file
                    tfile.nextDirectory();
                end
                tfile.writeDirectory();
            end
        end
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                errcode = 7;
                assert(false);
            end
        end
    end
end

for d = 1:depth
    tfile.setTag(tagstruct);
    tfile.write(data(:, :, :, d));
    if d ~= depth
       tfile.writeDirectory();
    end
end

tfile.close();
if exist('path_parent', 'var'), cd(path_parent); end

%tElapsed = toc(tStart);
%if options.message
   % display(sprintf('The file was saved successfully. Elapsed time : %.3f s.', tElapsed));
%end

catch exception
%% Exception management
    if exist('tfile', 'var'), tfile.close(); end
    switch errcode
        case 1
            if options.message, error '''data'' is empty.'; end;
        case 2
            if options.message, error 'Data dimension is too large.'; end;
        case 3
            if options.message, error 'Third dimesion (color depth) should be 3 or 4.'; end;
        case 4
            if options.message, error 'Color image cannot have int8, int16 or int32 format.'; end;
        case 5
            if options.message, error 'Unsupported Matlab data type. (char, logical, cell, struct, function_handle, class)'; end;
        case 6
            if options.message, error 'File already exists.'; end;
        case 7
            if options.message, error(['Failed to open the file ''' path '''.']); end;
        otherwise
            if exist('fname', 'var') && exist('fext', 'var')
                delete([fname fext]);
            end
            if exist('path_parent', 'var'), cd(path_parent); end
            rethrow(exception);
    end
    if exist('path_parent', 'var'), cd(path_parent); end
end
res = errcode;



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



function maxRelShift_Callback(hObject, eventdata, handles)
% hObject    handle to maxRelShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxRelShift as text
%        str2double(get(hObject,'String')) returns contents of maxRelShift as a double


% --- Executes during object creation, after setting all properties.
function maxRelShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxRelShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
