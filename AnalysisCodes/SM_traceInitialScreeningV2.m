function varargout = SM_traceInitialScreeningV2(varargin)
% SM_TRACEINITIALSCREENINGV2 MATLAB code for SM_traceInitialScreeningV2.fig
%      SM_TRACEINITIALSCREENINGV2, by itself, creates a new SM_TRACEINITIALSCREENINGV2 or raises the existing
%      singleton*.
%
%      H = SM_TRACEINITIALSCREENINGV2 returns the handle to a new SM_TRACEINITIALSCREENINGV2 or the handle to
%      the existing singleton*.
%
%      SM_TRACEINITIALSCREENINGV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_TRACEINITIALSCREENINGV2.M with the given input arguments.
%
%      SM_TRACEINITIALSCREENINGV2('Property','Value',...) creates a new SM_TRACEINITIALSCREENINGV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_traceInitialScreeningV2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_traceInitialScreeningV2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_traceInitialScreeningV2

% Last Modified by GUIDE v2.5 08-Sep-2022 14:19:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SM_traceInitialScreeningV2_OpeningFcn, ...
    'gui_OutputFcn',  @SM_traceInitialScreeningV2_OutputFcn, ...
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


% --- Executes just before SM_traceInitialScreeningV2 is made visible.
function SM_traceInitialScreeningV2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_traceInitialScreeningV2 (see VARARGIN)

% Choose default command line output for SM_traceInitialScreeningV2
handles.output = hObject;


handles.workingDir                           =pwd; % it prints the path of program directory


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_traceInitialScreeningV2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_traceInitialScreeningV2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadTrace.
function loadTrace_Callback(hObject, eventdata, handles)
workingDir='/Users/sujamac/Desktop/Best Analysis Codes As of 07062022/Files to run trial codes/trace files/movie - 1.ome.traces'; % where program is there
[fileName,pathName]=uigetfile('*.traces','Select ".traces" file',workingDir);
handles.pathName=pathName;
handles.fName=fileName;
guidata(hObject,handles)

handles=ReadDataFromFile(hObject, eventdata, handles);
guidata(hObject,handles)
CropData(hObject, eventdata, handles)
readPeaktable(hObject, eventdata, handles)
plotHistogram1(hObject, eventdata, handles)
plotHistogram2(hObject, eventdata, handles)
MolCountOutput=strcat('Showing:', num2str(handles.goodtraceCount),'/', num2str(handles.molCount),' molecules');
set(handles.MolCountOutput,'String', MolCountOutput);


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
handles=ReadDataFromFile(hObject, eventdata, handles);
guidata(hObject,handles)
CropData(hObject, eventdata, handles)
plotHistogram1(hObject, eventdata, handles)
plotHistogram2(hObject, eventdata, handles)
MolCountOutput=strcat('Showing:', num2str(handles.goodtraceCount),'/', num2str(handles.molCount),' molecules');
set(handles.MolCountOutput,'String', MolCountOutput);



function handles=ReadDataFromFile(hObject, eventdata, handles)
fullFileName=strcat(handles.pathName,handles.fName);
% read data
fid=fopen(fullFileName,'r');     %reading mode
length=fread(fid,1,'int32');     %length gives the number of frames
Ntraces=fread(fid,1,'int16');    % Ntraces/2 gives the number of molecules
handles.Nmol= Ntraces/2;
handles.goodMol=handles.Nmol;
rawData=fread(fid,(Ntraces+1)*length,'int16');
fclose(fid);
index=(1:(Ntraces+1)*length); %index is a vector starting from 1 to Ntraces+1*len;
Data=zeros(Ntraces+1,length);
Data(index)=rawData(index);

handles.donor=zeros(Ntraces/2,length);
handles.acceptor=zeros(Ntraces/2,length);

for i=1:(Ntraces/2)
    handles.donor(i,:)=Data(i*2,:);
    handles.acceptor(i,:)=Data(i*2+1,:);
end
handles.length=length;
handles.donorDup=handles.donor;
handles.acceptorDup=handles.acceptor;
handles.molCount=Ntraces/2;
handles.goodtraceCount=Ntraces/2;
guidata(hObject,handles)

function readPeaktable(hObject, eventdata, handles)
    fileName=fullfile(handles.pathName,handles.fName);
    Ntraces=2*handles.Nmol;
 % peak table reading
    filename2 = strcat(fileName(1:size(fileName,2)-7), '.pks'); %subtract ".traces" add ".pks"
    peaktable = load (filename2);
    for m=1 :(Ntraces/2)
         DonorX(m) = peaktable((2*m - 1), 2);
         DonorY(m) = peaktable((2*m - 1), 3);
         AcceptorX(m) = peaktable(2*m, 2);
         AcceptorY(m) = peaktable(2*m, 3);
    end 
    handles.peakTableMat=[DonorX' DonorY' AcceptorX' AcceptorY' ];
   guidata(hObject,handles)





function CropData(hObject, eventdata, handles)
HistStartFrame=str2double(get(handles.HistStartFrame,'string'));
HistEndFrame=str2double(get(handles.HistEndFrame,'string'));
handles.donor=handles.donorDup(:,HistStartFrame:HistEndFrame);
handles.acceptor=handles.acceptorDup(:,HistStartFrame:HistEndFrame);
guidata(hObject,handles)
function plotHistogram1(hObject, eventdata, handles)
subDonor=handles.donor;
subAcceptor=handles.acceptor;
%% Make Histogram:
FRET=subAcceptor./(subAcceptor+subDonor);
FRET = reshape(FRET',[size(FRET,1)*size(FRET,2),1]);
[count,center]=hist(FRET,0:0.02:1);
count=count/sum(count);
count=count(2:end-1);
center=center(2:end-1);
axes(handles.FRETHist);
bar(center,count);
xlim([-0 1])
set(handles.FRETHist,'XLim',[-0 1],'XGrid','on', 'YGrid','on', 'XTick',-0.2:0.2:1.2);
ylabel(handles.FRETHist,'Normalized count');
xlabel(handles.FRETHist,'FRET Efficiency');
function plotHistogram2(hObject, eventdata, handles)
subDonor=handles.donor;
subAcceptor=handles.acceptor;

%% Make Histogram:
subDonor=mean(subDonor,2);
subAcceptor=mean(subAcceptor,2);
Tot=(subAcceptor+subDonor);
% Plot Total Intensity
WhattoPlot=Tot;
Ax=handles.totInt;
c='k';
axes(Ax)
[count,center]=hist(WhattoPlot,min(WhattoPlot):(max(WhattoPlot)-min(WhattoPlot))/20:max(WhattoPlot));
count=count/sum(count);
count=count(1:end-1);
center=center(1:end-1);
bar(center,count,c);
set(Ax,'XLim',[min(WhattoPlot) max(WhattoPlot)],'XGrid','on', 'YGrid','on');
ylabel(Ax,'Normalized count');
xlabel(Ax,'Total Intensity');

set(handles.TLb,'String', round(min(WhattoPlot)));
set(handles.TUb,'String', round(max(WhattoPlot)));


% Plot Donor Intensity
WhattoPlot=subDonor;
Ax=handles.donorInt;
c='g';
axes(Ax)
[count,center]=hist(WhattoPlot,min(WhattoPlot):(max(WhattoPlot)-min(WhattoPlot))/20:max(WhattoPlot));
count=count/sum(count);
count=count(1:end-1);
center=center(1:end-1);
bar(center,count,c);
set(Ax,'XLim',[min(WhattoPlot) max(WhattoPlot)],'XGrid','on', 'YGrid','on');
xlabel(Ax,'Donor Intensity');
set(handles.DLb,'String', round(min(WhattoPlot)));
set(handles.DUb,'String', round(max(WhattoPlot)));

% Plot Acceptor Intensity
WhattoPlot=subAcceptor;
Ax=handles.acceptorInt;
c='r';
axes(Ax)
[count,center]=hist(WhattoPlot,min(WhattoPlot):(max(WhattoPlot)-min(WhattoPlot))/20:max(WhattoPlot));
count=count/sum(count);
count=count(1:end-1);
center=center(1:end-1);
bar(center,count,c);
set(Ax,'XLim',[min(WhattoPlot) max(WhattoPlot)],'XGrid','on', 'YGrid','on');
xlabel(Ax,'Acceptor Intensity');
set(handles.ALb,'String', round(min(WhattoPlot)));
set(handles.AUb,'String', round(max(WhattoPlot)));








% --- Executes on button press in useScreening.
function useScreening_Callback(hObject, eventdata, handles)

Tcb=get(handles.Tcb,'value');
Dcb=get(handles.Dcb,'value');
Acb=get(handles.Acb,'value');

%% First gather all user inputs
TLb= str2double(get(handles.TLb,'string'));
TUb= str2double(get(handles.TUb,'string'));
DLb= str2double(get(handles.DLb,'string'));
DUb= str2double(get(handles.DUb,'string'));
ALb= str2double(get(handles.ALb,'string'));
AUb= str2double(get(handles.AUb,'string'));

%Now find which traces are good
HistStartFrame=str2double(get(handles.HistStartFrame,'string'));
HistEndFrame=str2double(get(handles.HistEndFrame,'string'));
subDonor=handles.donorDup(:,HistStartFrame:HistEndFrame);
subAcceptor=handles.acceptorDup(:,HistStartFrame:HistEndFrame);
subDonor=mean(subDonor,2);
subAcceptor=mean(subAcceptor,2);
Tot=(subAcceptor+subDonor);

if Tcb==1
    C1=(Tot>TLb) & (Tot<TUb);
else
    C1=ones(size(Tot));
end



if Dcb==1
    C2=(subDonor>DLb) & (subDonor<DUb);
else
    C2=ones(size(Tot));
end

if Acb==1
    C3=(subAcceptor>ALb) & (subAcceptor<AUb);
else
    C3=ones(size(Tot));
end




goodtraceList=C1 & C2 & C3;


handles.donor=handles.donorDup(goodtraceList,HistStartFrame:HistEndFrame);
handles.acceptor=handles.acceptorDup(goodtraceList,HistStartFrame:HistEndFrame);
handles.goodtraceCount=length(handles.donor);

handles.goodtraceList=goodtraceList;
guidata(hObject,handles)

plotHistogram1(hObject, eventdata, handles)
plotHistogram2(hObject, eventdata, handles)

MolCountOutput=strcat('Showing:', num2str(handles.goodtraceCount),'/', num2str(handles.molCount),' molecules')
set(handles.MolCountOutput,'String', MolCountOutput);



% --- Executes on button press in performCorrection.
function performCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to performCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveNewTraces.
function saveNewTraces_Callback(hObject, eventdata, handles)
extention=get(handles.extention,'string');
mkdir([handles.pathName extention])
FileName=fullfile(handles.pathName,extention,handles.fName);
length = handles.length; %size of second dimension of traces
Ntraces = 2*handles.goodtraceCount;
outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '.traces'];
index=(1:(Ntraces+1)*length); %index is a vector starting from 1 to Ntraces+1*len;
Data=zeros(Ntraces+1,length);
newDonor=handles.donorDup(handles.goodtraceList,:);
newAcceptor=handles.acceptorDup(handles.goodtraceList,:);
for i=1:(Ntraces/2)
    Data(i*2,:)=newDonor(i,:);
    Data(i*2+1,:)=newAcceptor(i,:);
end
raw(index)= Data(index);


%% Saving .traces file

fid = fopen(outfilename,'w','ieee-le');
%'l' or 'ieee-le' Little-endian BYTE ordering
fwrite(fid,length,'int32');
fwrite(fid,Ntraces,'int16');
fwrite(fid,raw,'int16');
fclose(fid);

%% Saving .pks file
NewpeakTableMat=handles.peakTableMat(handles.goodtraceList,:);
DonorX=NewpeakTableMat(:,1)';
DonorY=NewpeakTableMat(:,2)';
AcceptorX=NewpeakTableMat(:,3)';
AcceptorY=NewpeakTableMat(:,4)';
peaktable=zeros(Ntraces,3);
for m=1 :(Ntraces/2)
    peaktable(m,1) = m; %pks file has just three columns, the first is number 1,2,3,,,,
    peaktable((2*m - 1), 2)=DonorX(m);
    peaktable((2*m - 1), 3)=DonorY(m);
    peaktable(2*m, 2)=AcceptorX(m);
    peaktable(2*m, 3)= AcceptorY(m);
end
outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '.pks'];
save(outfilename, 'peaktable', '-ascii');

%% Saving _img.tif file
outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '_img.tif'];
f1=fullfile(handles.pathName,handles.fName);
Img= imread([f1(1:max(strfind(f1,'.'))-1) '_img.tif']);
imwrite(Img,outfilename)









function HistStartFrame_Callback(hObject, eventdata, handles)
% hObject    handle to HistStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HistStartFrame as text
%        str2double(get(hObject,'String')) returns contents of HistStartFrame as a double


% --- Executes during object creation, after setting all properties.
function HistStartFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HistStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HistEndFrame_Callback(hObject, eventdata, handles)
% hObject    handle to HistEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HistEndFrame as text
%        str2double(get(hObject,'String')) returns contents of HistEndFrame as a double


% --- Executes during object creation, after setting all properties.
function HistEndFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HistEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TLb_Callback(hObject, eventdata, handles)
% hObject    handle to TLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TLb as text
%        str2double(get(hObject,'String')) returns contents of TLb as a double


% --- Executes during object creation, after setting all properties.
function TLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TUb_Callback(hObject, eventdata, handles)
% hObject    handle to TUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TUb as text
%        str2double(get(hObject,'String')) returns contents of TUb as a double


% --- Executes during object creation, after setting all properties.
function TUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DLb_Callback(hObject, eventdata, handles)
% hObject    handle to DLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DLb as text
%        str2double(get(hObject,'String')) returns contents of DLb as a double


% --- Executes during object creation, after setting all properties.
function DLb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DLb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DUb_Callback(hObject, eventdata, handles)
% hObject    handle to DUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DUb as text
%        str2double(get(hObject,'String')) returns contents of DUb as a double


% --- Executes during object creation, after setting all properties.
function DUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ALb_Callback(hObject, eventdata, handles)
% hObject    handle to ALb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ALb as text
%        str2double(get(hObject,'String')) returns contents of ALb as a double


% --- Executes during object creation, after setting all properties.
function ALb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ALb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AUb_Callback(hObject, eventdata, handles)
% hObject    handle to AUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AUb as text
%        str2double(get(hObject,'String')) returns contents of AUb as a double


% --- Executes during object creation, after setting all properties.
function AUb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AUb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tcb.
function Tcb_Callback(hObject, eventdata, handles)
% hObject    handle to Tcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tcb


% --- Executes on button press in Dcb.
function Dcb_Callback(hObject, eventdata, handles)
% hObject    handle to Dcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Dcb


% --- Executes on button press in Acb.
function Acb_Callback(hObject, eventdata, handles)
% hObject    handle to Acb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Acb



function extention_Callback(hObject, eventdata, handles)
% hObject    handle to extention (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extention as text
%        str2double(get(hObject,'String')) returns contents of extention as a double


% --- Executes during object creation, after setting all properties.
function extention_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extention (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
