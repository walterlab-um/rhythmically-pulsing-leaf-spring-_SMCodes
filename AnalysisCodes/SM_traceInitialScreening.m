function varargout = SM_traceInitialScreening(varargin)
% SM_TRACEINITIALSCREENING MATLAB code for SM_traceInitialScreening.fig
%      SM_TRACEINITIALSCREENING, by itself, creates a new SM_TRACEINITIALSCREENING or raises the existing
%      singleton*.
%
%      H = SM_TRACEINITIALSCREENING returns the handle to a new SM_TRACEINITIALSCREENING or the handle to
%      the existing singleton*.
%
%      SM_TRACEINITIALSCREENING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_TRACEINITIALSCREENING.M with the given input arguments.
%
%      SM_TRACEINITIALSCREENING('Property','Value',...) creates a new SM_TRACEINITIALSCREENING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_traceInitialScreening_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_traceInitialScreening_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_traceInitialScreening

% Last Modified by GUIDE v2.5 20-Oct-2017 17:13:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SM_traceInitialScreening_OpeningFcn, ...
                   'gui_OutputFcn',  @SM_traceInitialScreening_OutputFcn, ...
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


% --- Executes just before SM_traceInitialScreening is made visible.
function SM_traceInitialScreening_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_traceInitialScreening (see VARARGIN)

% Choose default command line output for SM_traceInitialScreening
handles.output = hObject;


handles.donor                                =[];
handles.acceptor                             =[];
handles.index                                = 1;
handles.timeUnit                             =0.1; % in sec 
handles.fName                                =[];
handles.pathName                             =[];
handles.workingDir                           =pwd; % it prints the path of program directory
handles.donorDup                             =[];
handles.acceptorDup                          =[];



% set axis properties for Intensities Plot
set(handles.Histogram,'XLim',[-0.2 1.2],'XGrid','on', 'YGrid','on', 'XTick',-0.2:0.2:1.2);
ylabel(handles.Histogram,'Normalized count (%)');
xlabel(handles.Histogram,'FRET Efficiency');


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_traceInitialScreening wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_traceInitialScreening_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function loadTrace_Callback(hObject, eventdata, handles)
workingDir=handles.workingDir; % where program is there
DefaultFilePath='C:\Users\Jagat\Dropbox\LAB_DATA_Analysis\Lab_DATA_Analysis';
[fileName,pathName]=uigetfile('*.traces','Select ".traces" file',DefaultFilePath);
if (fileName~=0)
    cd(pathName);
    handles.pathName=pathName;
    handles.fName=fileName;
    timeUnit=0.1;
    % read data
    fid=fopen(fileName,'r');     %reading mode
    len=fread(fid,1,'int32');    % len givs the number of data points in each file
    Ntraces=fread(fid,1,'int16');
    handles.NumberOfMolecule= Ntraces/2;
    handles.goodtraceList=handles.NumberOfMolecule;
    %mal=[num2str(Ntraces/2) ' molecules ' num2str(len) ' data points each ' folder ' ' fname];
    %disp(mal);
    handles.time=(0:(len-1))*timeUnit;% Total time is same for all molecules since movie stopped at the same time
    %raw=fread(fid,Ntraces*len,'int16');
    raw=fread(fid,(Ntraces+1)*len,'int16');
    fclose(fid);
    % peak table reading
    filename2 = strcat(fileName(1:size(fileName,2)-7), '.pks'); %subtract ".traces" add ".pks"
    peaktable = load (filename2);
    for m=1 :(Ntraces/2)
         handles.DonorX(m) = peaktable((2*m - 1), 2);
         handles.DonorY(m) = peaktable((2*m - 1), 3);
         handles.AcceptorX(m) = peaktable(2*m, 2);
         handles.AcceptorY(m) = peaktable(2*m, 3);
    end 
    
    index=(1:(Ntraces+1)*len); %index is a vector starting from 1 to Ntraces+1*len;
    Data=zeros(Ntraces+1,len);
    Data(index)=raw(index);
    min(min(Data))
    handles.AccDonor=[];
    handles.AccAcceptor=[];
    handles.molID=[];
   for i=1:(Ntraces/2)
        handles.donor{i}=Data(i*2,:);
        handles.donorDup{i}=handles.donor{i};
        handles.acceptor{i}=Data(i*2+1,:);
        handles.acceptorDup{i}=handles.acceptor{i};
    end
    cd (workingDir);
    
    guidata(hObject,handles)
    handles.donorH=handles.donor;
    handles.acceptorH=handles.acceptor;
    handles.molH=handles.NumberOfMolecule;
    guidata(hObject,handles)
    plotHistogram1(hObject, eventdata, handles)

    %% 9. Update Status:
statement=[num2str(handles.NumberOfMolecule) '  of  ' num2str(handles.NumberOfMolecule)];
set(handles.molNumberDisp,'String',statement);
set(handles.statementBox,'String','.traces file loaded successfully.');
drawnow;pause(0.1);
guidata(hObject, handles);
end



function plotHistogram1 (hObject, eventdata, handles)
HistStartFrame=str2double(get(handles.HistStartFrame,'string'));
HistEndFrame=str2double(get(handles.HistEndFrame,'string'));
handles.AccDonor=[];
handles.AccAcceptor=[];
for i=1:handles.molH
    handles.AccDonor=[handles.AccDonor handles.donorH{i}(HistStartFrame:HistEndFrame)];
    handles.AccAcceptor=[handles.AccAcceptor handles.acceptorH{i}(HistStartFrame:HistEndFrame)];
end
%% Make Histogram:

Fret=handles.AccAcceptor./(handles.AccAcceptor+handles.AccDonor);
[count,center]=hist(Fret,0:0.02:1);
perCount=100*count/sum(count);
perCount=perCount(2:end-1);
center=center(2:end-1);
axes(handles.Histogram);
bar(center,perCount);
xlim([-0 1])
set(handles.Histogram,'XLim',[-0 1],'XGrid','on', 'YGrid','on', 'XTick',-0.2:0.2:1.2);
ylabel(handles.Histogram,'Normalized count (%)');
xlabel(handles.Histogram,'FRET Efficiency');





function corrections_Callback(hObject, eventdata, handles)

if get(handles.GleakageCorrection,'value')==1
    DLeakage=str2double(get(handles.Gleakage,'string'))/100;
else
    DLeakage=0;
end

if get(handles.RleakageCorrection,'value')==1
    ALeakage=str2double(get(handles.Rleakage,'string'))/100;
else
    ALeakage=0;
end

if get(handles.GPower,'value')==1
    GPower=str2double(get(handles.GrPower,'string'));
    Gconst=402;
else
    GPower=0;
    Gconst=0;
end

if get(handles.RPower,'value')==1
    RPower=str2double(get(handles.RePower,'string'));
    Rconst=0;
else
    RPower=0;
    Rconst=0;
end

if get(handles.shift,'value')==1
    shiftSec=str2double(get(handles.shiftSec,'string'));
    timeResolution=str2double(get(handles.timeResolution,'string'));
    GreenFrameShift=round(shiftSec/(timeResolution/1000));
else
    GreenFrameShift=0;
end

if get(handles.backCorr,'value')==1
bgFrames=str2double(get(handles.bgFrames,'string'));
else
bgFrames=0;
end



for i=1: handles.NumberOfMolecule
   D=handles.donor{i};
   A=handles.acceptor{i};
   %% BG correction due to power
   cam1bg=GPower*0.92987+Gconst;
   cam2bg=RPower*0+ Rconst;
   D=D-cam1bg;
   A=A-cam2bg;
   %% Leakage correction improper dichroic.
   D1=D-A*ALeakage;
   A1=A-D*DLeakage;
   %% BG correction for each molecule
   if bgFrames>0
   Dsort=sort(D1);
   Asort=sort(A1);
   Dbg=mean(Dsort(1:bgFrames));
   Abg=mean(Asort(1:bgFrames));  
   D1=D1-Dbg;
   A1=A1-Abg;
   end
   
   %% Green signal relative shift.
   Dtemp=D1((1-GreenFrameShift):end);
   D2=[Dtemp mean(Dtemp(end+GreenFrameShift:end))*ones(1,abs(GreenFrameShift))];
   handles.donor{i}=D2;
   handles.acceptor{i}=A1;
end

  guidata(hObject,handles)
    handles.donorH=handles.donor;
    handles.acceptorH=handles.acceptor;
    handles.molH=handles.NumberOfMolecule;
    guidata(hObject,handles)

    plotHistogram1(hObject, eventdata, handles)


  %% 9. Update Status:
statement=[num2str(handles.NumberOfMolecule) '  of  ' num2str(handles.NumberOfMolecule)];
set(handles.molNumberDisp,'String',statement);
set(handles.statementBox,'String','Corrections are done.');
drawnow;pause(0.1);
guidata(hObject, handles);
guidata(hObject, handles);


function useScreening_Callback(hObject, eventdata, handles)
frameCheck=str2double(get(handles.frame,'string'));
if get(handles.tInt,'value')==1
    tIntVal=str2double(get(handles.tIntVal,'string'));
else
    tIntVal=0;
end
if get(handles.dInt,'value')==1
    dIntVal=str2double(get(handles.dIntVal,'string'));
else
    dIntVal=0;
end
if get(handles.aInt,'value')==1
    aIntVal=str2double(get(handles.aIntVal,'string'));
else
    aIntVal=0;
end
fretCheck= get(handles.Fret,'value');





if get(handles.dI,'value')==1
    DIval=str2double(get(handles.DIval,'string'));
else
    DIval=0;
end

if get(handles.aI,'value')==1
    AIval=str2double(get(handles.AIval,'string'));
else
    AIval=0;
end

frameLB=str2double(get(handles.frameLB,'string'));
frameUB=str2double(get(handles.frameUB,'string'));




goodtraceList=[];
for i=1: handles.NumberOfMolecule
    D=handles.donor{i};
    A=handles.acceptor{i};
    T=D+A;
    F=A./(D+A);
    
    %% Intinsity Based Criteria:
    sD=sort(D,'descend');
    sA=sort(A,'descend');
    sT=sort(T,'descend');
    sF=sort(F,'descend');
    %% Check for Condition 1
    if mean(sT(1:frameCheck))>tIntVal
        goodC1=1;
    else
        goodC1=0;
    end
    %% Check for Condition 2
    if mean(sD(1:frameCheck))>dIntVal
        goodC2=1;
    else
        goodC2=0;
    end
    %% Check for Condition 3
    if mean(sA(1:frameCheck))>aIntVal
        goodC3=1;
    else
        goodC3=0;
    end
     %% Check for Condition 4
    if mean(D(frameLB:frameUB))>DIval || mean(D(frameLB:frameUB))<0
        goodC4=1;
    else
        goodC4=0;
    end
    
      %% Check for Condition 5
    if mean(A(frameLB:frameUB))>AIval || mean(A(frameLB:frameUB))<0
        goodC5=1;
    else
        goodC5=0;
    end
    
    goodtrace=goodC1*goodC2*goodC3*goodC4*goodC5;
    goodtraceList=[goodtraceList;goodtrace];
end
%goodtracelist is containing ones for good tracesand zeros for excluded
%traces
handles.goodtraceList=goodtraceList;
guidata(hObject,handles)

   %% 9. Update Status:
statement=[num2str(length(goodtraceList(goodtraceList>0))) '  of  ' num2str(handles.NumberOfMolecule)];
set(handles.molNumberDisp,'String',statement);
set(handles.statementBox,'String','Selection Croteria Applied.');
drawnow;pause(0.1);
guidata(hObject, handles);

%% finally write into trace file again
j=1;
for i=1:handles.NumberOfMolecule
     D=handles.donor{i};
    A=handles.acceptor{i};
    DX=handles.DonorX(i);
    DY=handles.DonorY(i);
    AX=handles.AcceptorX(i);
    AY=handles.AcceptorY(i);
    if goodtraceList(i)==1
        handles.NewDonor{j}=D;
        handles.NewAcceptor{j}=A;
         handles.newDX(j)=DX;
         handles.newDY(j)=DY;
         handles.newAX(j)=AX;
        handles.newAY(j)=AY;
        j=j+1;
    end
end
handles.NewSize=j-1;



  guidata(hObject,handles)
    handles.donorH=handles.NewDonor;
    handles.acceptorH=handles.NewAcceptor;
    handles.molH=handles.NewSize;
     plotHistogram1(hObject, eventdata, handles)

    guidata(hObject,handles)





function saveTraces_Callback(hObject, eventdata, handles)


extention=get(handles.extention,'string');
mkdir([handles.pathName extention])
FileName=strcat(handles.pathName,extention,'\',handles.fName);



len = size(handles.NewDonor{1,1},2); %size of second dimension of traces
Ntraces = 2*handles.NewSize;
outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '.traces'];

index=(1:(Ntraces+1)*len); %index is a vector starting from 1 to Ntraces+1*len;
Data=zeros(Ntraces+1,len);
for i=1:(Ntraces/2)
    Data(i*2,:)=handles.NewDonor{i};
    Data(i*2+1,:)=handles.NewAcceptor{i};
end
raw(index)= Data(index);

fid = fopen(outfilename,'w','ieee-le');
%'l' or 'ieee-le' Little-endian BYTE ordering
fwrite(fid,len,'int32');
fwrite(fid,Ntraces,'int16');
fwrite(fid,raw,'int16');
fclose(fid);
%creating updated peaks table
peaktable = zeros(Ntraces,3);
for i = 1:Ntraces %new number of traces
    peaktable(i,1) = i; %pks file has just three columns, the first is number 1,2,3,,,,
    if mod(i,2) ~= 0 %mod(a,2) returns the remainder after division of a by 2, not zero - here are ODD numbers: acceptor
        peaktable(i,2:3) = cat(2,handles.newAX((i+1)/2),handles.newAY((i+1)/2));
    else % even numbers - donor
        peaktable(i,2:3) = cat(2, handles.newDX(i/2),handles.newDY(i/2));
    end
end


outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '.pks'];
save(outfilename, 'peaktable', '-ascii');


outfilename=[FileName(1:max(strfind(FileName,'.'))-1) extention '_img.tif'];
f1=strcat(handles.pathName,handles.fName);
Img= imread([f1(1:max(strfind(f1,'.'))-1) '_img.tif']);
imwrite(Img,outfilename)


   %% 9. Update Status:
set(handles.statementBox,'String','traces saved successfully.');
drawnow;pause(0.1);
guidata(hObject,handles)





function frameLB_Callback(hObject, eventdata, handles)
function frameUB_Callback(hObject, eventdata, handles)
function molNumberDisp_Callback(hObject, eventdata, handles)
function aI_Callback(hObject, eventdata, handles)
function AIval_Callback(hObject, eventdata, handles)
function dI_Callback(hObject, eventdata, handles)
function DIval_Callback(hObject, eventdata, handles)
function Fret_Callback(hObject, eventdata, handles)
function checkbox17_Callback(hObject, eventdata, handles)
function edit22_Callback(hObject, eventdata, handles)
function aInt_Callback(hObject, eventdata, handles)
function aIntVal_Callback(hObject, eventdata, handles)
function dInt_Callback(hObject, eventdata, handles)
function dIntVal_Callback(hObject, eventdata, handles)
function frame_Callback(hObject, eventdata, handles)
function tInt_Callback(hObject, eventdata, handles)
function tIntVal_Callback(hObject, eventdata, handles)
function timeResolution_Callback(hObject, eventdata, handles)
function shift_Callback(hObject, eventdata, handles)
function shiftSec_Callback(hObject, eventdata, handles)
function edit7_Callback(hObject, eventdata, handles)
function RPower_Callback(hObject, eventdata, handles)
function RePower_Callback(hObject, eventdata, handles)
function GmFactor_Callback(hObject, eventdata, handles)
function GPower_Callback(hObject, eventdata, handles)
function GrPower_Callback(hObject, eventdata, handles)
function extention_Callback(hObject, eventdata, handles)
function GleakageCorrection_Callback(hObject, eventdata, handles)
function Gleakage_Callback(hObject, eventdata, handles)
function RleakageCorrection_Callback(hObject, eventdata, handles)
function Rleakage_Callback(hObject, eventdata, handles)
function Rleakage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function GrPower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function GmFactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RePower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function shiftSec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function timeResolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tIntVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function frame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dIntVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function aIntVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DIval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function AIval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function frameLB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function frameUB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function molNumberDisp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function extention_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Gleakage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function statementBox_Callback(hObject, eventdata, handles)
% hObject    handle to statementBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statementBox as text
%        str2double(get(hObject,'String')) returns contents of statementBox as a double


% --- Executes during object creation, after setting all properties.
function statementBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statementBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in backCorr.
function backCorr_Callback(hObject, eventdata, handles)
% hObject    handle to backCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of backCorr



function bgFrames_Callback(hObject, eventdata, handles)
% hObject    handle to bgFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bgFrames as text
%        str2double(get(hObject,'String')) returns contents of bgFrames as a double


% --- Executes during object creation, after setting all properties.
function bgFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bgFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


