 function varargout = SM_Trace_Viewing_timecCalculations(varargin)
% SM_TRACE_VIEWING_TIMECCALCULATIONS MATLAB code for SM_Trace_Viewing_timecCalculations.fig
%      SM_TRACE_VIEWING_TIMECCALCULATIONS, by itself, creates a new SM_TRACE_VIEWING_TIMECCALCULATIONS or raises the existing
%      singleton*.
%
%      H = SM_TRACE_VIEWING_TIMECCALCULATIONS returns the handle to a new SM_TRACE_VIEWING_TIMECCALCULATIONS or the handle to
%      the existing singleton*.
%
%      SM_TRACE_VIEWING_TIMECCALCULATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_TRACE_VIEWING_TIMECCALCULATIONS.M with the given input arguments.
%
%      SM_TRACE_VIEWING_TIMECCALCULATIONS('Property','Value',...) creates a new SM_TRACE_VIEWING_TIMECCALCULATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_Trace_Viewing_timecCalculations_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_Trace_Viewing_timecCalculations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_Trace_Viewing_timecCalculations

% Last Modified by GUIDE v2.5 19-Sep-2022 14:28:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SM_Trace_Viewing_timecCalculations_OpeningFcn, ...
    'gui_OutputFcn',  @SM_Trace_Viewing_timecCalculations_OutputFcn, ...
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


% --- Executes just before SM_Trace_Viewing_timecCalculations is made visible.
function SM_Trace_Viewing_timecCalculations_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_Trace_Viewing_timecCalculations (see VARARGIN)

% Choose default command line output for SM_Trace_Viewing_timecCalculations
%

handles.output                               = hObject;
handles.time                                 =[];
handles.donor                                =[];
handles.acceptor                             =[];
handles.index                                = 1;
handles.NumberOfMolecule                     =[];
handles.ShowIntensitiesAndFret               =[];
handles.fName                                =[];
handles.fret                                 =[];
handles.pathName                             =[];
handles.X                                    =[];
handles.workingDir                           =pwd; % it prints the path of program directory
handles.backavgD                             =[];
handles.backavgA                             =[];
handles.donorDup                             =[];
handles.acceptorDup                          =[];
handles.minTotal                             =55;
handles.donorLowerCutoff                     =0;
handles.acceptorLowerCutoff                  =0;


%-----------------------------------------------------
%   inOrder to see molecules Info Such as helx_Molecule'xy'
handles.MoleculeInformation='Molecule Info';
set(handles.MoleculeInfo,'String',handles.MoleculeInformation)
%------------------------------------------------------
% set axis properties for Intensities Plot
set(handles.IntensitiesTrace,'YLim',[-50 700],'XGrid','on', 'YGrid','on');
ylabel(handles.IntensitiesTrace,'Intensity (a.u.)')

%handles.hdlzoom=zoom;
%setAllowAxesZoom(handles.hdlzoom,handles.IntensitiesTrace,true);
handles.f=gcf;

set(handles.IntensitiesTrace, 'ButtonDownFcn', {@IntensitiesTrace_ButtonDownFcn,handles});
%set(handles.IntensitiesTrace)

% This assigns a function handle to its ButtondownFcn property (created by using the @ symbol before the function name).
% set axis properties for Fret Plot
set(handles.FretTrace,'YLim',[-0.1 1.1],'XGrid','on', 'YGrid','on','YTick',-0.1:0.2:1.1);
xlabel(handles.FretTrace,'Time (s)');
ylabel(handles.FretTrace,'FRET Efficiency');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_Trace_Viewing_timecCalculations wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_Trace_Viewing_timecCalculations_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in loadFiles.
function loadFiles_Callback(hObject, eventdata, handles)
workingDir=handles.workingDir; % where program is there
[filelist, pathName] = uigetfile('*.dat','Select files on which you wish to perform time analysis.','MultiSelect','on');
if iscell(filelist) == 0
    filelist2{1} = filelist;
    filelist = filelist2;
    clear filelist2;
end
handles.pathName=pathName;
handles.filelist=filelist;
handles.NumberOfMolecule= length(filelist);
for i=1: length(filelist)
    handles.data{i}=load(strcat(pathName,'\',filelist{i}));    % read data
end   
 handles.DupData= handles.data;
    handles.index=1;
    guidata(hObject,handles) % pass information to funciton below
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotIntensitiesAndFret(hObject,eventdata,handles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    statement=['Molecule: ' filelist{handles.index}];
    set(handles.displayStatus,'String',statement);
    set(handles.jumpby,'String',num2str(handles.index));
    drawnow;pause(0.1);
    guidata(hObject,handles)
function Next_Callback(hObject, eventdata, handles)
handles.index=handles.index+1;
if handles.index>handles.NumberOfMolecule
statement='You are at the end of the movie. Can not go forward..';
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
handles.index=handles.NumberOfMolecule;
else
guidata(hObject,handles);
PlotIntensitiesAndFret(hObject, eventdata, handles)
 statement=['Molecule: ' handles.filelist{handles.index}];
 set(handles.displayStatus,'String',statement);
set(handles.jumpby,'String',num2str(handles.index));
drawnow;pause(0.1);
end
guidata(hObject,handles);
function Previous_Callback(hObject, eventdata, handles)
handles.index=handles.index-1;
if handles.index<1
    statement='You are at the bigining. Can not go back.';
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
handles.index=1;
else
guidata(hObject,handles);
PlotIntensitiesAndFret(hObject, eventdata, handles)
 statement=['Molecule: ' handles.filelist{handles.index}];
 set(handles.displayStatus,'String',statement);
set(handles.jumpby,'String',num2str(handles.index));
drawnow;pause(0.1);
end
guidata(hObject,handles);


function PlotIntensitiesAndFret(hObject, eventdata, handles)
data=handles.data{handles.index};
timeUB=max(data(:,1));
IntUB= max([data(:,3);data(:,2)]);
IntLB= min([data(:,3);data(:,2)]);
axes(handles.IntensitiesTrace)

if get(handles.isOnT,'value')==1
plot(data(:,1),data(:,2)+data(:,3),'Color',[0.8 0.8 0.8]);
hold on;
end
if get(handles.isOnR,'value')==1
plot(data(:,1),data(:,3),'r');
hold on;
end
if get(handles.isOnG,'value')==1
plot(data(:,1),data(:,2),'g');
end
hold off;
% set properties for 'IntensitiesTrace' axes everytime graph is plotted
set(handles.IntensitiesTrace,'XLim',[0 timeUB],'XGrid','on','XTickLabel',[]);
set(handles.IntensitiesTrace,'YLim',[IntLB-50 IntUB+50],'XGrid','on', 'YGrid','on');
ylabel(handles.IntensitiesTrace,'Intensity (a.u.)')
set(handles.IntensitiesTrace, 'ButtonDownFcn', {@IntensitiesTrace_ButtonDownFcn,handles});
h = zoom;
 % set(h,'Enable','on');
%   set(h,'Enable','off');
set(h,'ButtonDownFilter',@mycallback);
setAllowAxesZoom(h,handles.IntensitiesTrace,true);
setAxesZoomMotion(h,handles.IntensitiesTrace,'horizontal');
set(handles.IntensitiesTrace, 'ButtonDownFcn', {@IntensitiesTrace_ButtonDownFcn,handles});


% Plot Time Vs Fret
if get(handles.isOnF,'value')==1
plot(handles.FretTrace,data(:,1),data(:,4),'b');
end
set(handles.FretTrace,'XLim',[0 timeUB],'XGrid','on');
set(handles.FretTrace,'YLim',[-0.2 1.2],'YGrid','on','YTick',-0.1:0.2:1.1);
xlabel(handles.FretTrace,'Time (s)');
ylabel(handles.FretTrace,'FRET Efficiency');
%----------------------------------------
setAllowAxesZoom(h,handles.FretTrace,true);

setAxesZoomMotion(h,handles.FretTrace,'horizontal');
%--------------------------------------
linkaxes([handles.IntensitiesTrace handles.FretTrace],'x');
%--------------------------------------
% % To display Molecuse Info
% handles.fName;
% len=size(handles.fName,2);
% statement=[ handles.fName(1:len-7) 'Molecule'  num2str(handles.index)];% statement to display
% set(handles.MoleculeInfo,'String',statement)

guidata(hObject,handles)


function manBackSub_Callback(hObject, eventdata, handles)

dback=str2double(get(handles.DonorBackground,'String'));
aback=str2double(get(handles.AcceeptorBackground,'String'));
if isempty(dback)
    dback=0;
end
if isempty(aback)
    aback=0;
end
data=handles.data{handles.index};
data(:,2)=data(:,2)-dback;
data(:,3)=data(:,3)-aback;
data(:,4)=data(:,3)./(data(:,3)+data(:,2));
handles.data{handles.index}=data;

guidata(hObject,handles)
PlotIntensitiesAndFret(hObject, eventdata, handles)
statement=['Manual Background Subtration is done for Molecule # ' num2str(handles.index)];
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function SubtractBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubtractBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DonorBackground_Callback(hObject, eventdata, handles)
% hObject    handle to DonorBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DonorBackground as text
%        str2double(get(hObject,'String')) returns contents of DonorBackground as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DonorBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DonorBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AcceeptorBackground_Callback(hObject, eventdata, handles)
% hObject    handle to AcceeptorBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AcceeptorBackground as text
%        str2double(get(hObject,'String')) returns contents of AcceeptorBackground as a double
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function AcceeptorBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AcceeptorBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in feedToHistogramAndSave.
function feedToHistogramAndSave_Callback(hObject, eventdata, handles)
listA=get(handles.listA, 'value');
listB=get(handles.listB, 'value');
listC=get(handles.listC, 'value');

if listA==0 && listB==0 && listC==0
    f = msgbox("please select at least one time list","Error","error");
    E=1;
elseif  listA==0 && listC==0
    E=0;
    list=2;
elseif  listB==0 && listC==0
    E=0;
    list=1;
elseif  listA==0 && listB==0
    E=0;
    list=3;
else
    f = msgbox("please select only one time list","Error","error");
    E=1;
end

if E==0
    clear fretE X Y mousebutton templength ToothStartTime ToothStopTime RealToothStartTime RealToothStopTime RealSawtoothTimes;
    [handles.X,Y,mousebutton]=ginput;
    X=handles.X;
    if rem(length(X),2)==1
        msgbox("Please click even number of times","Error","error"); % check to see if you have clicked even number of times
    else
        if issorted(X)==0
            msgbox("Sequences are not in order!","Error","error"); % check to see if you have clicked in order
        else
            Xodd = X(1:2:end);  
            Xeven = X(2:2:end);  
            timeRecord=Xeven-Xodd;
        end
    end
end





%% Update list and save data
pathName=handles.pathName;
filelist=handles.filelist;
file=filelist{handles.index};

listBoxA = cellstr(get(handles.listboxA,'String'));
listBoxB = cellstr(get(handles.listboxB,'String'));
listBoxC = cellstr(get(handles.listboxC,'String'));
newData=cellstr(num2str(timeRecord));

    colA    = cell(length(X)/2,1);
    colA(:) = {file};
    data=[colA newData];

if list==1
    listBoxData = [listBoxA; newData];
    set(handles.listboxA,'String',listBoxData);
    fileName=strcat(pathName,'timeListA.xlsx');
    if isfile(fileName)
     oldData=readcell(fileName); 
     AllData=[oldData; data];
     writecell(AllData,fileName);
    else
       writecell(data,fileName);
    end
elseif list==2
    listBoxData = [listBoxB; newData];
    set(handles.listboxB,'String',listBoxData);
    
     fileName=strcat(pathName,'timeListB.xlsx');
    if isfile(fileName)
     oldData=readcell(fileName); 
     AllData=[oldData; data];
     writecell(AllData,fileName);
    else
       writecell(data,fileName);
    end
    
elseif list==3
    listBoxData = [listBoxC; newData];
    set(handles.listboxC,'String',listBoxData);
     fileName=strcat(pathName,'timeListC.xlsx');
    if isfile(fileName)
     oldData=readcell(fileName); 
     AllData=[oldData; data];
     writecell(AllData,fileName);
    else
       writecell(data,fileName);
    end
end
guidata(hObject,handles)




function jumpby_Callback(hObject, eventdata, handles)
desiredMolecule=str2double(get(handles.jumpby,'String'));
if (desiredMolecule<=handles.NumberOfMolecule && desiredMolecule>0)
    handles.index  = desiredMolecule;
    guidata(hObject,handles);
    
    PlotIntensitiesAndFret(hObject, eventdata, handles)
    
    statement=['Molecule: ' handles.filelist{handles.index}];
    set(handles.displayStatus,'String',statement);
    drawnow;pause(0.1);
else
    statement='Molecule number out of range, Plz enter again';
    set(handles.displayStatus,'String',statement);
    drawnow;pause(0.1);
end

% --- Executes during object creation, after setting all properties.
function jumpby_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jumpby (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on mouse press over axes background.
function IntensitiesTrace_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IntensitiesTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[handles.dbackAvg,handles.abackAvg]=ginput(2);
%handles.xy_g_input=xy_g_input;
%handles.IntensitiesTrace=hObject;
%handles.output=hObject;
%handles = handles.IntensitiesTrace;

%handles
%zoom on
clear X Y;
[X,Y]=ginput(2);
statement=['Background Subtration in Progress for molecule ' num2str(handles.index) 'of' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - '  num2str(handles.DonorY(handles.index))];
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
backgroundrange=(floor(X(1)/handles.timeUnit):floor(X(2)/handles.timeUnit));
%a=handles.index
%a=handles.donor{handles.index}

handles.backavgD=mean(handles.donor{handles.index}(backgroundrange));
handles.backavgA=mean(handles.acceptor{handles.index}(backgroundrange));
handles.donor{handles.index}=handles.donor{handles.index}-handles.backavgD;
handles.acceptor{handles.index}=handles.acceptor{handles.index}-handles.backavgA;
guidata(hObject,handles)
PlotIntensitiesAndFret(hObject, eventdata, handles)

statement=['Background Subtration is done for Molecule # ' num2str(handles.index)];
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
%guidata(src,handles);
guidata(hObject,handles);

function [flag] = mycallback(obj,event_obj,handles)

%disp(['Clicked ' get(obj,'Type') ' object'])
%gcf=handles.f;
% obj gives figure or axes etc when clicked on GUI

objType = get(obj,'Type');% type could be figure or axes or ... depending on where clicked
%objTag = get(obj,'Tag');% it tells me Tag of where I clicked
if strcmp(objType,'axes')
    flag = true;
else
    flag = false; % zoom works if I click on lines
end
function histStatus_Callback(hObject, eventdata, handles)
function histStatus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MoleculeInfo_Callback(hObject, eventdata, handles)
function MoleculeInfo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ForTDP.
function ForTDP_Callback(hObject, eventdata, handles)
% hObject    handle to ForTDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ForTDP
%isChecked=get (hObject,'Value');
%pause(0.2);
%g=handles.index

% after taking Molecule for TDP, plot next molecule

guidata(hObject,handles)


function SaveTrace_Callback(hObject, eventdata, handles)

isChecked= get(handles.SaveTrace,'Value');


if (isChecked==1)
    
    %makeing directory for saving Excellent traces
    
    dirForSavingBestTrace='AwesomeTrace';
    % handles.paht comes with '\' at the end
    dest_folder=[handles.pathName,dirForSavingBestTrace];
    if exist (dest_folder,'dir');
    else mkdir(handles.pathName,dirForSavingBestTrace);
    end
    %-----------------------------------------------------
    cd (dest_folder)
    handles.fName;
    length=size(handles.fName,2);
    fNamePrefix=[handles.fName(1:length-7) 'Molecule'  num2str(handles.index)];
    % -7 for '.trace'
    %----------------------------------------------------------------------
    fret=handles.acceptor{handles.index}./(handles.acceptor{handles.index}+handles.donor {handles.index});
    %----------------------------------------------------------------------
    dataToSave=[handles.time',handles.donor{handles.index}',handles.acceptor{handles.index}',fret'];
    fileName_Data=[fNamePrefix '.dat'];
    save(fileName_Data,'dataToSave','-ascii')
    
    
    
    timeUB=max(handles.time)+handles.timeUnit;
    IntUB= max(handles.acceptor{handles.index}+handles.donor{handles.index});
    hdl=figure('Visible','off');
    hdl_s1=subplot(211);
    plot(hdl_s1,handles.time,handles.donor{handles.index},'g',handles.time,handles.acceptor{handles.index},'r');
    set(hdl_s1,'YLim',[-50 IntUB]);
    set(hdl_s1,'XLim',[0 timeUB]);
    ylabel(hdl_s1,'Intensity (a.u.)');
    
    %---------------------------------------
    hdl_s2=subplot(212);
    
    plot(hdl_s2,handles.time,fret,'b');
    set(hdl_s2,'YLim',[-0.1 1.1],'YTick',-0.1:0.2:1.1);
    xlabel(hdl_s2,'Time (s)');
    ylabel(hdl_s2,'FRET Efficiency');
    %----------------------------------------
    
    fName_Picture=[fNamePrefix '.png'];
    saveas(hdl,fName_Picture,'png')
    %-----------------------------------------
    set(handles.SaveTrace,'Value',0)
    % after everything done return to working directory
    statement=['Trace and Data file saved for ' 'Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - '  num2str(handles.DonorY(handles.index))];
    set(handles.displayStatus,'String',statement);
    cd(handles.workingDir);
end

function GetBackTrace_Callback(hObject, eventdata, handles)
handles.data{handles.index}=handles.DupData{handles.index};
PlotIntensitiesAndFret(hObject, eventdata, handles);
statement=['Reset Molecule: ' handles.filelist{handles.index}];
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
guidata(hObject,handles)

function GetAverageFret_Callback(hObject, eventdata, handles)
% hObject    handle to GetAverageFret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear X Y;
[X,Y]=ginput(2);
initialIndex=floor(X(1)/handles.timeUnit);
finalIndex=floor(X(2)/handles.timeUnit);
rangeOfX=(initialIndex:finalIndex);
%--------------------------------------------------------------------------
% Calculat Fret Below
len=finalIndex-initialIndex;
%Y=handles.time' % handles.time is row vector
fretE=zeros(len,1);
donor=zeros(len,1);
acceptor=zeros(len,1);

for m=initialIndex:finalIndex
    if (handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m))==0
        fretE(m)=NaN; % This is to avoid undefined fretE interfering with future analysis
        %elseif (acceptor(i,m)<2*AcceptorTreshold) & (donor(i,m)<2*DonorTreshold)
        % fretE(m)=NaN;
    elseif handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m) <handles.minTotal %(acceptor(i,m)+donor(i,m) < 200 &  acceptor(i,m)+donor(i,m) > 1000)% minTotal
        fretE(m)=NaN;
    else
        fretE(m)=handles.acceptor{handles.index}(m)./(handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m));
        donor(m)=handles.donor {handles.index}(m);
        acceptor(m)=handles.acceptor{handles.index}(m);
    end
end
fretE=fretE';
donor=donor';
acceptor=acceptor';
%--------------------------------------------------------------------------
AverageFret=mean(fretE(rangeOfX));
Averagedonor=mean(donor(rangeOfX));
Averageacceptor=mean(acceptor(rangeOfX));

statement=['Average Donor= ' num2str(Averagedonor) ', Average Acceptor= ' num2str(Averageacceptor) ', Average FRET= ' num2str(AverageFret)];
set(handles.displayStatus,'String',statement);

function SmoothSignal_Callback(hObject, eventdata, handles)
[d,a] =slidingAvgForGUI(handles);
data=handles.data{handles.index};
data(:,2)=d;
data(:,3)=a;
handles.data{handles.index}=data;
guidata(hObject,handles);
PlotIntensitiesAndFret(hObject, eventdata, handles);

function sawtooth_fret_norm_GUI_WindowKeyPressFcn(hObject, eventdata, handles)


function [avgDonorSignal, avgAcceptorSignal]= slidingAvgForGUI(handles)
data=handles.data{handles.index};
d= data(:,2);     % d is for donor
a= data(:,3);  % a is for acceptor
N=2;              % N represnets size of window, in another word how many points to take for avg
if (isempty(d)) || (isempty(a)) || (N<=0)                                              % If the input array is empty or N is non-positive,
    disp(sprintf('SlidingAvg: (Error) donor/accepor data empty or N null.'));     % an error is reported to the standard output and the
    return;                                                               % execution of the routine is stopped.
end 
if (N==1)                                                              % If the number of neighbouring points over which the sliding
    avgDonorSignal     = d;                                               % average will be performed is '1', then no average actually occur and
    avgAcceptorSignal  = a;
    return;                                                               % OUTPUT_ARRAY will be the copy of INPUT_ARRAY and the execution of the routine
end                                                               % is stopped.
nx= length(d);                                                  % The length of the input data structure is acquired to later evaluate the 'mean' over the appropriate boundaries.
% length of donor is chosen since a and d are of same length.
if (N>=(2*(nx-1)))                                                     % If the number of neighbouring points over which the sliding
    avgDonorSignal      = mean(d)*ones(size(d));                                        % average will be performed is large enough, then the average actually covers all the points
    avgAcceptorSignal   = mean(a)*ones(size(a));
    return;                                                               % of INPUT_ARRAY, for each index of OUTPUT_ARRAY and some CPU time can be gained by such an approach.
end                                                                % The execution of the routine is stopped.
avgDonorSignal      = zeros(size(d));         % In all the other situations, the initialization of the output data structure is performed.
avgAcceptorSignal   = zeros(size(a));
if rem(N,2)~=1                 % When N is even, then we proceed in taking the half of it:
    m = N/2;                      % m = N     /  2.
else                           % Otherwise (N >= 3, N odd), N-1 is even ( N-1 >= 2) and we proceed taking the half of it:
    m = (N-1)/2;                  % m = (N-1) /  2.
end 
% nx=length of array
for i=1:nx                                                % For each element (i-th) contained in the input numerical array, a check must be performed:
    if ((i-m) < 1) && ((i+m) <= nx)                             % If not enough points are available on the left of the i-th element..
        avgDonorSignal(i)    = mean(d(1:i+m));
        avgAcceptorSignal(i) = mean(a(1:i+m));                       % then we proceed to evaluate the mean from the first element to the (i + m)-th.
    elseif ((i-m) >= 1) && ((i+m) <= nx)                        % If enough points are available on the left and on the right of the i-th element..
        avgDonorSignal(i)    = mean(d(i-m:i+m));
        avgAcceptorSignal(i) = mean(a(i-m:i+m));                     % then we proceed to evaluate the mean on 2*m elements centered on the i-th position.
    elseif ((i-m) >= 1) && ((i+m) > nx)                          % If not enough points are available on the rigth of the i-th element..
        avgDonorSignal(i)    = mean(d(i-m:nx));
        avgAcceptorSignal(i) = mean(a(i-m:nx)); % then we proceed to evaluate the mean from the element (i - m)-th to the last one.
    elseif ((i-m) < 1) && ((i+m) > nx)                          % If not enough points are available on the left and on the rigth of the i-th element..
        avgDonorSignal(i)     = mean(d(1:nx));
        avgAcceptorSignal(i)  = mean(a(1:nx));% then we proceed to evaluate the mean from the first element to the last.
    end 
end


% --- Executes on button press in isOnT.
function isOnT_Callback(hObject, eventdata, handles)
% hObject    handle to isOnT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isOnT


% --- Executes on button press in isOnG.
function isOnG_Callback(hObject, eventdata, handles)
% hObject    handle to isOnG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isOnG


% --- Executes on button press in isOnR.
function isOnR_Callback(hObject, eventdata, handles)
% hObject    handle to isOnR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isOnR


% --- Executes on button press in isOnF.
function isOnF_Callback(hObject, eventdata, handles)
% hObject    handle to isOnF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isOnF


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
 keyPressed = eventdata.Key;
 if strcmpi(keyPressed,'l')
     loadTrace_ClickedCallback(hObject, eventdata, handles);
 elseif strcmpi(keyPressed,'s')
     feedToHistogramAndSave_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'rightarrow')
     Next_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'leftarrow')
     Previous_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'r')
     GetBackTrace_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'m')
     manBackSub_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'w')
     SmoothSignal_Callback(hObject, eventdata, handles)
 elseif strcmpi(keyPressed,'t')
         set(handles.SaveTrace,'value',1)
         pause(0.2)
     SaveTrace_Callback(hObject, eventdata, handles)   
 end
 


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)


% --- Executes on selection change in listboxA.
function listboxA_Callback(hObject, eventdata, handles)
% hObject    handle to listboxA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxA


% --- Executes during object creation, after setting all properties.
function listboxA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function displayStatus_Callback(hObject, eventdata, handles)
function displayStatus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxB.
function listboxB_Callback(hObject, eventdata, handles)
% hObject    handle to listboxB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxB


% --- Executes during object creation, after setting all properties.
function listboxB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxC.
function listboxC_Callback(hObject, eventdata, handles)
% hObject    handle to listboxC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxC


% --- Executes during object creation, after setting all properties.
function listboxC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in listC.
function listC_Callback(hObject, eventdata, handles)
