 function varargout = SM_Trace_Viewing(varargin)
% SM_TRACE_VIEWING MATLAB code for SM_Trace_Viewing.fig
%      SM_TRACE_VIEWING, by itself, creates a new SM_TRACE_VIEWING or raises the existing
%      singleton*.
%
%      H = SM_TRACE_VIEWING returns the handle to a new SM_TRACE_VIEWING or the handle to
%      the existing singleton*.
%
%      SM_TRACE_VIEWING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_TRACE_VIEWING.M with the given input arguments.
%
%      SM_TRACE_VIEWING('Property','Value',...) creates a new SM_TRACE_VIEWING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_Trace_Viewing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_Trace_Viewing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_Trace_Viewing

% Last Modified by GUIDE v2.5 06-Jul-2022 17:36:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SM_Trace_Viewing_OpeningFcn, ...
    'gui_OutputFcn',  @SM_Trace_Viewing_OutputFcn, ...
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


% --- Executes just before SM_Trace_Viewing is made visible.
function SM_Trace_Viewing_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_Trace_Viewing (see VARARGIN)

% Choose default command line output for SM_Trace_Viewing
%

handles.output                               = hObject;
handles.time                                 =[];
handles.donor                                =[];
handles.acceptor                             =[];

handles.index                                = 1;
handles.NumberOfMolecule                     =[];
handles.ShowIntensitiesAndFret               =[];
handles.timeUnit                             =str2double(get(handles.timeUnit1,'string'))/1000;
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

handles.title         ='Normalized cumulative FRET-histogram';
set(handles.histStatus,'String',handles.title);
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

% setAllowAxesZoom(handles.h1,handles.FretTrace,false);
set(handles.Histogram,'XLim',[-0.2 1.2],'XTick',-0.2:0.2:1.2);
%set(handles.Histogram,'XLim',[-0.2 1.2],'XTick',0:0.2:45);

xlabel(handles.Histogram,'FRET Efficiency');
ylabel(handles.Histogram,'NormCount(%)');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_Trace_Viewing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_Trace_Viewing_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function loadTrace_ClickedCallback(hObject, eventdata, handles)
workingDir=handles.workingDir; % where program is there
DefaultFilePath='C:\Users\Jagat\Dropbox\LAB_DATA_Analysis\Lab_DATA_Analysis';
[fileName,pathName]=uigetfile('*.traces','Select ".traces" file',DefaultFilePath);
handles.fileName=fileName;
handles.pathName=pathName;
handles.timeUnit=str2double(get(handles.timeUnit1,'string'))/1000;
guidata(hObject,handles)
[~,~,ext] = fileparts(fileName);
comp=strcmp(ext,'.traces');
if (fileName~=0)
    if (comp==1) 
    cd(pathName);
    handles.pathName=pathName;
    handles.fName=fileName;
       % read data
    fid=fopen(fileName,'r');     %reading mode
    len=fread(fid,1,'int32');    % len givs the number of data points in each file
    Ntraces=fread(fid,1,'int16');
    handles.NumberOfMolecule= Ntraces/2;
    handles.time=(0:(len-1))*handles.timeUnit;% Total time is same for all molecules since movie stopped at the same time
    %raw=fread(fid,Ntraces*len,'int16');
    raw=fread(fid,(Ntraces+1)*len,'int16');
    
    % raw=abs(raw); %
    fclose(fid);
    % peak table reading
    filename2 = strcat(fileName(1:size(fileName,2)-7), '.pks');
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
    for i=1:(Ntraces/2)
        handles.donor{i}=Data(i*2,:);
        handles.acceptor{i}=Data(i*2+1,:);
        handles.donorDup{i}=handles.donor{i};
        handles.acceptorDup{i}=handles.acceptor{i};
    end
    
    cd (workingDir);
    % display intensity profiles of first molecule below
    handles.index=1;
    
    guidata(hObject,handles) % pass information to funciton below
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PlotIntensitiesAndFret(hObject,eventdata,handles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    statement=['Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - ' num2str(handles.DonorY(handles.index))];
    set(handles.displayStatus,'String',statement);
    set(handles.jumpby,'String',num2str(handles.index));
    drawnow;pause(0.1);
    guidata(hObject,handles)
    end  
 end

function displayStatus_Callback(hObject, eventdata, handles)
function displayStatus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
statement=['Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - ' num2str(handles.DonorY(handles.index))];
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
statement=['Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - ' num2str(handles.DonorY(handles.index))];
set(handles.displayStatus,'String',statement);
set(handles.jumpby,'String',num2str(handles.index));
drawnow;pause(0.1);
end
guidata(hObject,handles);


function PlotIntensitiesAndFret(hObject, eventdata, handles)
timeUB=max(handles.time)+handles.timeUnit;
IntUB= max(handles.acceptor{handles.index}+handles.donor{handles.index});
IntLB= min(0,min(min(handles.acceptor{handles.index}),min(handles.donor{handles.index})));
axes(handles.IntensitiesTrace)

if get(handles.isOnT,'value')==1
plot(handles.time,handles.acceptor{handles.index}+handles.donor{handles.index},'Color',[0.8 0.8 0.8]);
hold on;
end
if get(handles.isOnR,'value')==1
plot(handles.time,handles.acceptor{handles.index},'r');
hold on;
end
if get(handles.isOnG,'value')==1
plot(handles.time,handles.donor{handles.index},'g');
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

len=length(handles.time);
%fret=zeros(len,1);
 for m=1:len
%     if (handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m))==0
%         fret(m)=NaN; 
%     elseif handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m) <handles.minTotal %(acceptor(i,m)+donor(i,m) < 200 &  acceptor(i,m)+donor(i,m) > 1000)% minTotal
%         fret(m)=NaN;
%     else
        fret(m)=handles.acceptor{handles.index}(m)./(handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m));
%     end
 end


% Plot Time Vs Fret
if get(handles.isOnF,'value')==1
plot(handles.FretTrace,handles.time',fret,'b');
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
% To display Molecuse Info
handles.fName;
len=size(handles.fName,2);
statement=[ handles.fName(1:len-7) 'Molecule'  num2str(handles.index)];% statement to display
set(handles.MoleculeInfo,'String',statement)

guidata(hObject,handles)


function manBackSub_Callback(hObject, eventdata, handles)

statement='Background Subtration in Progress';
set(handles.displayStatus,'String',statement);
drawnow;pause(0.1);
dback        =str2double(get(handles.DonorBackground,'String'));
aback        =str2double(get(handles.AcceeptorBackground,'String'));
if isempty(dback)
    dback=0;
end
if isempty(aback)
    aback=0;
end
handles.donor{handles.index}=handles.donor{handles.index}-dback;
handles.acceptor{handles.index}=handles.acceptor{handles.index}-aback;
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
% hObject    handle to feedToHistogramAndSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear fretE X Y mousebutton templength ToothStartTime ToothStopTime RealToothStartTime RealToothStopTime RealSawtoothTimes;
[handles.X,Y,mousebutton]=ginput;
X=handles.X;
handles.fName;
len=size(handles.fName,2);
fname=handles.fName(1:len-7);% -7 for '.trace'

%----------------------------------------------------------------------
len=length(handles.time);
%Y=handles.time' % handles.time is row vector
fretE=zeros(len,1);
for m=1:len,
    if (handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m))==0
        fretE(m)=NaN; % This is to avoid undefined fretE interfering with future analysis
        %elseif (acceptor(i,m)<2*AcceptorTreshold) & (donor(i,m)<2*DonorTreshold)
        % fretE(m)=NaN;
    elseif handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m) <handles.minTotal %(acceptor(i,m)+donor(i,m) < 200 &  acceptor(i,m)+donor(i,m) > 1000)% minTotal
        fretE(m)=NaN;
    else
        fretE(m)=handles.acceptor{handles.index}(m)./(handles.acceptor{handles.index}(m)+handles.donor {handles.index}(m));
    end
end
fretE=fretE';
%----------------------------------------------------------------------
%fretE=handles.acceptor{handles.index}./(handles.acceptor{handles.index}+handles.donor {handles.index});

templength=length(X);

if (mod(templength,2)==0) && (templength > 0) % this condition guarentees that you clicked left and right equal number of times
    
    StartTimeCounter =0;
    StopTimeCounter  =0;
    
    for zz=1:templength
        
        if mousebutton(zz)==1
            % when mousebutton(zz) is 1, it means that the corresponding
            % time (X) is the starting point for the sawtooth
            StartTimeCounter = StartTimeCounter + 1;
            ToothStartTime(StartTimeCounter) = round(X(zz)/handles.timeUnit);
            RealToothStartTime(StartTimeCounter)=X(zz);
        else       % when mousebutton(zz) is 3, end time for sawtooth
            StopTimeCounter = StopTimeCounter + 1;
            ToothStopTime(StopTimeCounter) = round(X(zz)/handles.timeUnit);
            RealToothStopTime(StopTimeCounter)=X(zz);
        end
        
    end
    cd(handles.pathName);
    
    if StartTimeCounter==StopTimeCounter
        tooth_no = 0 ;
        Scale = 100 ./ ( ToothStopTime - ToothStartTime) ; %scaling factor array for the sawtooths in this trace
        %---------------------------------------------------------------
        RealSawtoothTimes= [ RealToothStopTime - RealToothStartTime ]';
        outname=[ fname ' timeinterval' num2str(handles.index) '_00.dat' ];
        if exist(outname)==2;
            PreviousSawtoothTimes=load(outname,'-ascii');
            RealSawtoothTimes=[PreviousSawtoothTimes' RealSawtoothTimes']' ;
        end
        output= RealSawtoothTimes ;
        save(outname, 'output', '-ascii');
        %-----------------------------------------------------------------
        %-----------------------------------------------------------------
        
        %-------------------------------------------------------------
        for kk = 1:StopTimeCounter % Number of sawtooth patterns from this molecule is same as the last value of StopTimeCounter
            %clear ToothProfile ToothProfileforHammy ToothWNumOfPoints;
            tooth_no = tooth_no + 1 ;
            temptime = [ ToothStartTime(kk):ToothStopTime(kk) ] - ToothStartTime(kk);
            ScaledTime = round( Scale(kk)*temptime );
            %index=[ToothStartTime(kk):ToothStopTime(kk)]
            tempfretE = fretE( ToothStartTime(kk):ToothStopTime(kk) );
            NaNTester=isnan(tempfretE);
            LengthOfScaledTime=length(ScaledTime);
            clear ToothProfile ToothProfileforvBFret ToothWNumOfPoints;
            % REALLY IMP
            %a=handles.donor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))';
            %b=handles.acceptor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))'
            %c=tempfretE; a and b are column vector while c
            % is row vector of same size
            
            %size(a),size(b),size(c)
            ToothProfile=[handles.time((ToothStartTime(kk):ToothStopTime(kk)))'  handles.donor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))' handles.acceptor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))' tempfretE' ];
            ToothProfileforvBFret=[handles.donor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))' handles.acceptor{handles.index}(ToothStartTime(kk):ToothStopTime(kk))'];
            %----------
            %---------------------------------------------------------------------
            outname = [fname '_' 'molecule-' num2str(handles.index) '_region-' num2str(tooth_no) '_00.dat'];
            save(outname,'ToothProfile','-ascii','-append') ;% here, each region profile that are selected after one 'w', are written into a separate file.
            
            outname = [fname '_' 'molecule-' num2str(handles.index) '_region-' num2str(tooth_no) '_00_forvBFRET.dat'];
            
            save(outname,'ToothProfileforvBFret','-ascii') ;
            %----------------------------------------------------------------------
            % Create one more variable combinedFret And All
            % segmetn of single trace will be treated as
            % singel molecule
            if kk==1
                combinedFret=tempfretE;
                combinedToothProfile=ToothProfile;
            else
                combinedFret=[combinedFret tempfretE];
                combinedToothProfile=[combinedToothProfile; ToothProfile]; % stack ToothProfile below combinedToothProfile
            end
        end
        
        
        topEdge      = 1.2;  % Max Fret Lim
        bottomEdge   = -0.2; % Min Fret Lim
        numBin       = 70;   % bin interval= 0.02   delta_freq=0.02;
        binEdges     = linspace(bottomEdge,topEdge,numBin+1);  %creating bins
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        grand_outname = [fname '_grand_00.dat'];  %%%% Notice the change from 'outname' to 'grand_outname'
        if exist(grand_outname)==2
            AllToothWNumOfPoints=load(grand_outname,'-ascii');  %%% Selected AllToothWNumOfPoints will  inside the
        else
            AllToothWNumOfPoints=[];
        end
        NumOfPoints  = size(combinedToothProfile,1);
        
        indices=find(combinedToothProfile(:,4)>-0.21 & combinedToothProfile(:,4)<=1.21); % fourth Column is fret
        
        NumOfValidPoints=length(indices); %This is a normalization factor, so it should be the total count of fret points that goes into the histogram. some fret values may fall outside our histogram range,
        % note that I changed the histogram range from 0:1 to -0.2:1.2
        NormalizedContr=zeros(NumOfPoints,1);
        NormalizedContr(indices)=1/NumOfValidPoints;
        
        combinedToothProfile=[combinedToothProfile NormalizedContr]; % size is consistent; now has five columns
        AllToothWNumOfPoints=[AllToothWNumOfPoints ; combinedToothProfile];
        save(grand_outname,'AllToothWNumOfPoints','-ascii');
        
        clear AllToothWNumOfPoints % after saving file, clear memory
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        grand_histogram_outname = [fname '_grand_histogram_00.dat'];
        moleculeFretBinCount    = [fname 'Single_histograms_00.dat'];
        if exist (moleculeFretBinCount)==2
            single_histogram   = load (moleculeFretBinCount);
        else
            single_histogram   = binEdges'; % assign as column
        end
        
        if exist(grand_histogram_outname)==2
            GrandHistogram=load(grand_histogram_outname);
        else
            GrandHistogram(:,1)=binEdges;
            GrandHistogram(:,2)=zeros(71,1);
        end
        %---------------------------------
        n         = histc(combinedFret,binEdges);%it gives count in each bin, any fret out of range discarted
        n         = n./sum(n);   % Normalize n w.r.t total valid points
        GrandHistogram(:,2)=GrandHistogram(:,2)+n';
        %---------------------------------
        single_histogram=[single_histogram n']; % to save each molecule fret contribution as column
        
        GrandHistogram=[binEdges' GrandHistogram(:,2)];
        
        %figure(4);
        bar(handles.Histogram,binEdges,GrandHistogram(:,2),'c');
        %---------------------------------------------------
        set(handles.Histogram,'XLim',[-0.2 1.2],'XTick',-0.2:0.2:1.2);
        xlabel(handles.Histogram,'FretEfficiency');
        ylabel(handles.Histogram,'NormCount(%)');
        %---------------------------------------------------
        %title 'Cumulative FRET Histogram Normalized with Dwell Time';
        TotalMoleculesForHist    =sum(GrandHistogram(:,2));
        histStatement=['Number of Molecule in Hist= ' num2str(TotalMoleculesForHist)];
        set(handles.histStatus,'String',histStatement);
        Normalization            =sum(n);
        save(grand_histogram_outname,'GrandHistogram','-ascii');
        save(moleculeFretBinCount,'single_histogram','-ascii')
        cd(handles.workingDir);
        statement=['Chk Histogram!!! Files Recorded for Molecule # ' num2str(handles.index)];
        set(handles.displayStatus,'String',statement);
        drawnow;pause(0.1);
        %---------------------------------------------------------------------------------------
        % For TDP and Saving Pictures
        isChecked=get(handles.ForTDP,'Value');
        if isChecked==1
            %--------------------------------------------
            newdir_forTDP='forTDP';   newdir_forSavingTrace='GoodTracePic';
            dest_folder=[handles.pathName,newdir_forTDP];
            dest_folder1=[handles.pathName,newdir_forSavingTrace];
            
            if exist (dest_folder,'dir');
            else mkdir(handles.pathName,newdir_forTDP);
            end
            if exist (dest_folder1,'dir');
            else mkdir(handles.pathName,newdir_forSavingTrace);
            end
            %---------------------------------------------
            handles.fName;
            len      =size(handles.fName,2);
            fname    =handles.fName(1:len-7);% -7 for removing '.trace'
            %---------------------------------------------
            timeDivider=handles.X;
            len=length(timeDivider);
            numOfFilesToWrite= len/2;
            %-------------------------------------------
            %save file below
            cd(dest_folder)
            for ii=1:numOfFilesToWrite,
                startTime=round(timeDivider(2*ii-1)/handles.timeUnit);
                endTime  =round(timeDivider(2*ii)/handles.timeUnit);
                dataToSaveForvBFRET=[handles.time(startTime:endTime)'  handles.donor{handles.index}(startTime:endTime)' handles.acceptor{handles.index}(startTime:endTime)'];
                
                outname = [fname '_molecule-' num2str(handles.index) '_region-' num2str(ii) '_forvBFRET.dat'];
                save (outname,'dataToSaveForvBFRET','-ascii');
            end
            % change folder and Save Pictures
            cd(dest_folder1)
            for ii=1:numOfFilesToWrite,
                startTime =round(timeDivider(2*ii-1)/handles.timeUnit);
                endTime   =round(timeDivider(2*ii)/handles.timeUnit);
                %dataToSaveForvBFRET=[handles.time(startTime:endTime)'  handles.donor{handles.index-1}(startTime:endTime)' handles.acceptor{handles.index-1}(startTime:endTime)'];
                
                hdl=figure('Visible','off');
                hdl_s1=subplot(211);
                plot(hdl_s1,handles.time(startTime:endTime),handles.donor{handles.index}(startTime:endTime),'g',handles.time(startTime:endTime),handles.acceptor{handles.index}(startTime:endTime),'r');
                a=max(handles.donor{handles.index}(startTime:endTime));
                b=max(handles.acceptor{handles.index}(startTime:endTime));
                if a>b
                    MaxYLim=a;
                else
                    MaxYLim=b;
                end
                
                set(hdl_s1,'YLim',[-100 MaxYLim+400]);
                ylabel(hdl_s1,'Intensity (a.u.)');
                
                %---------------------------------------
                hdl_s2=subplot(212);
                fret=handles.acceptor{handles.index}(startTime:endTime)./(handles.acceptor{handles.index}(startTime:endTime)+handles.donor{handles.index}(startTime:endTime));
                plot(hdl_s2,handles.time(startTime:endTime),fret,'b');
                set(hdl_s2,'YLim',[0 1],'YTick',0:0.2:1.0);
                xlabel(hdl_s2,'Time (s)');
                ylabel(hdl_s2,'FRET');
                %----------------------------------------
                outname = [fname '_molecule-' num2str(handles.index) '_region-' num2str(ii) '.png'];
                saveas(hdl,outname,'png')
                close(hdl)
                %-----------------------------------------
                
            end
            % after everything done return to working directory
            statement=['Data Saved For vBFret Program for ' 'Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule)];
            set(handles.displayStatus,'String',statement);
            % pause(1);
            set(handles.ForTDP,'Value',0)
            
            cd(handles.workingDir);
            
            
            % this line is important since if this invisible figure window is not closed zoom will be affected
            % in gui and program crashes
            %     handles.index   =handles.index+1;
            
            %     guidata(hObject,handles);
            %     PlotIntensitiesAndFret(hObject, eventdata, handles);
            
            
            
        end % End of Is checked
        %---------------------------------------------------------------------------------------
        
        handles.index   =handles.index+1;
        
        guidata(hObject,handles);
        PlotIntensitiesAndFret(hObject, eventdata, handles);
        
        statement=['Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - '  num2str(handles.DonorY(handles.index)) ];
        set(handles.displayStatus,'String',statement); 
        drawnow;pause(0.1);
    else
        
        statement='Plz Click the left and right mouse buttons same number of times';
        set(handles.displayStatus,'String',statement);
        drawnow;pause(0.1);
        cd(handles.workingDir);
        guidata(hObject,handles);
        
        PlotIntensitiesAndFret(hObject, eventdata, handles)
        
        
    end
    
else
    statement='You clicked odd number of times, please try again';
    set(handles.displayStatus,'String',statement);
    drawnow;pause(0.1);
    
    cd(handles.workingDir);
    guidata(hObject,handles);
    
    PlotIntensitiesAndFret(hObject, eventdata, handles)
end
guidata(hObject,handles)




function jumpby_Callback(hObject, eventdata, handles)
desiredMolecule=str2double(get(handles.jumpby,'String'));
if (desiredMolecule<=handles.NumberOfMolecule && desiredMolecule>0)
    handles.index  = desiredMolecule;
    guidata(hObject,handles);
    
    PlotIntensitiesAndFret(hObject, eventdata, handles)
    
    statement=['Jumped to Molecule # ' num2str(handles.index) ' of ' num2str(handles.NumberOfMolecule) ' at ' num2str(handles.DonorX(handles.index)) ' - '  num2str(handles.DonorY(handles.index))];
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
handles.donor{handles.index}=handles.donorDup{handles.index}; % Reset donor and
handles.acceptor{handles.index}=handles.acceptorDup{handles.index};% acceptor values here
guidata(hObject,handles);
PlotIntensitiesAndFret(hObject, eventdata, handles);
statement=['Reset Molecule: ' num2str(handles.index)];
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
handles.donor{handles.index}     =d; % 'd' and 'a' are averaged donor and averaged acceptor
handles.acceptor{handles.index}  =a; %
guidata(hObject,handles);
PlotIntensitiesAndFret(hObject, eventdata, handles);
function sawtooth_fret_norm_GUI_WindowKeyPressFcn(hObject, eventdata, handles)


function [avgDonorSignal, avgAcceptorSignal]= slidingAvgForGUI(handles)
d= handles.donor{handles.index} ;     % d is for donor
a= handles.acceptor{handles.index} ;  % a is for acceptor
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
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)


function showImage_ClickedCallback(hObject, eventdata, handles)
fileName=handles.fileName;
pathName=handles.pathName;
fileName2=strcat(fileName(1:end-7),'_img.tif');
I=imread(strcat(pathName,fileName2));
fileName3=strcat(fileName(1:end-7),'.pks');
peaktable = load (strcat(pathName,fileName3));
handles.index
figure(2);
imshow(I);
hold on;
plot(peaktable(1:2:end-1,2), peaktable(1:2:end-1,3), 'yo');
plot(peaktable(2:2:end,2)+size(I,2)/2, peaktable(2:2:end,3), 'co');

pause(1)
plot(peaktable(2*handles.index-1,2), peaktable(2*handles.index-1,3), 'go');
plot(peaktable(2*handles.index,2)+size(I,2)/2, peaktable(2*handles.index,3), 'ro');

hold off;



function timeUnit1_Callback(hObject, eventdata, handles)
% hObject    handle to timeUnit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timeUnit1 as text
%        str2double(get(hObject,'String')) returns contents of timeUnit1 as a double


% --- Executes during object creation, after setting all properties.
function timeUnit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeUnit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
