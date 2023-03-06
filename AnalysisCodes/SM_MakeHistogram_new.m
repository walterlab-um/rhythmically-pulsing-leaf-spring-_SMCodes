function varargout = SM_MakeHistogram_new(varargin)
% SM_MAKEHISTOGRAM_NEW MATLAB code for SM_MakeHistogram_new.fig
%      SM_MAKEHISTOGRAM_NEW, by itself, creates a new SM_MAKEHISTOGRAM_NEW or raises the existing
%      singleton*.
%
%      H = SM_MAKEHISTOGRAM_NEW returns the handle to a new SM_MAKEHISTOGRAM_NEW or the handle to
%      the existing singleton*.
%
%      SM_MAKEHISTOGRAM_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SM_MAKEHISTOGRAM_NEW.M with the given input arguments.
%
%      SM_MAKEHISTOGRAM_NEW('Property','Value',...) creates a new SM_MAKEHISTOGRAM_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SM_MakeHistogram_new_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SM_MakeHistogram_new_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_MakeHistogram_new

% Last Modified by GUIDE v2.5 06-Jun-2022 17:30:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SM_MakeHistogram_new_OpeningFcn, ...
                   'gui_OutputFcn',  @SM_MakeHistogram_new_OutputFcn, ...
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


% --- Executes just before SM_MakeHistogram_new is made visible.
function SM_MakeHistogram_new_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SM_MakeHistogram_new (see VARARGIN)

% Choose default command line output for SM_MakeHistogram_new
handles.output = hObject;
handles.histogramData=[];
handles.fitFlag=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_MakeHistogram_new wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SM_MakeHistogram_new_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

handles.workingDir=pwd; % it prints the path of program directory
handles.pathName =[];
handles.fileName=[];
handles.faceColor=[0.9020 0.9020 0.9020];
handles.edgeColor=[0.6510 0.6510 0.6510];
handles.fitColor=[1 0 0];
handles.statement='';
varargout{1} = handles.output;
guidata(hObject, handles);


%% User active buttons...
% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
workingDir=handles.workingDir; % where program is there
[filelist pathName fi] = uigetfile('*_forvBFRET*.dat','Select files on which you wish to perform drift correction.','MultiSelect','on', workingDir);
if iscell(filelist) == 0
    filelist2{1} = filelist;
    filelist = filelist2;
    clear filelist2;
end
handles.pathName=pathName;
handles.filelist=filelist;
guidata(hObject, handles);
handles.histogramData=CalculateHistogram(hObject, eventdata, handles);
handles.fitFlag=0;
guidata(hObject, handles);
plotHistogram(handles.histogramPlot,hObject, eventdata, handles)

% --- Executes on button press in rePlot.
function rePlot_Callback(hObject, eventdata, handles)
plotHistogram(handles.histogramPlot,hObject, eventdata, handles)




function [R2_value,parameter]=fitGaussMath(hObject, eventdata, handles)
centers=handles.histogramData(1,:);
norm_T_histCount=handles.histogramData(2,:);
in_Centers=str2double(split(get(handles.in_Centers,'String'),","));
in_sigma=str2double(split(get(handles.in_sigma,'String'),","));
initialGuesses=[in_Centers in_sigma];
startingGuesses = reshape(initialGuesses', 1, []);
%%%%%%%%%%%%%%%%%%%%%%%%%%
global  c NumTrials TrialError
% 	warning off

% Initializations
NumTrials = 0;  % Track trials
TrialError = 0; % Track errors
% t and y must be row vectors.
tFit = reshape(centers, 1, []);
y = reshape(norm_T_histCount, 1, []);
%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
options.TolFun = 1e-4;
options.TolX = 1e-4;
options.MaxIter = 100000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
[parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate residuals
means = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);
numGaussians = length(c);
thisEstimatedCurve=zeros(length(tFit),numGaussians);
for k = 1 : numGaussians
    % Get each curve.
    thisEstimatedCurve(:,k) = c(k) .* gaussian(tFit, means(k), widths(k));
end
% Overall curve estimate is the sum of the component curves.
yhat = sum(thisEstimatedCurve,2);

% residual sum of squares:SS_res=sum((y-yhat)^2)
SS_res= sum((y-yhat').^2);
SS_tot=sum((y-mean(y)).^2);
R2_value = 1-SS_res/SS_tot;
guidata(hObject, handles);
%% -------------------------------------

function CurveFitData=fitGaussPlots(figureHandle,R2_value,parameter,hObject, eventdata, handles)
global  c
numGaussians = length(c);

centers=handles.histogramData(1,:);
tFit = reshape(centers, 1, []);
means = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);

%% Now plot results.

tFit2 = linspace(tFit(1),tFit(end),1000);
thisEstimatedCurve=zeros(length(tFit2),numGaussians);
CM = jet(numGaussians);
	for k = 1 : numGaussians
		% Get each curve.
		thisEstimatedCurve(:,k) = c(k) .* gaussian(tFit2, means(k), widths(k));
        hold on;
        plot(figureHandle,tFit2, thisEstimatedCurve(:,k),'color',CM(k,:), 'LineWidth', str2double(get(handles.fitThickness,'String'))-1);
	    statement{k}= strcat('peak ', num2str(k),': ',num2str(means(k),1),char(177),num2str(widths(k),1) );
    end
% Overall curve estimate is the sum of the component curves.
yhat2 = sum(thisEstimatedCurve,2);

hold on;
plot(figureHandle, tFit2, yhat2, 'color',handles.fitColor, 'LineWidth', str2double(get(handles.fitThickness,'String')));
hold off;
%%%%%%%%%%%%%%%%%%%%%
uistack(figureHandle, 'top')


CurveFitData=[tFit2' thisEstimatedCurve yhat2];

handles.statement=[strcat('N=',num2str(length(handles.filelist))),statement,{strcat('R2 value: ', num2str(R2_value))}];

textPos=str2double(split(get(handles.textPos,'String'),","));
handles.text=text(figureHandle,textPos(1),textPos(2),handles.statement,'FontSize',str2double(get(handles.axisSize,'String'))-1);
guidata(hObject, handles);



% --- Executes on button press in fitGauss.
function fitGauss_Callback(hObject, eventdata, handles)

plotHistogram(handles.histogramPlot,hObject, eventdata, handles);
[R2_value,parameter]=fitGaussMath(hObject, eventdata, handles);
handles.R2_value=R2_value;
handles.parameter=parameter;
CurveFitData=fitGaussPlots(handles.histogramPlot,R2_value,parameter,hObject, eventdata, handles);
handles.CurveFitData=CurveFitData;
handles.fitFlag=1;
guidata(hObject, handles);

















% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
    dirForSavingHistogram='Histogram';
    % handles.paht comes with '\' at the end
    
    dest_folder=[handles.pathName,dirForSavingHistogram];
    if exist (dest_folder,'dir');
    else mkdir(handles.pathName,dirForSavingHistogram);
    end
    %-----------------------------------------------------
    NamePrefix=get(handles.plotTitle, 'String');
    NamePrefix = regexprep(NamePrefix,'[ ]','');
    NamePrefix = regexprep(NamePrefix,'[,.;:/\~!@#$%^&*()_+]','_');
    
    fName_Picture=[dest_folder '\' NamePrefix];

%% redraw the entire figure again! Stupid!! but that is how you do it.
F=figure();

plotHandle = axes('Parent',F);
plotHistogram(plotHandle,hObject, eventdata, handles);
 
 if handles.fitFlag==1
   CurveFitData=fitGaussPlots(plotHandle,handles.R2_value,handles.parameter,hObject, eventdata, handles);
   save(strcat(fName_Picture,'_CurveFitData.dat'), 'CurveFitData', '-ascii');
 end
 
 box off;
 %% Now plot results.
saveas(gcf,strcat(fName_Picture,'.jpg'));
saveas(gcf,strcat(fName_Picture,'.eps'));


HistData=handles.histogramData';



save(strcat(fName_Picture,'_Histdata.dat'), 'HistData', '-ascii');

 close(F);




















%% ========================  helper functions   ====================================

function histogramData=CalculateHistogram(hObject, eventdata, handles)
BinSize =str2double(get(handles.binSize,'String'));
edges = -0.2+BinSize/2:BinSize:1.2+BinSize/2;
centers=edges(1)+BinSize/2:BinSize:edges(end)-BinSize/2;
pathName=handles.pathName;
filelist=handles.filelist;
T_histCount=zeros(size(centers));
for i=1:length(filelist)
data=importdata(strcat(pathName,filelist{i}));
%FRET=data(:,3)./(data(:,2)+data(:,3));
FRET=data(:,2)./(data(:,1)+data(:,2));
histCount = histcounts(FRET,edges);
histCount=histCount./size(data,1);
T_histCount=T_histCount+histCount;
end
norm_T_histCount=T_histCount./length(filelist);
histogramData=[centers; norm_T_histCount];
%=======================================================================================================================================================
function plotHistogram(plotHandle,hObject, eventdata, handles)
centers=handles.histogramData(1,:);
norm_T_histCount=handles.histogramData(2,:);
axes(plotHandle)
bar(centers, norm_T_histCount, 'FaceColor',handles.faceColor, 'EdgeColor',handles.edgeColor, 'LineWidth',str2double(get(handles.edgeWidth,'String')));
set(plotHandle, 'XLim', [-0.2 1.2], 'XTick', -0.2:0.2:1.2,'XTickLabel', -0.2:0.2:1.2);
plotHandle.FontSize =str2double(get(handles.tickSize,'String'));
ylabel(plotHandle,'Normalized Count', 'FontSize',str2double(get(handles.axisSize,'String')) );
xlabel(plotHandle,'FRET Efficiency', 'FontSize',str2double(get(handles.axisSize,'String')) );
title(plotHandle,get(handles.plotTitle,'String'), 'FontSize', str2double(get(handles.titleSize,'String')));
guidata(hObject, handles);
%=======================================================================================================================================================
function theError = fitgauss(lambda, t, y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006

global c NumTrials TrialError
try
	
	A = zeros(length(t), round(length(lambda) / 2));
	for j = 1 : length(lambda) / 2
		A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';
	end
	
	c = A \ y';
	z = A * c;
	theError = norm(z - y');
	
	% Penalty so that heights don't become negative.
	if sum(c < 0) > 0
		theError = theError + 1000000;
	end
	
	NumTrials = NumTrials + 1;
	TrialError(NumTrials) = theError;
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
%=======================================================================================================================================================
function callStackString = GetCallStack(errorObject)
try
	theStack = errorObject.stack;
	callStackString = '';
	stackLength = length(theStack);
	% Get the date of the main, top level function:
	% 	d = dir(theStack(1).file);
	% 	fileDateTime = d.date(1:end-3);
	if stackLength <= 3
		% Some problem in the OpeningFcn
		% Only the first item is useful, so just alert on that.
		[folder, baseFileName, ext] = fileparts(theStack(1).file);
		baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
		callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(1).name, theStack(1).line);
	else
		% Got past the OpeningFcn and had a problem in some other function.
		for k = 1 : length(theStack)-3
			[folder, baseFileName, ext] = fileparts(theStack(k).file);
			baseFileName = sprintf('%s%s', baseFileName, ext);	% Tack on extension.
			callStackString = sprintf('%s in file %s, in the function %s, at line %d\n', callStackString, baseFileName, theStack(k).name, theStack(k).line);
		end
	end
catch ME
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\nError Message:\n%s', ...
		mfilename, ME.message);
	WarnUser(errorMessage);
end
%=======================================================================================================================================================
% Pops up a warning message, and prints the error to the command window.
function WarnUser(warningMessage)
if nargin == 0
	return; % Bail out if they called it without any arguments.
end
try
	fprintf('%s\n', warningMessage);
	uiwait(warndlg(warningMessage));
	% Write the warning message to the log file
	folder = 'C:\Users\Public\Documents\MATLAB Settings';
	if ~exist(folder, 'dir')
		mkdir(folder);
	end
	fullFileName = fullfile(folder, 'Error Log.txt');
	fid = fopen(fullFileName, 'at');
	if fid >= 0
		fprintf(fid, '\nThe error below occurred on %s.\n%s\n', datestr(now), warningMessage);
		fprintf(fid, '-------------------------------------------------------------------------------\n');
		fclose(fid);
	end
catch ME
	message = sprintf('Error in WarnUser():\n%s', ME.message);
	fprintf('%s\n', message);
	uiwait(warndlg(message));
end
%=======================================================================================================================================================
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
%=======================================================================================================================================================
function yhat = PlotComponentCurves(x, y, t, c, parameter)
try
	fontSize = 20;
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
	hFig2 = figure;
	hFig2.Name = 'Fitted Component Curves';
	% 	plot(x, y, '--', 'LineWidth', 2)
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
	legendStrings = cell(numGaussians + 2, 1);
	for k = 1 : numGaussians
		% Get each component curve.
		thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
		% Plot component curves.
		plot(x, thisEstimatedCurve, '-', 'LineWidth', 2);
		hold on;
		% Overall curve estimate is the sum of the component curves.
		yhat = yhat + thisEstimatedCurve;
		legendStrings{k} = sprintf('Estimated Gaussian %d', k);
	end
	% Plot original summation curve, that is the actual curve.
	plot(x, y, 'r-', 'LineWidth', 1)
	% Plot estimated summation curve, that is the estimate of the curve.
	plot(x, yhat, 'k--', 'LineWidth', 2)
	grid on;
	xlabel('X', 'FontSize', fontSize)
	ylabel('Y', 'FontSize', fontSize)
	caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
	title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
	grid on
	legendStrings{numGaussians+1} = sprintf('Actual original signal');
	legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
	legend(legendStrings);
	xlim(sort([x(1) x(end)]));
	hFig2.WindowState = 'maximized';
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
% --- Executes on button press in fitColor.
function fitColor_Callback(hObject, eventdata, handles)
handles.fitColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);












%% Functions that are not actively needed....
% --- Executes on button press in faceColor.
function faceColor_Callback(hObject, eventdata, handles)
handles.faceColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);
% --- Executes on button press in edgeColor.
function edgeColor_Callback(hObject, eventdata, handles)
handles.edgeColor = uisetcolor([0.6 0.8 1]);
guidata(hObject, handles);


function titleSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function titleSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to titleSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function axisSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function axisSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tickSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function tickSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fitThickness_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function fitThickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function textSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function textPos_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function textPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function savingfileName_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function savingfileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savingfileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in_Centers_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function in_Centers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in_Centers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function in_sigma_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function in_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to in_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edgeWidth_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edgeWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgeWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTitle_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function plotTitle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function binSize_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function binSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
