function varargout = stormtracker(varargin)
% stormtracker MATLAB code for stormtracker.fig
%      stormtracker, by itself, creates a new stormtracker or raises the existing
%      singleton*.
%
%      H = stormtracker returns the handle to a new stormtracker or the handle to
%      the existing singleton*.
%
%      stormtracker('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in stormtracker.M with the given input arguments.
%
%      stormtracker('Property','Value',...) creates a new stormtracker or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stormtracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stormtracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stormtracker

% Last Modified by GUIDE v2.5 11-Apr-2013 11:34:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @stormtracker_OpeningFcn, ...
    'gui_OutputFcn',  @stormtracker_OutputFcn, ...
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


% --- Executes just before stormtracker is made visible.
function stormtracker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stormtracker (see VARARGIN)

% Choose default command line output for stormtracker
handles.output = hObject;
% use guiMain to initialise appData
guiMain('init',handles, varargin);
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = stormtracker_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%% create functions

function trackingWindow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTracksMinSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotTracksColour_CreateFcn(hObject, ~, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sigmaNoise_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bacLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bacWidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rangeD_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DhistMinSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function initGuess_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CDFfittingOptions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function localizationThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DhistThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function curveColor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function testDetectionFrame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderTestDetectionFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





%%% callbacks

function tracking_Callback(hObject, eventdata, handles)
guiMain('tracking_Callback',handles);

function trackingWindow_Callback(hObject, eventdata, handles)
guiMain('trackingWindow_Callback',handles);

function plotTracks_Callback(hObject, eventdata, handles)
guiMain('plotTracks_Callback',handles);

function plotTracksMinSteps_Callback(hObject, eventdata, handles)
guiMain('plotTracksMinSteps_Callback',handles);

function plotTracksColour_Callback(hObject, eventdata, handles)
guiMain('plotTracksColour_Callback',handles);

function checkboxSaveTracks_Callback(hObject, eventdata, handles)
guiMain('checkboxSaveTracks_Callback',handles);

function MSDcurve_Callback(hObject, eventdata, handles)
guiMain('MSDcurve_Callback',handles);

function pixel_Callback(hObject, eventdata, handles)
guiMain('pixel_Callback',handles);

function dT_Callback(hObject, eventdata, handles)
guiMain('dT_Callback',handles);

function sigmaNoise_Callback(hObject, eventdata, handles)
guiMain('sigmaNoise_Callback',handles);

function bacLength_Callback(hObject, eventdata, handles)
guiMain('bacLength_Callback',handles);

function bacWidth_Callback(hObject, eventdata, handles)
guiMain('bacWidth_Callback',handles);

function maxStep_Callback(hObject, eventdata, handles)
guiMain('maxStep_Callback',handles);

function clearData_Callback(hObject, eventdata, handles)
guiMain('clearData_Callback',handles);

function histD_Callback(hObject, eventdata, handles)
guiMain('histD_Callback',handles);

function rangeD_Callback(hObject, eventdata, handles)
guiMain('rangeD_Callback',handles);

function DhistMinSteps_Callback(hObject, eventdata, handles)
guiMain('DhistMinSteps_Callback',handles);

function CDFcurve_Callback(hObject, eventdata, handles)
guiMain('CDFcurve_Callback',handles);

function initGuess_Callback(hObject, eventdata, handles)
guiMain('initGuess_Callback',handles);

function CDFfittingOptions_Callback(hObject, eventdata, handles)
guiMain('CDFfittingOptions_Callback',handles);

function plotLocalizations_Callback(hObject, eventdata, handles)
guiMain('plotLocalizations_Callback',handles);

function selectROIs_Callback(hObject, eventdata, handles)
guiMain('selectROIs_Callback',handles);

function localization_Callback(hObject, eventdata, handles)
guiMain('localization_Callback',handles);

function localizationThresh_Callback(hObject, eventdata, handles)
guiMain('localizationThresh_Callback',handles);

function findBindingDiffusionTracks_Callback(hObject, eventdata, handles)
guiMain('findBindingDiffusionTracks_Callback',handles);

function DThresh_Callback(hObject, eventdata, handles)
guiMain('DThresh_Callback',handles);

function TwoSpeciesMSDThreshold_Callback(hObject, eventdata, handles)
guiMain('TwoSpeciesMSDThreshold_Callback',handles);

function DhistThresh_Callback(hObject, eventdata, handles)
guiMain('DhistThresh_Callback',handles);

function holdFigureCheckbox_Callback(hObject, eventdata, handles)
guiMain('holdFigureCheckbox_Callback',handles);

function curveColor_Callback(hObject, eventdata, handles)
guiMain('curveColor_Callback',handles);

function combineTracks_Callback(hObject, eventdata, handles)
guiMain('combineTracks_Callback',handles);

function plotBrightfield_Callback(hObject, eventdata, handles)
guiMain('plotBrightfield_Callback',handles);

function countMoleculesInROIs_Callback(hObject, eventdata, handles)
guiMain('countMoleculesInROIs_Callback',handles);

function localizationProfile_Callback(hObject, eventdata, handles)
guiMain('localizationProfile_Callback',handles);

function DThresholdInROIs_Callback(hObject, eventdata, handles)
guiMain('DThresholdInROIs_Callback',handles);

function selectCellROIs_Callback(hObject, eventdata, handles)
guiMain('selectCellROIs_Callback',handles);

function testDetection_Callback(hObject, eventdata, handles)
guiMain('testDetection_Callback',handles);

function testDetectionFrame_Callback(hObject, eventdata, handles)
guiMain('testDetectionFrame_Callback',handles);

function sliderTestDetectionFrame_Callback(hObject, eventdata, handles)
guiMain('sliderTestDetectionFrame_Callback',handles);

function transformGreenBlue_Callback(hObject, eventdata, handles)
guiMain('transformGreenBlue_Callback',handles);

function coLocalisation_Callback(hObject, eventdata, handles)
guiMain('coLocalisation_Callback',handles);

function closeAll_Callback(hObject, eventdata, handles)
guiMain('closeAll_Callback',handles);

function editLocalisations_Callback(hObject, eventdata, handles)
guiMain('editLocalisations_Callback',handles);

function histogram2dOverlay_Callback(hObject, eventdata, handles)
guiMain('histogram2dOverlay_Callback',handles);

function temporaryFunction_Callback(hObject, eventdata, handles)
guiMain('temporaryFunction_Callback',handles);


%------------------------------------------------------------------------
%--------application specific code
function guiMain(param, handles,varargin)
% function guiMain(param, handles)
% ----------Main control function ---------------------------------------

if strcmp(param, 'init') % on initialise, create the appData variable & initialise fields
    
    appData.trackParams.mem = 0;
    appData.trackParams.dim = 2;
    appData.trackParams.good = 0;
    appData.trackParams.quiet = 0;
    appData.trackParams.maxDisp = 5;
    
    set(handles.trackingWindow,'String',num2str(appData.trackParams.maxDisp));
    
    appData.plotTracksMinSteps = 1;
    set(handles.plotTracksMinSteps,'String',num2str(appData.plotTracksMinSteps));
    
    appData.plotTracksColour = 'random';
    set(handles.plotTracksColour,'String',appData.plotTracksColour);
    
    setappdata(handles.figure1, 'appData', appData);
    
    appData.checkboxSaveTracks = 0;
    set(handles.checkboxSaveTracks,'Value',appData.checkboxSaveTracks);
    
    appData.pixel = 0.1145;
    set(handles.pixel,'String',num2str(appData.pixel));
    
    appData.dT = 0.01526;
    set(handles.dT,'String',num2str(appData.dT));
    
    appData.sigmaNoise = 0.04/appData.pixel;
    set(handles.sigmaNoise,'String',num2str(appData.sigmaNoise*appData.pixel));
    
    appData.bacLength = 0.6; % length of rectangular section, not including round caps
    set(handles.bacLength,'String',num2str(appData.bacLength));
    
    appData.bacWidth = 0.6;
    set(handles.bacWidth,'String',num2str(appData.bacWidth));
    
    appData.bacArea = appData.bacLength * appData.bacWidth / (appData.pixel^2);
    
    appData.maxStep = 8;
    set(handles.maxStep,'String',num2str(appData.maxStep));
    
    appData.rangeDString = '-0.225:0.05:2'; % D range histD
    appData.rangeD = str2num(appData.rangeDString);
    set(handles.rangeD,'String',appData.rangeDString);
    
    appData.DhistMinSteps = 4;
    set(handles.DhistMinSteps,'String',num2str(appData.DhistMinSteps));
    
    appData.initGuess = 1;
    set(handles.initGuess,'String',num2str(appData.initGuess));
    
    appData.CDFfittingOption = 2;
    set(handles.CDFfittingOptions,'Value',appData.CDFfittingOption);
    
    appData.localizationWindow = 7;
    appData.localizationThresh = 12;
    set(handles.localizationThresh,'String',num2str(appData.localizationThresh));
    
    appData.nFiles = 1;
    
    appData.DThreshString = '0.15 0.075'; %findBindingDiffusionTracks Threshold
    appData.DThresh = str2num(appData.DThreshString);
    set(handles.DThresh,'String',appData.DThreshString);
    
    appData.DhistThresh = 0.1; %TwoSpeciesMSDThreshold Threshold
    set(handles.DhistThresh,'String',num2str(appData.DhistThresh));
    
    appData.holdFigureCheckbox = 0;
    set(handles.holdFigureCheckbox,'Value',appData.holdFigureCheckbox);
    
    appData.curveColor = 'b';
    set(handles.curveColor,'String',appData.curveColor);
    
    appData.testDetectionFrame = 1;
    set(handles.testDetectionFrame,'String',num2str(appData.testDetectionFrame));
    
    appData.sliderTestDetectionFrame = 1;
    set(handles.sliderTestDetectionFrame,'value',appData.sliderTestDetectionFrame);
    
    appData.sliderTestDetectionFrameMax = 10000;
    set(handles.sliderTestDetectionFrame,'max',appData.sliderTestDetectionFrameMax);
    
    appData.sliderTestDetectionFrameMin = 1;
    set(handles.sliderTestDetectionFrame,'min',appData.sliderTestDetectionFrameMin);
    
    appData.sliderTestDetectionFrameStep = [.0001 .01];
    set(handles.sliderTestDetectionFrame,'SliderStep',appData.sliderTestDetectionFrameStep);
    
else  % subsequent runs, retrieve appData
    if isappdata(handles.figure1, 'appData')
        appData = getappdata(handles.figure1, 'appData');
    else
        error('appData, main data variable structure for avgGui not initialised');
    end
end


switch param
    
    case 'localization_Callback'
        
        [movieFilename, moviePathname] =...
            uigetfile('*.fits', 'MatFile data:','MultiSelect', 'on');
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
            
            info = whos('movieFilename');
            if strcmp(info.class,'char')
                appData.nFiles = 1;
            else
                appData.nFiles = numel(movieFilename);
            end
            
            for ii = 1:appData.nFiles
                
                if appData.nFiles == 1
                    loadname  = [moviePathname movieFilename];
                else
                    loadname = [moviePathname movieFilename{1,ii}];
                end
                
                gaussStorm(loadname,...
                    appData.localizationThresh, appData.localizationWindow);
                
            end
            
            if appData.nFiles == 1
                appData.data = importdata([moviePathname movieFilename(1:end-5),...
                    '_thresh',num2str(appData.localizationThresh),'.gaussstorm.pos.out']);
            end
            
        end
        
        
    case 'localizationThresh_Callback'
        
        localizationThresh = str2double(get(handles.localizationThresh,'String'));
        
        if isnan(localizationThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.localizationThresh = localizationThresh;
        
        
    case 'tracking_Callback'
        
        [dataFilename, dataPathname] =...
            uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
        
        if ~(isnumeric(dataFilename)&&dataFilename==0) %check the user has not pressed cancel
            
            info = whos('dataFilename');
            if strcmp(info.class,'char')
                appData.nFiles = 1;
            else
                appData.nFiles = numel(dataFilename);
            end
            
            for ii = 1:appData.nFiles
                
                if appData.nFiles == 1
                    loadname  = [dataPathname dataFilename];
                else
                    loadname = [dataPathname dataFilename{1,ii}];
                end
                
                newData = importdata(loadname);
                if isstruct(newData)
                    appData.data = newData.data;
                else
                    appData.data = newData;
                end
                clear newData;
                if ~isempty(appData.data)
                    pos = zeros(length(appData.data(:,1)),3);

                    % standard indexing
                    pos(:,1) = appData.data(:,2);
                    pos(:,2) = appData.data(:,3);
                    pos(:,3) = appData.data(:,1);

                    % Select for localised

                    appData.trackParams.mem = 1;

                    tracks = trackWithDummy(pos, appData.trackParams);
                    nMolecules = max(tracks(:,4));
                    appData.tracks = tracks;
                    appData.nMolecules = nMolecules;

                    set(handles.results,'String',['nMolecules = ' num2str(nMolecules)]);


                    if appData.checkboxSaveTracks == 1
                        save([loadname '.disp' num2str(appData.trackParams.maxDisp) ...
                            '.mem' num2str(appData.trackParams.mem) ...
                            '.tracks'], 'tracks');
                    end
                end
                
            end
            
        end
        
        
    case 'trackingWindow_Callback'
        
        trackingWindow = str2double(get(handles.trackingWindow,'String'));
        
        if isnan(trackingWindow)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.trackParams.mem = 0;
        appData.trackParams.dim = 2;
        appData.trackParams.good = 0;
        appData.trackParams.quiet = 0;
        appData.trackParams.maxDisp = trackingWindow;
        
        
    case 'plotTracks_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
                appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        mymap = colormap(jet(5000));
        
            trackStart = [];
            
          figure(1)
          hold all
        molecInds = unique(appData.tracks(:,4));
        colormap gray
         
              
        if  isfield(appData,'tracks')
            pos = zeros(max(appData.tracks(:,4)),2);
           for tNum = 1:length(molecInds) % for each track
               if mod(tNum,100)==0
               A = ceil(tNum*100/length(molecInds))
               end
                xx = find(appData.tracks(:,4)==molecInds(tNum));
                 if length(xx)>4
%                 plot(appData.tracks(xx,1),appData.tracks(xx,2),'-')
%                 end
%                 if ((mean(appData.tracks(xx,1))>90) && (mean(appData.tracks(xx,2)<60)))
                plot(mean(appData.tracks(xx,1)),mean(appData.tracks(xx,2)),'.',...
           'MarkerFaceColor',mymap(ceil(appData.tracks(xx(1),3)),:),...
           'MarkerEdgeColor',mymap(ceil(appData.tracks(xx(1),3)),:),'MarkerSize',10)
            text(mean(appData.tracks(xx,1)),mean(appData.tracks(xx,2)),...
                [num2str(appData.tracks(xx(1),3)),'-',num2str(appData.tracks(xx(end),3))],'FontSize',1)
                end
                if ~isempty(xx) 
                    if length(xx)>=4
                    trackStart = [trackStart; appData.tracks(xx(1),3)];
                    end
                end
            end
            figure
            [n, bin] = histc(trackStart,1:50:5000);
            xBar = (1:50:5000)*.015;
            bar(xBar,n)
            xlabel('Time [seconds]')
            ylabel('Number of Tracks (molecules)')
            
            hold off
            
           
%             set(handles.results,'String',['nMolecules = ' num2str(max(appData.tracks(:,4)))]);
%             plotTracksColourCoded(appData.tracks,appData)
        end
        
        
    case 'plotTracksMinSteps_Callback'
        
        plotTracksMinSteps = str2double(get(handles.plotTracksMinSteps,'String'));
        
        if isnan(plotTracksMinSteps)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.plotTracksMinSteps = plotTracksMinSteps;
        
        
    case 'plotTracksColour_Callback'
        
        plotTracksColour = get(handles.plotTracksColour,'String');
        
        appData.plotTracksColour = plotTracksColour;
        
        
    case 'checkboxSaveTracks_Callback'
        
        appData.checkboxSaveTracks = get(handles.checkboxSaveTracks,'Value');
        
        
    case 'MSDcurve_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
                appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        
        if  isfield(appData,'tracks')
            
            [appData.D appData.locD appData.locFrameD appData.locFrameConfineD] =...
                averageFrameMSD(appData.tracks, appData);
            
            assignin('base', 'MSDresults', [appData.D appData.locD appData.locFrameD appData.locFrameConfineD]);
            
            disp_str = {['D = ' num2str(appData.D)];...
                ['locD = ' num2str(appData.locD)];...
                ['locFrameD = ' num2str(appData.locFrameD)];...
                ['locFrameConfineD = ' num2str(appData.locFrameConfineD)]};
            
            set(handles.results,'String',disp_str);
            
        end
        
        
    case 'pixel_Callback'
        
        pixel = str2double(get(handles.pixel,'String'));
        
        if isnan(pixel)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.pixel = pixel;
        
        
    case 'dT_Callback'
        
        dT = str2double(get(handles.dT,'String'));
        
        if isnan(dT)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.dT = dT;
        
        
    case 'sigmaNoise_Callback'
        
        sigmaNoise = str2double(get(handles.sigmaNoise,'String'))/appData.pixel;
        
        if isnan(sigmaNoise)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.sigmaNoise = sigmaNoise;
        
        
    case 'bacLength_Callback'
        
        bacLength = str2double(get(handles.bacLength,'String'));
        
        if isnan(bacLength)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.bacLength = bacLength;
        appData.bacArea = appData.bacWidth*appData.bacLength / (appData.pixel^2);
        
        
    case 'bacWidth_Callback'
        
        bacWidth = str2double(get(handles.bacWidth,'String'));
        
        if isnan(bacWidth)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.bacWidth = bacWidth;
        %also update bacArea
        appData.bacArea = appData.bacWidth*appData.bacLength / (appData.pixel^2);
        
        
    case 'maxStep_Callback'
        
        maxStep = str2double(get(handles.maxStep,'String'));
        
        if isnan(maxStep)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.maxStep = maxStep;
        
        
    case 'histD_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
                appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')
            
            appData.histD = histD2(appData.tracks,appData);
            hold on;
            
        end
        
        
    case 'rangeD_Callback'
        
        rangeD = str2num(get(handles.rangeD,'String'));
        
        if isnan(rangeD)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.rangeD = rangeD;
        
        
    case 'DhistMinSteps_Callback'
        
        DhistMinSteps = str2double(get(handles.DhistMinSteps,'String'));
        
        if isnan(DhistMinSteps)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DhistMinSteps = DhistMinSteps;
        
        
    case 'CDFcurve_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
                appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')
            [appData.CDFresults] =...
                CDFAnalysis(appData.tracks, appData, appData.initGuess, appData.CDFfittingOption);
            
            assignin('base', 'CDFresults', appData.CDFresults);
            
            disp_str = {['CDFresults = ' num2str(appData.CDFresults)]};
            
            set(handles.results,'String',disp_str);
        end
        
        
    case 'CDFfittingOptions_Callback'
        
        appData.CDFfittingOption = get(handles.CDFfittingOptions,'Value');
        
        if appData.CDFfittingOption == 1
            appData.initGuess = [];
            set(handles.initGuess,'String',num2str(appData.initGuess));
            disp_str = 'CDF curve: no fit - no initial guess';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 2
            appData.initGuess = 1;
            set(handles.initGuess,'String',num2str(appData.initGuess));
            disp_str = 'CDF curve: initial guess for D [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 3
            appData.initGuess = [0.5 0.1 1];
            set(handles.initGuess,'String',...
                [num2str(appData.initGuess(1)) ' ' num2str(appData.initGuess(2)) ' ' num2str(appData.initGuess(3))]);
            disp_str = 'CDF curve: initial guess for a, fixed values for D1, D2 [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        elseif appData.CDFfittingOption == 4
            appData.initGuess = [0.5 0.1 1];
            set(handles.initGuess,'String',...
                [num2str(appData.initGuess(1)) ' ' num2str(appData.initGuess(2)) ' ' num2str(appData.initGuess(3))]);
            disp_str = 'CDF curve: initial guess for a, D1, D2 [um^2/s]';
            set(handles.initGuessText,'String',disp_str);
            
        end
        
    case 'initGuess_Callback'
        
        initGuess = str2num(get(handles.initGuess,'String'));
        
        if isnan(initGuess)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.initGuess = initGuess;
        
        
    case 'plotLocalizations_Callback'
        
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
                newData = importdata([appData.dataPathname appData.dataFilename]);
                if isstruct(newData)
                    appData.data = newData.data;
                else
                    appData.data = newData;
                end
                clear newData;
            end
        end
        
        if isfield(appData,'data')
            
            % standard indexing
            pos(:,1) = appData.data(:,2);
            pos(:,2) = appData.data(:,3);
            
            % Clustering algorithm
            %clusterThreshold = 0.04/appData.pixel;
            %[nPositionsInCluster, distancesWithinClusters, distancesBetweenClusters] = nearestNeighbourClustering(pos, clusterThreshold, appData);
            
            if appData.holdFigureCheckbox
                figure(appData.figureHandle);
            else
                figure;
            end
            hold all
              
            hs = plot(pos(:,1),pos(:,2),'b.','MarkerSize',1);
            
            %hs = scatter(pos(:,1),pos(:,2),'b.')
            
            xlabel('x [pixels]');
            ylabel('y [pixels]');
            axis image;
            
        end
        
        
    case 'selectROIs_Callback'
        
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
                newData = importdata([appData.dataPathname appData.dataFilename]);
                if isstruct(newData)
                    appData.data = newData.data;
                else
                    appData.data = newData;
                end
                clear newData;
            end
        end
        
        if isfield(appData,'data')
            
            ROIData = selectCellROIs(appData.data,appData);
            
            appData.data = ROIData;
            save([appData.dataPathname appData.dataFilename 'mem1.ROIData'],'ROIData')
        end
        
        
    case 'findBindingDiffusionTracks_Callback'
        
        if  ~isfield(appData,'tracks')
            [appData.tracksFilename, appData.tracksPathname] = ...
                uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
            
            if ~(isnumeric(appData.tracksFilename)&&appData.tracksFilename==0) %check the user has not pressed cancel
                appData.tracks = importdata([appData.tracksPathname appData.tracksFilename]);
            end
            
        end
        
        if  isfield(appData,'tracks')
            findBindingDiffusionTracks(appData.tracks, appData)
        end
        
        
    case 'DThresh_Callback'
        
        DThresh = str2num(get(handles.DThresh,'String'));
        
        if isnan(DThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DThresh = DThresh;
        
        
        %     case 'TwoSpeciesMSDThreshold_Callback'
        %
        %         if  ~isfield(appData,'dataStruct')
        %             [appData.dataStructFilename, appData.dataStructPathname] = ...
        %                 uigetfile('*.dataStruct', 'MatFile data:','MultiSelect', 'on');
        %
        %             if ~(isnumeric(appData.dataStructFilename)&&appData.dataStructFilename==0) %check the user has not pressed cancel
        %                 appData.dataStruct = importdata([appData.dataStructPathname appData.dataStructFilename]);
        %             end
        %
        %         end
        %
        %         if  isfield(appData,'dataStruct')
        %
        %             [appData.diffusionFraction, appData.D1, appData.D2] =...
        %                 twoSpeciesMSDThreshold(appData.dataStruct.ROIData(1,1).tracks, appData);
        %
        %             assignin('base', 'MSDThreshResults', [appData.diffusionFraction, appData.D1, appData.D2]);
        %
        %             disp_str = {['D histogram threshold results = '...
        %                 num2str(appData.diffusionFraction) ', ' num2str(appData.D1) ', '  num2str(appData.D2)]};
        %
        %             set(handles.results,'String',disp_str);
        %
        %         end
        
        
    case 'TwoSpeciesMSDThreshold_Callback'
        
        [tracksFilename, tracksPathname] = ...
            uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
        
        if ~(isnumeric(tracksFilename)&&tracksFilename==0) %check the user has not pressed cancel
            
            info = whos('tracksFilename');
            if strcmp(info.class,'char')
                appData.nFiles = 1;
            else
                appData.nFiles = numel(tracksFilename);
            end
            
            diffusionFraction = zeros(appData.nFiles,1);
            %D1 = zeros(appData.nFiles,1);
            %D2 = zeros(appData.nFiles,1);
            
           D1All= [];
           D2All = [];
           groupD1 = [];
           groupD2 = [];
           meanD1 = [];
           meanD2 = [];
           binCenterOfPeakAll = [];
           meanSqDiffAll = [];
           
    
            for ii = 1:appData.nFiles
                
                if appData.nFiles == 1
                    loadname  = [tracksPathname tracksFilename];
                else
                    loadname = [tracksPathname tracksFilename{1,ii}];
                end
                
                newData = importdata(loadname);
                if isstruct(newData)
                    appData.tracks = newData.tracks;
                else
                    appData.tracks = newData;
                end
                clear newData;
                
                [diffusionFraction(ii), D1, D2, binCenterOfPeak, meanSqDiff] =...
                    twoSpeciesMSD2Threshold_vStephan(appData.tracks, appData);
                
                disp_str = {['D histogram threshold results = '...
                    num2str(diffusionFraction(ii)) ', ' num2str(D1(ii)) ', '  num2str(D2(ii))]};
                
                binCenterOfPeakAll = [binCenterOfPeakAll binCenterOfPeak];
                
                D1All = [D1All; D1];
                groupD1 = [groupD1; ii*ones(length(D1),1)];
                meanD1 = [meanD1; mean(D1)];
                
                D2All = [D2All; D2]; 
                groupD2 = [groupD2; ii*ones(length(D2),1)];
                meanD2 = [meanD2; mean(D2)]; 
                
                meanSqDiffAll = [meanSqDiffAll meanSqDiff];
                
                set(handles.results,'String',disp_str);
                
            end
            
            %assignin('base', 'MSDThreshResults', [diffusionFraction, D1, D2]);
            if length(D1All)>1
                
                figure(100),
                boxplot(D2All,groupD2)
            
                figure(101)
                plot(meanD1, '*')
             
                figure(102)
                plot(binCenterOfPeakAll,'*')
                
                figure(103)
                plot(meanSqDiffAll,'*')
            
            end
            
        end
        
        
    case 'DhistThresh_Callback'
        
        DhistThresh = str2num(get(handles.DhistThresh,'String'));
        
        if isnan(DhistThresh)
            set(hObject, 'String', 0);
            errordlg('Input must be a number','Error');
        end
        
        appData.DhistThresh = DhistThresh;
        
        
    case 'holdFigureCheckbox_Callback'
        
        appData.holdFigureCheckbox = get(handles.holdFigureCheckbox,'Value');
        if appData.holdFigureCheckbox
            appData.figureHandle = figure;
        end
        
        
    case 'curveColor_Callback'
        
        appData.curveColor = get(handles.curveColor,'String');
        
        
    case 'combineTracks_Callback'
        
        appData.tracks = combineTracks;
        
        
    case 'plotBrightfield_Callback'
        
        reply = input('Do you want to plot an averaged fits movie (press enter) or plot a tif image (press t)? a/b: ', 's');
        if isempty(reply)
            
            [brightfieldFilename, brightfieldPathname] = ...
                uigetfile('*.fits', 'fits data:');
            
            if ~(isnumeric(brightfieldFilename)&&brightfieldFilename==0) %check the user has not pressed cancel
                ImageInfo = fits_read_header([brightfieldPathname brightfieldFilename]);
                imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
                
                if any(strcmp(fieldnames(ImageInfo),'NAXIS3'))
                    nFrames = ImageInfo.NAXIS3;
                else
                    nFrames = 1;
                end
                dataIn=ImageStack([brightfieldPathname brightfieldFilename], imageLim);
                averageImage = double(getFrame(dataIn,1));
                
                if nFrames >= 2
                    
                    for ii = 2:nFrames
                        image = double(getFrame(dataIn,ii));
                        averageImage = averageImage + image;
                    end
                    
                    averageImage = averageImage / double(nFrames);
                    
                end
                
                averageImage = flipud(averageImage);
                
                if appData.holdFigureCheckbox;
                    figure(appData.figureHandle);
                    hold all;
                else
                    figure;
                end
                
                imshow(averageImage,[min(min(averageImage)) max(max(averageImage))]);
            end
            
        elseif reply == 't'
            
            [brightfieldFilename, brightfieldPathname] = ...
                uigetfile('*.tif', 'tif data:');
            
            if ~(isnumeric(brightfieldFilename)&&brightfieldFilename==0)
                image = importdata([brightfieldPathname brightfieldFilename]);
                image = flipud(image);
                if appData.holdFigureCheckbox;
                    figure(appData.figureHandle);
                    hold all;
                else
                    figure;
                end
                
                imshow(image,[min(min(image)) max(max(image))])
                
            end
            
        end
        
        
        
    case 'countMoleculesInROIs_Callback'
        
        countMoleculesInROIs(appData);
        
        
    case 'localizationProfile_Callback'
        
        localizationProjection_v2(appData)
        
        
    case 'DThresholdInROIs_Callback'
        
        DThresholdInROIs(appData);
        
        
    case 'selectCellROIs_Callback'
        reply = input('Do you want to use Schnitzcell (press enter) or load an ROIData File (press r/R)? ', 's');
        if isempty(reply)
            UseSchnitzcellForSegmentation
        elseif (reply == 'r' || reply=='R')
            if -1
                % NOT RUN
                if ~isfield(appData,'data')
                    [appData.dataFilename, appData.dataPathname] = uigetfile('*.out', 'MatFile data:','MultiSelect', 'on');
                    
                    if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
                        newData = importdata([appData.dataPathname appData.dataFilename]);
                        if isstruct(newData)
                            appData.data = newData.data;
                        else
                            appData.data = newData;
                        end
                        clear newData;
                    end
                end
                
                if isfield(appData,'data')
                    
                    ROIData = selectCellROIs(appData.data,appData);
                    
                    save([appData.dataFilename '.ROIData_SU'], 'ROIData');
                    
                end
            else
                [ROIData, appData.dataFilename] = selectCellROIsImpoly(appData);
                save([appData.dataFilename '.ROIData_SU'], 'ROIData');
            end
        end
        
        
    case 'testDetection_Callback'
        
        %From the Brightfield
        if ~isfield(appData,'testFile')
            [movieFilename, moviePathname] = ...
                uigetfile('*.*', 'MatFile data:','MultiSelect', 'off');
            
            appData.testFile.movieFilename = movieFilename;
            appData.testFile.moviePathname = moviePathname;
        else
            movieFilename = appData.testFile.movieFilename;
            moviePathname = appData.testFile.moviePathname;
        end
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
            ImageInfo = fits_read_header([moviePathname movieFilename]);
            imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
            
            nFrames = ImageInfo.NAXIS3;
            dataIn=ImageStack([moviePathname movieFilename], imageLim);
            
            singleImage = flipud(double(getFrame(dataIn,appData.testDetectionFrame)));
            
            if appData.holdFigureCheckbox;
                figure(appData.figureHandle);
                hold all;
            else
                figure;
            end
            
            imshow(singleImage,[min(min(singleImage)) max(max(singleImage))]);
            
            %From the Localisation
            if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
                
                
                loadname  = [moviePathname movieFilename];
                
                [pos(:,1),pos(:,2)] = gaussStormSingle(loadname,...,
                    appData.localizationThresh, appData.localizationWindow,...,
                    appData.testDetectionFrame);
                
            end
            
            if ~any(isnan(pos))
                hold all
                plot(pos(:,1),pos(:,2),'ro','MarkerSize',10);
                xlabel('x [pixels]');
                ylabel('y [pixels]');
                axis image;
            end
            
        else isfield(appData,'testFile')
            appData = rmfield(appData,'testFile');
        end
        
        
    case 'testDetectionFrame_Callback'
        
        testDetectionFrame = str2double(get(handles.testDetectionFrame,'String'));
        
        
        if isnan(testDetectionFrame)
            set(hObject, 'String', 1);
            errordlg('Input must be a number','Error');
        end
        
        appData.testDetectionFrame = testDetectionFrame;
        appData.sliderTestDetectionFrame = testDetectionFrame;
        set(handles.sliderTestDetectionFrame,'value',testDetectionFrame);
        
    case 'sliderTestDetectionFrame_Callback'
        L = round(get(handles.sliderTestDetectionFrame,'value'));
        set(handles.testDetectionFrame,'string',num2str(L));
        appData.testDetectionFrame = L;
        
        %From the Brightfield
        if ~isfield(appData,'testFile')
            [movieFilename, moviePathname] = ...
                uigetfile('*.*', 'MatFile data:','MultiSelect', 'off');
            
            appData.testFile.movieFilename = movieFilename;
            appData.testFile.moviePathname = moviePathname;
        else
            movieFilename = appData.testFile.movieFilename;
            moviePathname = appData.testFile.moviePathname;
        end
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
            ImageInfo = fits_read_header([moviePathname movieFilename]);
            imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
            
            nFrames = ImageInfo.NAXIS3;
            dataIn=ImageStack([moviePathname movieFilename], imageLim);
            
            singleImage = flipud(double(getFrame(dataIn,appData.testDetectionFrame)));
            
            if appData.holdFigureCheckbox;
                figure(appData.figureHandle);
                hold all;
            else
                figure;
            end
            
            imshow(singleImage,[min(min(singleImage)) max(max(singleImage))]);
            
            %From the Localisation
            if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
                
                
                loadname  = [moviePathname movieFilename];
                
                [pos(:,1),pos(:,2)] = gaussStormSingle(loadname,...,
                    appData.localizationThresh, appData.localizationWindow,...,
                    appData.testDetectionFrame);
                
            end
            
            if ~any(isnan(pos))
                hold all
                plot(pos(:,1),pos(:,2),'ro','MarkerSize',10);
                xlabel('x [pixels]');
                ylabel('y [pixels]');
                axis image;
            end
            
        else isfield(appData,'testFile')
            appData = rmfield(appData,'testFile');
        end
        
        
        
    case 'transformGreenBlue_Callback'
        
        transformLocalizations;
        
    case 'coLocalisation_Callback'
        % mod.pos.out Are the files from the modified localisations
        if ~isfield(appData,'data')
            [appData.dataFilename, appData.dataPathname] = uigetfile('*mod.pos.out', 'MatFile data:','MultiSelect', 'off');
            
            if ~(isnumeric(appData.dataFilename)&&appData.dataFilename==0) %check the user has not pressed cancel
                newData = importdata([appData.dataPathname appData.dataFilename]);
                if isstruct(newData)
                    appData.data = newData.data;
                else
                    appData.data = newData;
                end
                clear newData;
            end
        end
        
        if isfield(appData,'data')
            % Position information of the replication fork
            pos(:,1) = appData.data(:,2);
            pos(:,2) = appData.data(:,3);
        end
        
        % .ROIData contains the ROI of the cells and their tracks in side
        % the ROI
        if ~exist('ROIData')
            [ROIDataFilename, ROIDataPathname] = ...
                uigetfile('*.ROIData', 'MatFile data:','MultiSelect', 'off');
            if ~(isnumeric(ROIDataFilename)&&ROIDataFilename==0) %check the user has not pressed cancel
                
                info = whos('ROIDataFilename');
                if strcmp(info.class,'char')
                    appData.nFiles = 1;
                else
                    appData.nFiles = numel(ROIDataFilename);
                end
                
                for ii = 1:appData.nFiles
                    
                    if appData.nFiles == 1
                        loadname  = [ROIDataPathname ROIDataFilename];
                    else
                        loadname = [ROIDataPathname ROIDataFilename{1,ii}];
                    end
                    
                    ROIData = importdata(loadname);
                    
                    coLocalisationDistance_v2( pos, ROIData, ROIDataFilename, appData )
                end
            end
        end
        
        
        
    case 'editLocalisations_Callback'
        

        if ~isfield(appData,'testFile')
            [movieFilename, moviePathname] = ...
                uigetfile('*.fits', 'MatFile data:','MultiSelect', 'off');
            
            appData.testFile.movieFilename = movieFilename;
            appData.testFile.moviePathname = moviePathname;
        else
            movieFilename = appData.testFile.movieFilename;
            moviePathname = appData.testFile.moviePathname;
        end
        
        if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
            ImageInfo = fits_read_header([moviePathname movieFilename]);
            imageLim = [1 ImageInfo.NAXIS1 1 ImageInfo.NAXIS2];
            
            %nFrames = ImageInfo.NAXIS3;
            dataIn=ImageStack([moviePathname movieFilename], imageLim);
            
            singleImage = flipud(double(getFrame(dataIn,appData.testDetectionFrame)));
            
            if appData.holdFigureCheckbox;
                fig = figure(appData.figureHandle);
                hold all;
            else
                fig = figure;
            end
            
            imshow(singleImage,[min(min(singleImage)) max(max(singleImage))]);
            

            if ~(isnumeric(movieFilename)&&movieFilename==0) %check the user has not pressed cancel
                
                loadname  = [moviePathname movieFilename];
                
                [pos(:,1),pos(:,2)] = gaussStormSingle(loadname,...,
                    appData.localizationThresh, appData.localizationWindow,...,
                    appData.testDetectionFrame);
                
            end
            
            newPos = addRemoveLocalisations(pos,fig,loadname,...,
                appData.localizationThresh,...,
                appData.testDetectionFrame);
            
            [~,~] = gaussStormSingleParticle(loadname,...,
                appData.localizationThresh, 3,...,
                appData.testDetectionFrame,newPos,1);
            
        else isfield(appData,'testFile')
            appData = rmfield(appData,'testFile');
        end
        
    case 'closeAll_Callback'
        
        figHandles = get(0,'Children');
        close(figHandles(figHandles~=handles.figure1))
        
    case 'clearData_Callback'
        
        if isfield(appData,'data')
            appData = rmfield(appData,'data');
        end
        if isfield(appData,'tracks')
            appData = rmfield(appData,'tracks');
        end
        if isfield(appData,'dataStruct')
            appData = rmfield(appData,'dataStruct');
        end
        if isfield(appData,'testFile')
            appData = rmfield(appData,'testFile');
        end
        set(handles.results,'String','');
        
    case 'histogram2dOverlay_Callback'
        
        histogram2DOverlay

    case 'temporaryFunction_Callback'       
                 
        
            appData.localizationWindow = 1;
            [tracksFilename, tracksPathname] = ...
            uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
        
        if ~(isnumeric(tracksFilename)&&tracksFilename==0) %check the user has not pressed cancel
            
            info = whos('tracksFilename');
            if strcmp(info.class,'char')
                appData.nFiles = 1;
            else
                appData.nFiles = numel(tracksFilename);
            end
            
            diffusionFraction = zeros(appData.nFiles,1);
            D1 = zeros(appData.nFiles,1);
            D2 = zeros(appData.nFiles,1);
            
            for ii = 1:appData.nFiles
                
                if appData.nFiles == 1
                    loadname  = [tracksPathname tracksFilename];
                else
                    loadname = [tracksPathname tracksFilename{1,ii}];
                end
                
                newData = importdata(loadname);
                if isstruct(newData)
                    appData.tracks = newData.tracks;
                else
                    appData.tracks = newData;
                end
                clear newData;
                
               
                    boundTracks = twoSpeciesMSD2Threshold_forTraces(appData.tracks, appData);
                
                
            end
            
            assignin('base', 'MSDThreshResults', [diffusionFraction, D1, D2]);
            
        end
        
        % Track the bound tracks
         
        boundTracksMeans = [];
        molecInds = unique(boundTracks(:,4));
        altAppData = appData;
        
        for tNum = 1:length(molecInds) % for each track
    
            xx = find(boundTracks(:,4)==molecInds(tNum));
            boundTracksMeans = [boundTracksMeans; mean(boundTracks(xx,1)), mean(boundTracks(xx,2)),...
                boundTracks(xx(1),3), boundTracks(xx(1),4), boundTracks(xx(end),3)];
                        
        end
        
        altAppData.trackParams.mem = appData.trackParams.mem;
        altAppData.trackParams.maxDisp = 5;
        longTracks = trackWithDummy(boundTracks(:,1:3),altAppData.trackParams);
        save([loadname '.OnlyBound.tracks'], 'boundTracks');
        

        
        
%          appData.localizationWindow = 1;
%             [tracksFilename, tracksPathname] = ...
%             uigetfile('*.tracks', 'MatFile data:','MultiSelect', 'on');
%         
%         if ~(isnumeric(tracksFilename)&&tracksFilename==0) %check the user has not pressed cancel
%             
%             info = whos('tracksFilename');
%             if strcmp(info.class,'char')
%                 appData.nFiles = 1;
%             else
%                 appData.nFiles = numel(tracksFilename);
%             end
%             
%             diffusionFraction = zeros(appData.nFiles,1);
%             D1 = zeros(appData.nFiles,1);
%             D2 = zeros(appData.nFiles,1);
%             
%             for ii = 1:appData.nFiles
%                 
%                 if appData.nFiles == 1
%                     loadname  = [tracksPathname tracksFilename];
%                 else
%                     loadname = [tracksPathname tracksFilename{1,ii}];
%                 end
%                 
%                 newData = importdata(loadname);
%                 if isstruct(newData)
%                     appData.tracks = newData.tracks;
%                 else
%                     appData.tracks = newData;
%                 end
%                 clear newData;
%                 
%                
%                     boundTracks = twoSpeciesMSD2Threshold_forTraces(appData.tracks, appData);
%                 
%                 
%             end
%             
%             assignin('base', 'MSDThreshResults', [diffusionFraction, D1, D2]);
%             
%         end
%          [movieFilename, moviePathname] =...
%             uigetfile('*.fits', 'MatFile data:','MultiSelect', 'on');
%             
%             info = whos('movieFilename');
%             
%             for ii = 1:appData.nFiles
%                 
%                 if appData.nFiles == 1
%                     loadname  = [moviePathname movieFilename];
%                 else
%                     loadname = [moviePathname movieFilename{1,ii}];
%                 end
%                 
%                 gaussStorm_forTraces(loadname,...
%                     appData.localizationThresh, appData.localizationWindow,...
%                     boundTracks);
%                 
%             end
%             
    
end
%update handles
setappdata(handles.figure1, 'appData', appData);
guidata(handles.figure1, handles);



