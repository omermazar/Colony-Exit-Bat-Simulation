function varargout = BatGUI1xxxx(varargin)
% BATGUI1XXXX MATLAB code for BatGUI1xxxx.fig
%      BATGUI1XXXX, by itself, creates a new BATGUI1XXXX or raises the existing
%      singleton*.
%
%      H = BATGUI1XXXX returns the handle to a new BATGUI1XXXX or the handle to
%      the existing singleton*.
%
%      BATGUI1XXXX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATGUI1XXXX.M with the given input arguments.
%
%      BATGUI1XXXX('Property','Value',...) creates a new BATGUI1XXXX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BatGUI1xxxx_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BatGUI1xxxx_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BatGUI1xxxx

% Last Modified by GUIDE v2.5 30-Jun-2022 15:45:22

% Begin initialization code - DO NOT EDIT

% Changes By Omer 07/10/2017
%%%%% Changes By Omer Nov2021 Restart Values and Acoustics


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BatGUI1xxxx_OpeningFcn, ...
                   'gui_OutputFcn',  @BatGUI1xxxx_OutputFcn, ...
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


% --- Executes just before BatGUI1xxxx is made visible.
function BatGUI1xxxx_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BatGUI1xxxx (see VARARGIN)

%%%%%%%%
%% Defaualt parameters for The Simulation
% FileNameOfDefaultParams = 'DATA\DefaultParameters.mat'; 
FileNameOfDefaultParams =  'DATA\DefaultParamsTable.xlsx';
[handles.SimParams, handles.BatFlightParams, handles.BatSonarParams, ...
    handles.TerrainParams, handles.PreyFlightParams] = InitiateWithDefaultParams(FileNameOfDefaultParams);

%     load('C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\DATA\DefaultParameters');
%     % DefaultParameters is a MAT-file with the parameters neede to run the
%     % simulation
%     handles.SimParams = Params.SimParams;
%     handles.BatFlightParams = Params.BatFlight;
%     handles.BatSonarParams = Params.BatSonar;
%     handles.TerrainParams = Params.Terrain;
%     handles.PreyFlightParams = Params.PreyFlight;

%%% update the GUI with the default parameters %%%
% SimParams
set(handles.textDataAvailabilty,'String','DATA: NOT Available');
set(handles.editSampleTime,'String',handles.SimParams.SampleTime);
set(handles.editXYRes,'String',handles.SimParams.xyResolution);
set(handles.editSimulationTime,'String',handles.SimParams.SimulationTime);
set(handles.checkboxAvoidObs,'Value',handles.BatFlightParams.AvoidObsFlag);
set(handles.checkboxIsPrey,'Value',handles.SimParams.IsPreyFlag);
set(handles.editPreysNumber,'Value',handles.SimParams.TotalPreysNumber);
set(handles.editBatsNumber,'Value',handles.SimParams.TotalBatsNumber);

set(handles.popupmenu_TestMode, 'Value', find(strcmp(handles.popupmenu_TestMode.String, handles.SimParams.TestMode)) );

set(handles.editXroom,'String',handles.TerrainParams.Xmax);
set(handles.editYroom,'String',handles.TerrainParams.Ymax);
set(handles.editZroom,'String',handles.TerrainParams.Zmax);
set(handles.edit_swarm_start_pos_std, 'String', handles.SimParams.swarm_initial_ybat_std);

% set(handles.popupmenuEnvironment,'Value',2); % Empty Room
[handles.Terrain, handles.TerrainParams.Xmax,...
    handles.TerrainParams.Ymax, handles.TerrainParams.Zmax,...
    handles.TerrainParams.Xmin, handles.TerrainParams.Ymin, handles.TerrainParams.Zmin] = ...
    BuildEnvironment(handles.TerrainParams, handles.SimParams.xyResolution,handles.TerrainParams.Type);
set(handles.editZroom,'String',handles.TerrainParams.Zmax);
set(handles.editXroom,'String',handles.TerrainParams.Xmax);
set(handles.editYroom,'String',handles.TerrainParams.Ymax);
set(handles.popupmenuEnvironment, 'Value',find(strcmp(handles.popupmenuEnvironment.String, handles.TerrainParams.Type)));
 
if strcmp(handles.TerrainParams.Type, 'Free Space')
    set(handles.editXroom,'Enable', 'off');
    set(handles.editYroom,'Enable', 'off');
    set(handles.editZroom,'Enable', 'off');
end % if
  
set(handles.editBatsNumber,  'String',handles.SimParams.TotalBatsNumber);
set(handles.editPreysNumber, 'String',handles.SimParams.TotalPreysNumber);

% Bat Behavior
set(handles.checkboxDetectPrey,    'Value',handles.SimParams.DetectPrey);
set(handles.checkboxDetectClutter, 'Value',handles.SimParams.DetectObs);
set(handles.checkboxDetectConsps,  'Value',handles.SimParams.DetectConsps);

% BatFlighParams
set(handles.editMaxVelocity,    'String',handles.BatFlightParams.MaxVelocity);
set(handles.editNominalVelocity,'String',handles.BatFlightParams.NominalVelocity);
set(handles.editMaxAccel,       'String',handles.BatFlightParams.MaxAccelaration);
set(handles.editNominalAccel,   'String',handles.BatFlightParams.NominalAccel);
% set(handles.checkboxPreyHunting,'Value',handles.BatFlightParams.PreyHuntingFlag);            
set(handles.popupmenuFlightType,'Value',find(strcmp(handles.popupmenuFlightType.String, handles.BatFlightParams.FlightTypeFlag)) );
 

% BatSonarParams
set(handles.editBeamWidth,'String',handles.BatSonarParams.BatBeamWidth);
set(handles.editDetectRange,'String',handles.BatSonarParams.BatDetectionRange);
set(handles.editTransmitPower,'String',handles.BatSonarParams.PulsePower);
set(handles.editPRFMin,'String',handles.BatSonarParams.BatNominalPRF);
set(handles.editPRFMax,'String',handles.BatSonarParams.BatMaxPRF);
set(handles.editPulseTimeMin,'String',handles.BatSonarParams.PulseTimeShort);
set(handles.editPulseTimeMax,'String',handles.BatSonarParams.PulseTimeLong);
set(handles.editCenterFreq,'String',handles.BatSonarParams.CenterFreq);
set(handles.edit58_PulseDetectionTH,'String',handles.BatSonarParams.PulseDetectionTH);
% set(handles.editPulsePower,'String',handles.BatSonarParams.PulsePower);
set(handles.edit59_NoiseLeveldB,'String',handles.BatSonarParams.NoiseLeveldB);
set(handles.edit42TargetWingLength,'String',handles.BatSonarParams.TargetWingLength);
set(handles.checkbox_MissDetectionProbFlag,'Value',handles.BatSonarParams.MissDetectionProbFlag);
set(handles.edit43_ProbToDetect,'String',handles.BatSonarParams.ProbToDetect);
% set(handles.edit45DetectionTHdB,'String',handles.BatSonarParams.PulseDetectionTH);
% set(handles.edit46_MemorySize,'String',handles.BatFlightParams.HuntedMemorySize);
% set(handles.edit47_MemoryKHunted,'String',handles.BatFlightParams.NumOfSamePreyForAppraoach);
set(handles.edit48_SIRneededTH,'String',handles.BatSonarParams.Signal2InterferenceRatio);
set(handles.edit49_FwdMaskingTH,'String',handles.BatSonarParams.MaskingFwdSIRth);
set(handles.edit50_FwdMakingTime,'String',handles.BatSonarParams.MaskingFwdTime);
set(handles.edit51_MinPulsePower,'String',handles.BatSonarParams.PulseMinPower);
set(handles.edit52_Dist2Buzz,'String',handles.BatFlightParams.DistanceToBuzz);
set(handles.edit53_Dist2Appraoch,'String',handles.BatFlightParams.DistanceFromPreyToApproach);
set(handles.edit54_BatMinVelocity,'String',handles.BatFlightParams.MinVelocity);
set(handles.checkbox29_LocalizationErrFlag, 'Value', handles.BatSonarParams.LocalizationErrFlag);

set(handles.checkbox37_MaskingByConsps, 'Value', handles.SimParams.MaskingByConsps);
set(handles.checkbox38_MaskingByClutter, 'Value', handles.SimParams.MaskingByClutter);
set(handles.radiobutton28_Acoustics, 'Value', handles.SimParams.AcousticsCalcFlag);

% handles.popupmenuReceiverType.Value =  find(strcmp(handles.popupmenuReceiverType.String, handles.BatSonarParams.ReceiverType) );
handles.popupmenu_receiverType.Value =  find(strcmp(handles.popupmenu_receiverType.String, handles.BatSonarParams.ReceiverType) );
% set(handles.popupmenuReceiverType,'String',handles.BatSonarParams.ReceiverType);


% PreyFlightParams
set(handles.editPreyMaxVelocity,'String',handles.PreyFlightParams.PmaxVelocity);
set(handles.editPreyNominalVelocity,'String',handles.PreyFlightParams.PNominalVelocity);
set(handles.editPreyMaxAccelaration,'String',handles.PreyFlightParams.PMaxAccelaration);
set(handles.editPreyNominalAccel,'String',handles.PreyFlightParams.PNominalAccel);
set(handles.popupmenuPreyFlightType, 'Value', find(strcmp(handles.popupmenuPreyFlightType.String, handles.PreyFlightParams.PFlightTypeFlag)) );    
set(handles.checkbox_JamMothFlg, 'Value', handles.PreyFlightParams.React2BatsFlag);
set(handles.edit_jamMothTxLvl, 'String', handles.PreyFlightParams.JamTxLvl);
% multi Prey an Bats
handles.DataToPlot.PlotAllPreysFlag = handles.checkbox19AllPreysPlot.Value;
handles.DataToPlot.PlotAllBatsFlag = handles.checkboxAllBatsPlot.Value;
handles.edit18BatNumToPlot.Value = 1;
handles.DataToPlot.BatNumToPlot = handles.edit18BatNumToPlot.Value;
handles.edit40PreyNumToPlot.Value = 1;
handles.DataToPlot.PreyNumToPlot = handles.edit40PreyNumToPlot.Value;
% X-Y Plot Paramsha
 handles.DataToPlot.PreyPositionFlag = handles.radiobuttonPeryPosFLag.Value;
 handles.DataToPlot.TerrainFlag = handles.radiobuttonTerrain.Value; 
 handles.DataToPlot.BatPosFlag = handles.radiobuttonBatPos.Value;
 handles.DataToPlot.SonarPosFlag = handles.radiobuttonSonarPos.Value;
 handles.DataToPlot.ManPosFlag = handles.radiobuttonManuevresPos.Value;
 handles.DataToPlot.TerrainFindFlag = handles.radiobuttonTerrainFind.Value;
 handles.DataToPlot.EndTimeFilter = handles.editEndTime.Value;
 handles.DataToPlot.StartTimeFilter = handles.editStartTime.Value;
 handles.DataToPlot.TimeFilterFlag = handles.checkboxTimeFilter.Value;
 handles.DataToPlot.CatchesPosFlag = handles.radiobuttonCatchesPositions.Value;
 handles.DataToPlot.JamPreyPosFlag = handles.radiobutton_PreyJam.Value;
 set(handles.checkboxHoldOn,'Value',1);
 handles.handleXYBatPosLine =[];
 % Time Plot Params
 handles.DataToPlot.HuntingManuevers =  handles.radiobuttonHuntMan.Value;
 handles.DataToPlot.FlightVelocity =  handles.radiobuttonFlightVelocity.Value; 
 handles.DataToPlot.SonarTimes =  handles.radiobuttonSonarTimes.Value;
 handles.DataToPlot.PulseDuration=  handles.radiobuttonPulseDur.Value;
 handles.DataToPlot.SonarPulses =  handles.radiobuttonRadialACcel.Value;
 handles.DataToPlot.PRFVec = handles.radiobuttonPRFVec.Value;
 handles.DataToPlot.Bat2PreyDistFlag = handles.radiobuttonBat2PreyDistFlag.Value ;
 handles.DataToPlot.Bat2PreyAngleFlag =  handles.radiobutton15Bat2PreyAngleFlag.Value;
 handles.DataToPlot.PreyFindsFlag = handles.radiobuttonPreyFindsFlag.Value;
 handles.DataToPlot.BatAllPulsesFlag = handles.radiobuttonAllSonarPulses.Value;
 handles.DataToPlot.ObsEchosFlag = handles.radiobuttonObsEchos.Value;
 handles.DataToPlot.InterferencePlotFlag =  handles.radiobutton20_PlotInterferenceFlag.Value;
 handles.DataToPlot.PreyEchosFlag = handles.radiobuttonPreyEchos.Value;
 handles.DataToPlot.FlightVelocity = handles.radiobuttonFlightVelocity.Value;
 set(handles.checkboxTimeHoldOn,'Value',1); 
 hold(handles.TimePlot, 'on')

 % Animation
handles.SaveAminationFlag = 0;
handles.AnimationTimeRate = 75;
handles.AnimmateCurrentOnly = 0;
handles.DataToPlot.AnimateStartTime = 0; % seconds
handles.DataToPlot.AnimateEndTime = handles.SimParams.SimulationTime; %seconds
guidata(hObject,handles); % save the changes to the structure 

% Data to Analyze
handles.DataToAnalyze.BatPos = [];
handles.DataToAnalyze.BatDirection = [];
handles.DataToAnalyze.BatDirectionChange = [];
handles.DataToAnalyze.SonarTimesVec = [];
handles.DataToAnalyze.ObsManueverFlag = [];
handles.DataToAnalyze.HuntingFlag = [];
handles.DataToAnalyze.ObsticlesFound = [];
handles.DataToAnalyze.PRFVec = [];
handles.DataToAnalyze.PulseWidthVec = [];
handles.DataToAnalyze.CatchPreyPos = []; 
handles.ZoomHandle = zoom;
guidata(hObject,handles); % save the changes to the structure 


% Choose default command line output for BatGUI1xxxx
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);






% --- Outputs from this function are returned to the command line.
function varargout = BatGUI1xxxx_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edit_menu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_flose_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_flose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_saveAs_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_saveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function view_menu_Callback(hObject, eventdata, handles)
% hObject    handle to view_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonRunAnimation.
function pushbuttonRunAnimation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunAnimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% BatAnimation(handles.DataToAnalyze, handles.SimParams );

MyBatAnimation(handles.DataToAnalyze, handles.DataToPlot, handles.SaveAminationFlag, handles.Terrain, handles.AnimationTimeRate, handles.AnimmateCurrentOnly);


% --- Executes on button press in radiobuttonTerrain.
function radiobuttonTerrain_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonTerrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonTerrain
 handles.DataToPlot.TerrainFlag =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure  
        
% --- Executes on button press in radiobuttonBatPos.
function radiobuttonBatPos_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonBatPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.DataToPlot.BatPosFlag =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure  
        
% Hint: get(hObject,'Value') returns toggle state of radiobuttonBatPos


% --- Executes on button press in radiobuttonSonarPos.
function radiobuttonSonarPos_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonSonarPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.SonarPosFlag =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure  
        
% Hint: get(hObject,'Value') returns toggle state of radiobuttonSonarPos


% --- Executes on button press in radiobuttonManuevresPos.
function radiobuttonManuevresPos_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonManuevresPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles.DataToPlot.ManPosFlag =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure  
        
% Hint: get(hObject,'Value') returns toggle state of radiobuttonManuevresPos


% --- Executes on button press in radiobuttonTerrainFind.
function radiobuttonTerrainFind_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonTerrainFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.TerrainFindFlag =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure  
        
% Hint: get(hObject,'Value') returns toggle state of radiobuttonTerrainFind


% --- Executes on button press in radiobuttonFlightVelocity.
function radiobuttonFlightVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonFlightVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.FlightVelocity =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonFlightVelocity


% --- Executes on button press in radiobuttonSonarTimes.
function radiobuttonSonarTimes_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonSonarTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.SonarTimes =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonSonarTimes



function editSampleTime_Callback(hObject, eventdata, handles)
% hObject    handle to editSimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSimulationTime as text
%        str2double(get(hObject,'String')) returns contents of editSimulationTime as a double
ed1 = get(hObject,'String');
handles.SimParams.SampleTime = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editSampleTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% hObject.String = 5;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editResolution_Callback(hObject, eventdata, handles)
% hObject    handle to editResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editResolution as text
%        str2double(get(hObject,'String')) returns contents of editResolution as a double
ed1 = get(hObject,'String');
handles.SimParams.xyResolution = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkboxAvoidObs.
function checkboxAvoidObs_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAvoidObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAvoidObs
ChBx1 = get(hObject,'Value');
handles.BatFlightParams.AvoidObsFlag = ChBx1;
guidata(hObject,handles); % save the changes to the structure 


% function editResolution_Callback(hObject, eventdata, handles)
% hObject    handle to editResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editResolution as text
%        str2double(get(hObject,'String')) returns contents of editResolution as a double
% ed1 = get(hObject,'String');
% handles.SimParams.xyResolution = str2double(ed1);



function editSimulationTime_Callback(hObject, eventdata, handles)
% hObject    handle to editResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editResolution as text
%        str2double(get(hObject,'String')) returns contents of editResolution as a double
ed1 = get(hObject,'String');
handles.SimParams.SimulationTime  = str2double(ed1);
handles.DataToPlot.AnimateEndTime = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editSimulationTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function editTotalTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTotalTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editXYRes_Callback(hObject, eventdata, handles)
% hObject    handle to editXYRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXYRes as text
%        str2double(get(hObject,'String')) returns contents of editXYRes as a double
ed1 = get(hObject,'String');
handles.SimParams.xyResolution = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editXYRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXYRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuEnvironment.
function popupmenuEnvironment_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEnvironment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEnvironment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEnvironment
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
handles.TerrainParams.Type = popupValue;
[handles.Terrain, handles.TerrainParams.Xmax,...
    handles.TerrainParams.Ymax, handles.TerrainParams.Zmax,...
    handles.TerrainParams.Xmin, handles.TerrainParams.Ymin, handles.TerrainParams.Zmin] =...
    BuildEnvironment(handles.TerrainParams, handles.SimParams.xyResolution,...
    handles.TerrainParams.Type);
% Control t size room buttons
% if strcmp(handles.TerrainParams.Type, 'Free Space')
%     set(handles.editXroom,'Enable', 'off');
%     set(handles.editYroom,'Enable', 'off');
%     set(handles.editZroom,'Enable', 'off');
% else
%    set(handles.editXroom,'Enable', 'on');
%    set(handles.editYroom,'Enable', 'on');
%    set(handles.editZroom,'Enable', 'on'); 
% end % if


handles.TerrainParams.Type = popupValue;
guidata(hObject, handles);
                       
% end % popupmenuEnvironment_Callback

% --- Executes during object creation, after setting all properties.
function popupmenuEnvironment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuEnvironment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editXroom_Callback(hObject, eventdata, handles)
% hObject    handle to editXroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXroom as text
%        str2double(get(hObject,'String')) returns contents of editXroom as a double
ed1 = get(hObject,'String');
handles.TerrainParams.Xmax = str2double(ed1);
[handles.Terrain, handles.TerrainParams.Xmax,...
    handles.TerrainParams.Ymax, handles.TerrainParams.Zmax,...
    handles.TerrainParams.Xmin, handles.TerrainParams.Ymin, handles.TerrainParams.Zmin] =...
    BuildEnvironment(handles.TerrainParams, handles.SimParams.xyResolution,...
    handles.TerrainParams.Type);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editXroom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editZroom_Callback(hObject, eventdata, handles)
% hObject    handle to editZroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editZroom as text
%        str2double(get(hObject,'String')) returns contents of editZroom as a double
ed1 = get(hObject,'String');
handles.TerrainParams.Zmax = str2double(ed1);
[handles.Terrain, handles.TerrainParams.Xmax,...
    handles.TerrainParams.Ymax, handles.TerrainParams.Zmax,...
    handles.TerrainParams.Xmin, handles.TerrainParams.Ymin, handles.TerrainParams.Zmin] =...
    BuildEnvironment(handles.TerrainParams, handles.SimParams.xyResolution,...
    handles.TerrainParams.Type);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editZroom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editYroom_Callback(hObject, eventdata, handles)
% hObject    handle to editYroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editYroom as text
%        str2double(get(hObject,'String')) returns contents of editYroom as a double
ed1 = get(hObject,'String');
handles.TerrainParams.Ymax = str2double(ed1);
[handles.Terrain, handles.TerrainParams.Xmax,...
    handles.TerrainParams.Ymax, handles.TerrainParams.Zmax,...
    handles.TerrainParams.Xmin, handles.TerrainParams.Ymin, handles.TerrainParams.Zmin] =...
    BuildEnvironment(handles.TerrainParams, handles.SimParams.xyResolution,...
    handles.TerrainParams.Type);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editYroom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editYroom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonPlotxy.
function pushbuttonPlotxy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ZoomHandle,'ActionPostCallback',...
    {@MyZoom1,handles.xyPlot, handles.TimePlot, handles.DataToPlot, handles.DataToAnalyze,...
     handles.SimParams, handles.TerrainParams, handles.Terrain});
ZoomFlag= 0;
IDn = [];
h1 = handles.xyPlot; % the axes of the current fig
c1 = get(handles.checkboxHoldOn,'Value'); % the flag wether to hold on
CheckAndHoldOnOff(h1,c1);

if handles.checkbox41_NewFig.Value == 0
    axHandle = handles.xyPlot;
else
    figure(); hold on;
    axHandle = gca;
end % if handles.checkbox41_NewFig.Value == 1

MyXYPlot(axHandle,handles.DataToAnalyze,handles.DataToPlot,...
    handles.SimParams, handles.TerrainParams, handles.Terrain, ZoomFlag, IDn);
guidata(hObject,handles);


% --- Executes on mouse press over axes background.
function xyPlot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to xyPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function TimePlot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to TimePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function xyPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xyPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate xyPlot



function editMaxAccel_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxAccel as text
%        str2double(get(hObject,'String')) returns contents of editMaxAccel as a double
ed1 = get(hObject,'String');
handles.BatFlightParams.MaxAccelaration = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editMaxAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNominalVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to editNominalVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNominalVelocity as text
%        str2double(get(hObject,'String')) returns contents of editNominalVelocity as a double
ed1 = get(hObject,'String');
handles.BatFlightParams.NominalVelocity = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editNominalVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNominalVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ed1 = get(hObject,'String');
handles.BatFlightParams.MaxVelocity = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% Hints: get(hObject,'String') returns contents of editMaxVelocity as text
%        str2double(get(hObject,'String')) returns contents of editMaxVelocity as a double


% --- Executes during object creation, after setting all properties.
function editMaxVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxAvoidObs.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAvoidObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAvoidObs
ChBx1 = get(hObject,'Value');
handles.SimParams.AvoidObsticleFlag = ChBx1;
guidata(hObject,handles) % save the changes to the structure 



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to editSimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSimulationTime as text
%        str2double(get(hObject,'String')) returns contents of editSimulationTime as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSimulationTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to editSampleTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSampleTime as text
%        str2double(get(hObject,'String')) returns contents of editSampleTime as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSampleTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to editXYRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXYRes as text
%        str2double(get(hObject,'String')) returns contents of editXYRes as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXYRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFlightType.
function popupmenuFlightType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFlightType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
handles.BatFlightParams.FlightTypeFlag = popupValue;
guidata(hObject,handles);% save the struct

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuFlightType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFlightType


% --- Executes during object creation, after setting all properties.
function popupmenuFlightType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFlightType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNominalAccel_Callback(hObject, eventdata, handles)
% hObject    handle to editNominalAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNominalAccel as text
%        str2double(get(hObject,'String')) returns contents of editNominalAccel as a double
ed1 = get(hObject,'String');
handles.BatFlightParams.NominalAccel = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editNominalAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNominalAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uibuttongroup2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxVarianceSonar.
function checkboxVarianceSonar_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxVarianceSonar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxVarianceSonar



function editBeamWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editBeamWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBeamWidth as text
%        str2double(get(hObject,'String')) returns contents of editBeamWidth as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.BatBeamWidth = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editBeamWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBeamWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDetectRange_Callback(hObject, eventdata, handles)
% hObject    handle to editDetectRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDetectRange as text
%        str2double(get(hObject,'String')) returns contents of editDetectRange as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.BatDetectionRange = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editDetectRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDetectRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTransmitPower_Callback(hObject, eventdata, handles)
% hObject    handle to editTransmitPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTransmitPower as text
%        str2double(get(hObject,'String')) returns contents of editTransmitPower as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.PulsePower = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editTransmitPower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTransmitPower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPRFMin_Callback(hObject, eventdata, handles)
% hObject    handle to editPRFMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPRFMin as text
%        str2double(get(hObject,'String')) returns contents of editPRFMin as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.BatNominalPRF = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPRFMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPRFMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPRFMax_Callback(hObject, eventdata, handles)
% hObject    handle to editPRFMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPRFMax as text
%        str2double(get(hObject,'String')) returns contents of editPRFMax as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.BatMaxPRF = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPRFMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPRFMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCenterFreq_Callback(hObject, eventdata, handles)
% hObject    handle to editCenterFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCenterFreq as text
%        str2double(get(hObject,'String')) returns contents of editCenterFreq as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.CenterFreq = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editCenterFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCenterFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% function edit23PulseDetectionTH_Callback(hObject, eventdata, handles)
% % hObject    handle to edit23PulseDetectionTH (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit23PulseDetectionTH as text
% %        str2double(get(hObject,'String')) returns contents of edit23PulseDetectionTH as a double
% ed1 = get(hObject,'String');
% handles.BatSonarParams.PulseDetectionTH = str2double(ed1);
% guidata(hObject,handles); % save the changes to the structure 
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit23PulseDetectionTH_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit23PulseDetectionTH (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPulseTimeMin_Callback(hObject, eventdata, handles)
% hObject    handle to editPulseTimeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPulseTimeMin as text
%        str2double(get(hObject,'String')) returns contents of editPulseTimeMin as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.PulseTimeShort = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPulseTimeMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPulseTimeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPulseTimeMax_Callback(hObject, eventdata, handles)
% hObject    handle to editPulseTimeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPulseTimeMax as text
%        str2double(get(hObject,'String')) returns contents of editPulseTimeMax as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.PulseTimeLong = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 



% --- Executes during object creation, after setting all properties.
function editPulseTimeMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPulseTimeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadDefaultParams.
function pushbuttonLoadDefaultParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadDefaultParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BatGUI1xxxx_OpeningFcn();

% --- Executes on button press in pushbuttonRunSim.
function pushbuttonRunSim_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set(hObject,'String','WAIT');
set(hObject,'Enable','off');
AllParameters.SimParams = handles.SimParams;
AllParameters.BatFlightParams = handles.BatFlightParams;
AllParameters.BatSonarParams = handles.BatSonarParams;
AllParameters.TerrainParams = handles.TerrainParams;
AllParameters.PreyFlightParams = handles.PreyFlightParams;
set(handles.textDataAvailabilty,'String','DATA on Progress'); 
handles.DataToAnalyze = BatFlightForGui(AllParameters, handles.Terrain);

set(hObject,'Enable','on');
set(hObject,'String','RUN');
set(handles.textDataAvailabilty,'String','DATA: from running Simulation'); % 

%%% Flight Summary and statitics for GUI
SummaryString1 = ['Number Of Cataches: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanNumberOfCatches), ...
    '   (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDNumberOfCatches,3), ' )'];
set(handles.textSummaryGroup1,'String',SummaryString1); % 

SummaryString2 = ['Number Of Prey Detections: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanNumberOfDetections,4), ...
    '   (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDNumberOfDetections,3), ' )'];
set(handles.textSummaryGroup2,'String',SummaryString2); % 

SummaryString3 = ['Approach Masking Ratio Hunted: ',...
    num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanApproach_MaskingRatioToHunted,2), ...
    ' (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDApproach_MaskingRatioToHunted,2), ' )'];
set(handles.textSummaryGroup3,'String',SummaryString3); % 

SummaryString4 = ['Buzz Masking Ratio Hunted: ', ...
    num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanBuzz_MaskingRatioToHunted,2), ...
    ' (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDBuzz_MaskingRatioToHunted,2), ' )'];
set(handles.textSummaryGroup4,'String',SummaryString4); % 

SummaryString5 = ['Search Masking Ratio Hunted: ',...
    num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanSearch_MaskingRatioToHunted,2) ...
    ' (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDSearch_MaskingRatioToHunted,2), ' )'];
set(handles.textSummaryGroup5,'String',SummaryString5); % 

% SummaryString6 = ['Interfernce while Search: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanInterferfenceWhileSearching,3)...
%     '   (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDInterferfenceWhileSearching,3), ' )'];
% set(handles.textSummaryGroup6,'String',SummaryString6); % 

SummaryString7 = ['Interfernce to Detection Ratio: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.TotalSearchAndApproachInterfernceRatio,3)];
set(handles.textSummaryGroup7,'String',SummaryString7); % 

SummaryString8 = ['Interfernce to Hunted Ratio: ', ...
num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanInterferenceToHuntedRatio,3)];
set(handles.textSummaryGroup8,'String',SummaryString8); % 
% 
% SummaryString9 = ['Succesful Search Stages: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanNumberOfSuccesulSearchStages,3), ...
%     '   (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDNumberOfSuccesulSearchStages,3), ' )'];
% set(handles.textSummaryGroup9,'String',SummaryString9); % 

% SummaryString10 = ['Jamming  of Search Stages Ratio: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.JammedToSuccessfulSearchStagesRatio,3)];
% set(handles.textSummaryGroup10,'String',SummaryString10); %
% 
% SummaryString11 = ['Succesful Approach Stages: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.MeanNumberOfSuccesulApproachStages,3), ...
%     '   (' , num2str(handles.DataToAnalyze.FlightInterferceSummary.Statistics.STDNumberOfSuccesulApproachStages,3), ' )'];
% set(handles.textSummaryGroup11,'String',SummaryString11); %
% 
% SummaryString12 = ['Jamming  of Approach Stages Ratio: ', num2str(handles.DataToAnalyze.FlightInterferceSummary.JammedToSuccessfulApproachStagesRatio,3)];
% set(handles.textSummaryGroup12,'String',SummaryString12); %


% SummaryString9 = ['Interfernce to Detection Ratio: ', num2str(handles.DataToAnalyze.FilghtInterferceSummary.InterferenceToHuntedRatio,3)];
% set(handles.textSummaryGroup7,'String',SummaryString9); % 

guidata(hObject,handles); % save the changes to the structure 




% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbuttonRunSim.
function pushbuttonRunSim_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbuttonRunSim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonRunSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonSTOP.
function pushbuttonSTOP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSTOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return;

% --- Executes during object creation, after setting all properties.
function pushbuttonSTOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonSTOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkboxHoldOn.
function checkboxHoldOn_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxHoldOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxHoldOn
h1 = handles.xyPlot;
c1 = get(hObject,'Value');
CheckAndHoldOnOff(h1,c1);
% if c1
%     hold(h1,'on');
% else
%     hold(h1,'off');
% end
guidata(hObject, handles);

% --- Executes on button press in pushbuttonClearFigxy.
function pushbuttonClearFigxy_Callback(hObject, ~, handles)
% hObject    handle to pushbuttonClearFigxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h1 = handles.xyPlot;
cla(h1);


% --- Executes on button press in pushbuttonPlotTimeGraph.
function pushbuttonPlotTimeGraph_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotTimeGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ZoomHandle,'ActionPostCallback',...
    {@MyZoom1,handles.xyPlot, handles.TimePlot, handles.DataToPlot, handles.DataToAnalyze,...
     handles.SimParams, handles.TerrainParams, handles.Terrain});
ZoomFlag = 0;
IDn = [];
h1 = handles.xyPlot; % the axes of the current fig
c1 = get(handles.checkboxTimeHoldOn,'Value'); % the flag wether to hold on
CheckAndHoldOnOff(h1,c1);

if handles.checkbox41_NewFig.Value == 0
    axHandle = handles.TimePlot;
else
    figure();
    axHandle = gca;
end % if handles.checkbox41_NewFig.Value == 1

MyTimePlot(axHandle, handles.DataToAnalyze,handles.DataToPlot,...
    handles.SimParams, handles.TerrainParams, handles.Terrain, ZoomFlag, IDn) 

guidata(hObject,handles) % save the changes to the structure 

% --- Executes on button press in pushbuttonClearTimeFig.
function pushbuttonClearTimeFig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClearTimeFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h2 = handles.TimePlot;
cla(h2);

% --- Executes on button press in radiobuttonPulseDur.
function radiobuttonPulseDur_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonPulseDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.PulseDuration =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonPulseDur


% --- Executes on button press in radiobuttonHuntMan.
function radiobuttonHuntMan_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonHuntMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.HuntingManuevers =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonHuntMan


% --- Executes on button press in radiobuttonRadialACcel.
function radiobuttonRadialACcel_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonRadialACcel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.RadialAccel =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonRadialACcel


% --- Executes on button press in checkboxTimeHoldOn.
function checkboxTimeHoldOn_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTimeHoldOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h2 = handles.TimePlot;
c2 = get(hObject,'Value');
CheckAndHoldOnOff(h2,c2);

guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkboxTimeHoldOn


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in pushbuttonSendData.
function pushbuttonSendData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSendData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BatDATA = handles.DataToAnalyze;
assignin('base','BatDATA',BatDATA);


% --- Executes on button press in radiobuttonPRFVec.
function radiobuttonPRFVec_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonPRFVec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonPRFVec
handles.DataToPlot.PRFVec =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 


% --- Executes on button press in checkboxPreyHunting.
% function checkboxPreyHunting_Callback(hObject, eventdata, handles)
% % hObject    handle to checkboxPreyHunting (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% ChBx1 = get(hObject,'Value');
% handles.BatFlightParams.PreyHuntingFlag = ChBx1;
% guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkboxPreyHunting



%


% --- Executes on button press in checkboxTimeFilter.
function checkboxTimeFilter_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTimeFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTimeFilter
ChBx1 = get(hObject,'Value');
handles.DataToPlot.TimeFilterFlag = ChBx1;
guidata(hObject,handles) 


function editStartTime_Callback(hObject, eventdata, handles)
% hObject    handle to editStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStartTime as text
%        str2double(get(hObject,'String')) returns contents of editStartTime as a double
ed1 = get(hObject,'String');
handles.DataToPlot.StartTimeFilter = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure

% --- Executes during object creation, after setting all properties.
function editStartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEndTime_Callback(hObject, eventdata, handles)
% hObject    handle to editEndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEndTime as text
%        str2double(get(hObject,'String')) returns contents of editEndTime as a double
ed1 = get(hObject,'String');
handles.DataToPlot.EndTimeFilter = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure


% --- Executes during object creation, after setting all properties.
function editEndTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuPreyFlightType.
function popupmenuPreyFlightType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPreyFlightType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuPreyFlightType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuPreyFlightType
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
handles.PreyFlightParams.PFlightTypeFlag = popupValue;
guidata(hObject,handles);% save the struct


% --- Executes during object creation, after setting all properties.
function popupmenuPreyFlightType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPreyFlightType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPreyMaxVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to editPreyMaxVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreyMaxVelocity as text
%        str2double(get(hObject,'String')) returns contents of editPreyMaxVelocity as a double
ed1 = get(hObject,'String');
handles.PreyFlightParams.PMaxVelocity = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPreyMaxVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreyMaxVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPreyNominalVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to editPreyNominalVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreyNominalVelocity as text
%        str2double(get(hObject,'String')) returns contents of editPreyNominalVelocity as a double
ed1 = get(hObject,'String');
handles.PreyFlightParams.PNominalVelocity = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPreyNominalVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreyNominalVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPreyNominalAccel_Callback(hObject, eventdata, handles)
% hObject    handle to editPreyNominalAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreyNominalAccel as text
%        str2double(get(hObject,'String')) returns contents of editPreyNominalAccel as a double
ed1 = get(hObject,'String');
handles.PreyFlightParams.PNominalAccelaration = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editPreyNominalAccel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreyNominalAccel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadData.
function pushbuttonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [FileName,PathName,FilterIndex] = uigetfile('C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\DATA\*.mat');
[FileName,PathName,FilterIndex] = uigetfile('DATA\*.mat');

PathFileName = [PathName,FileName];
load(PathFileName);
LoadedBatDATA = load(PathFileName);
try
    handles.DataToAnalyze = LoadedBatDATA.BatDATA;
catch
    handles.DataToAnalyze = LoadedBatDATA.DataToAnalyze; % for Preleimnary results that are saves as (xxx).DataToAnalyze
end % try
set(handles.textDataAvailabilty,'String',['DATA: ', FileName]);
guidata(hObject,handles) % save the changes to the structure 

% % % %     


% --- Executes on button press in pushbuttonSaveData.
function pushbuttonSaveData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
BatDATA = handles.DataToAnalyze;

StrDate = datestr(datetime('now'));
SaveFileName = strrep(StrDate,' ','_');
SaveFileName = strrep(SaveFileName,':','');
SaveFileName = ['BatData_',SaveFileName,'.mat'];
FileFilter = '*.mat';
% [SaveFileName,PathName,FilterIndex] = uiputfile(['C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\DATA\', SaveFileName],FileFilter);
[SaveFileName,PathName,FilterIndex] = uiputfile(['DATA\', SaveFileName],FileFilter);
% PathFileName = [PathName,FileName];

% save(['C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\DATA\',SaveFileName], 'BatDATA');
save([PathName,SaveFileName], 'BatDATA');

% % % % % % % % % --- Executes during object creation, after setting all properties.
% % % % % % % % function editPreyMaxVelocity_CreateFcn(hObject, eventdata, handles)
% % % % % % % % % hObject    handle to editPreyMaxVelocity (see GCBO)
% % % % % % % % % eventdata  reserved - to be defined in a future version of MATLAB
% % % % % % % % % handles    empty - handles not created until after all CreateFcns called
% % % % % % % % 
% % % % % % % % % Hint: edit controls usually have a white background on Windows.
% % % % % % % % %       See ISPC and COMPUTER.
% % % % % % % % if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
% % % % % % % %     set(hObject,'BackgroundColor','white');
% % % % % % % % end



function editPreyMaxAccelaration_Callback(hObject, eventdata, handles)
% hObject    handle to editPreyMaxAccelaration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreyMaxAccelaration as text
%        str2double(get(hObject,'String')) returns contents of editPreyMaxAccelaration as a double
ed1 = get(hObject,'String');
handles.PreyFlightParams.PMaxAccelaration = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editPreyMaxAccelaration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreyMaxAccelaration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkboxIsPrey.
function checkboxIsPrey_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxIsPrey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxIsPrey
ChBx1 = get(hObject,'Value');
handles.SimParams.IsPreyFlag = ChBx1;
guidata(hObject,handles) % save the changes to the structure 


% --- Executes on button press in radiobuttonPeryPosFLag.
function radiobuttonPeryPosFLag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonPeryPosFLag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.PreyPositionFlag =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonPeryPosFLag


% --- Executes on button press in radiobuttonPreyFindsFlag.
function radiobuttonPreyFindsFlag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonPreyFindsFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.PreyFindsFlag =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonPreyFindsFlag


% --- Executes on button press in radiobuttonBat2PreyDistFlag.
function radiobuttonBat2PreyDistFlag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonBat2PreyDistFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.Bat2PreyDistFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobuttonBat2PreyDistFlag


% --- Executes on button press in radiobutton15Bat2PreyAngleFlag.
function radiobutton15Bat2PreyAngleFlag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton15Bat2PreyAngleFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.Bat2PreyAngleFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton15Bat2PreyAngleFlag


% --- Executes during object creation, after setting all properties.
function text26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text9.
function text9_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function CheckAndHoldOnOff(Axes, HoldFlag)
% this function decides weather the axes(Axes) should be hold on or off
% according to the flag HoldFlag (o = hold off, 1= hold on)

if HoldFlag
    hold(Axes,'on');
else
    hold(Axes,'off');
end

% --- Executes on button press in radiobuttonCatchesPositions.
function radiobuttonCatchesPositions_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCatchesPositions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.CatchesPosFlag =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure  
% Hint: get(hObject,'Value') returns toggle state of radiobuttonCatchesPositions


function [ SimParams, BatFlightParams, BatSonarParams,...
    TerrainParams, PreyFlightParams ] = InitiateWithDefaultParams(FilePathName) 
%this functios load the parmaters saved in FilePAthe name to handles
%Structure
%
%     load(FilePathName);
    Params = ReadDefaultParamsTable(FilePathName);
    % DefaultParameters is a MAT-file with the parameters neede to run the
    % simulation
    SimParams = Params.SimParams;
    BatFlightParams = Params.BatFlight;
    BatSonarParams = Params.BatSonar;
    TerrainParams = Params.Terrain;
    PreyFlightParams = Params.PreyFlight;
%     guidata(hObject,handles)


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkboxChirp.
function checkboxChirp_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxChirp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    ChBx1 = get(hObject,'Value');
    handles.BatSonarParams.ChirpFlag = ChBx1;
    guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkboxChirp



function editChirpSpan_Callback(hObject, eventdata, handles)
% hObject    handle to editChirpSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editChirpSpan as text
%        str2double(get(hObject,'String')) returns contents of editChirpSpan as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.ChirpSpan = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function editChirpSpan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editChirpSpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over text15.
function text15_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to text15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editPulsePower_Callback(hObject, eventdata, handles)
% hObject    handle to editPulsePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPulsePower as text
%        str2double(get(hObject,'String')) returns contents of editPulsePower as a double


% --- Executes during object creation, after setting all properties.
function editPulsePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPulsePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonSaveParameters.
function pushbuttonSaveParameters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SavePath = 'DATA\DefaultParameters';
% SavePath = 'C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\DATA\DefaultParameters';
Params.SimParams = handles.SimParams;
Params.BatFlight = handles.BatFlightParams;
Params.BatSonar = handles.BatSonarParams;
Params.Terrain= handles.TerrainParams;
Params.PreyFlight = handles.PreyFlightParams;
save(SavePath,'Params');


% --- Executes on button press in radiobuttonAllSonarPulses.
function radiobuttonAllSonarPulses_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonAllSonarPulses (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.BatAllPulsesFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobuttonAllSonarPulses


% --- Executes on button press in radiobuttonPreyEchos.
function radiobuttonPreyEchos_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonPreyEchos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.PreyEchosFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobuttonPreyEchos


% --- Executes on button press in radiobuttonObsEchos.
function radiobuttonObsEchos_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonObsEchos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.ObsEchosFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobuttonObsEchos


% --- Executes on button press in pushbuttonPlotSonogram.
function pushbuttonPlotSonogram_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlotSonogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MySonogramPlot(handles.DataToAnalyze, handles.DataToPlot);
% MySonogramPlot(handles.DataToAnalyze.BatSonarEchosMat(2,:),...
%     handles.DataToAnalyze.BatSonarEchosMat(1,:), handles.DataToAnalyze.BatSonogram);
    


function editBatsNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editBatsNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ed1 = get(hObject,'String');
handles.SimParams.TotalBatsNumber = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 
% Hints: get(hObject,'String') returns contents of editBatsNumber as text
%        str2double(get(hObject,'String')) returns contents of editBatsNumber as a double


% --- Executes during object creation, after setting all properties.
function editBatsNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBatsNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPreysNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editPreysNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreysNumber as text
%        str2double(get(hObject,'String')) returns contents of editPreysNumber as a double
ed1 = get(hObject,'String');
handles.SimParams.TotalPreysNumber = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function editPreysNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreysNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxAllBatsPlot.
function checkboxAllBatsPlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAllBatsPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChBx1 = get(hObject,'Value');
handles.DataToPlot.PlotAllBatsFlag = ChBx1;
guidata(hObject,handles) 
% Hint: get(hObject,'Value') returns toggle state of checkboxAllBatsPlot


% --- Executes on selection change in listboxPLotBatNum.
function listboxPLotBatNum_Callback(hObject, eventdata, handles)
% hObject    handle to listboxPLotBatNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxPLotBatNum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxPLotBatNum


% --- Executes during object creation, after setting all properties.
function listboxPLotBatNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxPLotBatNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox19AllPreysPlot.
function checkbox19AllPreysPlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19AllPreysPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChBx1 = get(hObject,'Value');
handles.DataToPlot.PlotAllPreysFlag = ChBx1;
guidata(hObject,handles) 
% Hint: get(hObject,'Value') returns toggle state of checkbox19AllPreysPlot



function edit40PreyNumToPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit40PreyNumToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit40PreyNumToPlot as text
%        str2double(get(hObject,'String')) returns contents of edit40PreyNumToPlot as a double
ed1 = get(hObject,'String');
handles.DataToPlot.PreyNumToPlot = str2double(ed1);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit40PreyNumToPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit40PreyNumToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18BatNumToPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit18BatNumToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hied1 = get(hObject,'String');
ed1 = get(hObject,'String');
handles.DataToPlot.BatNumToPlot = str2double(ed1);
guidata(hObject,handles) % save the changes to the structurets: get(hObject,'String') returns contents of edit18BatNumToPlot as text
%        str2double(get(hObject,'String')) returns contents of edit18BatNumToPlot as a double


% --- Executes during object creation, after setting all properties.
function edit18BatNumToPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18BatNumToPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% % function edit41NoiseLevel_Callback(hObject, eventdata, handles)
% % % hObject    handle to edit41NoiseLevel (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    structure with handles and user data (see GUIDATA)
% % 
% % % Hints: get(hObject,'String') returns contents of edit41NoiseLevel as text
% % %        str2double(get(hObject,'String')) returns contents of edit41NoiseLevel as a double
% % ed1 = get(hObject,'String');
% % handles.BatSonarParams.NoiseLeveldB = str2double(ed1);
% % guidata(hObject,handles)
% % 
% % % --- Executes during object creation, after setting all properties.
% % function edit41NoiseLevel_CreateFcn(hObject, eventdata, handles)
% % % hObject    handle to edit41NoiseLevel (see GCBO)
% % % eventdata  reserved - to be defined in a future version of MATLAB
% % % handles    empty - handles not created until after all CreateFcns called
% % 
% % % Hint: edit controls usually have a white background on Windows.
% % %       See ISPC and COMPUTER.
% % if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
% %     set(hObject,'BackgroundColor','white');
% % end



function edit42TargetWingLength_Callback(hObject, eventdata, handles)
% hObject    handle to edit42TargetWingLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42TargetWingLength as text
%        str2double(get(hObject,'String')) returns contents of edit42TargetWingLength as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.TargetWingLength = str2double(ed1);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit42TargetWingLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42TargetWingLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton20_PlotInterferenceFlag.
function radiobutton20_PlotInterferenceFlag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton20_PlotInterferenceFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.InterferencePlotFlag =  get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton20_PlotInterferenceFlag


% --- Executes on button press in checkbox_MissDetectionProbFlag.
function checkbox_MissDetectionProbFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MissDetectionProbFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChBx1 = get(hObject,'Value');
handles.BatSonarParams.MissDetectionProbFlag = ChBx1;
guidata(hObject,handles) % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of checkbox_MissDetectionProbFlag



function edit43_ProbToDetect_Callback(hObject, eventdata, handles)
% hObject    handle to edit43_ProbToDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43_ProbToDetect as text
%        str2double(get(hObject,'String')) returns contents of edit43_ProbToDetect as a double

ed1 = get(hObject,'String');
handles.BatSonarParams.ProbToDetect = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit43_ProbToDetect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43_ProbToDetect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17_ChooseTest.
% hObject    handle to pushbutton17_ChooseTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % % [FileName,PathName,FilterIndex] = uigetfile('C:\Users\DELL\Documents\òåîø\ìéîåãéí\îç÷ø\MATLAB Bat Simulation\Experiment Scripts\*.m');
% % % PathFileName = [PathName,FileName];
% % % % load(PathFileName);
% % % % LoadedExperimentScript = load(PathFileName);
% % % handles.Experiment.ScriptTOrun = PathFileName;
% % % set(handles.text60_ScriptName,'String', FileName); % text feedback
% % % guidata(hObject,handles) % save the changes to the structure 

% --- Executes on button press in pushbutton18_ChooseSaveDir.
function pushbutton18_ChooseSaveDir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18_ChooseSaveDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20_RunExp.
function pushbutton20_RunExp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20_RunExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathFileName = handles.Experiment.ScriptTOrun ;
run(PathFileName)
function edit_NumOfTrials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NumOfTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NumOfTrials as text
%        str2double(get(hObject,'String')) returns contents of edit_NumOfTrials as a double


% --- Executes during object creation, after setting all properties.
function edit_NumOfTrials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NumOfTrials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function textSummaryGroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textSummaryGroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkboxSaveAnimation.
function checkboxSaveAnimation_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSaveAnimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChBx1 = get(hObject,'Value');
handles.SaveAminationFlag = ChBx1;
guidata(hObject,handles); % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of checkboxSaveAnimation


% --- Executes on button press in pushbuttonMyBatUmweltPlot.
function pushbuttonMyBatUmweltPlot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMyBatUmweltPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.DataToAnalyze, handles.DataToPlot
MyBatUmweltOnLine( handles.DataToAnalyze.BAT(handles.DataToPlot.BatNumToPlot) ,...
    handles.DataToAnalyze.AllParams)



function edit45DetectionTHdB_Callback(hObject, eventdata, handles)
% hObject    handle to edit45DetectionTHdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ed1 = get(hObject,'String');
handles.BatSonarParams.PulseDetectionTH = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 



% --- Executes during object creation, after setting all properties.
function edit45DetectionTHdB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45DetectionTHdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuReceiverType.
function popupmenuReceiverType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuReceiverType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuReceiverType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuReceiverType
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
% handles.BatSonarParams.ReceiverType  = popupValue;
%%%% Somthing Wrong with CorrealtionDetector .. so that's the solution
switch popupValue
    case 'MaxPowerTH'
        handles.BatSonarParams.ReceiverType = 'MaxPowerTH';
    
    case 'FilterBank'
        handles.BatSonarParams.ReceiverType = 'FilterBank';

    case 'Correlation'
        handles.BatSonarParams.ReceiverType = 'CorrelationDetector';
    
    case 'TimeSequenceTH'
        handles.BatSonarParams.ReceiverType = 'TimeSequenceTH';
        
    case 'NoMasking'
        handles.BatSonarParams.ReceiverType = 'NoMasking';
 
end % switch popupValue
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function popupmenuReceiverType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuReceiverType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% set(hObject,'String',{'MaxPowerTH' ; 'Correlation' ; 'Binaural'});



% function edit46_MemorySize_Callback(hObject, eventdata, handles)
% % hObject    handle to edit46_MemorySize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit46_MemorySize as text
% %        str2double(get(hObject,'String')) returns contents of edit46_MemorySize as a double
% ed1 = get(hObject,'String');
% handles.BatFlightParams.HuntedMemorySize = str2double(ed1);
% guidata(hObject,handles) % save the changes to the structure 
% 
% % --- Executes during object creation, after setting all properties.
% function edit46_MemorySize_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit46_MemorySize (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



function edit47_MemoryKHunted_Callback(hObject, eventdata, handles)
% hObject    handle to edit47_MemoryKHunted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47_MemoryKHunted as text
%        str2double(get(hObject,'String')) returns contents of edit47_MemoryKHunted as a double
ed1 = get(hObject,'String');
handles.BatFlightParams.NumOfSamePreyForAppraoach = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit47_MemoryKHunted_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47_MemoryKHunted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit48_SIRneededTH_Callback(hObject, eventdata, handles)
% hObject    handle to edit48_SIRneededTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48_SIRneededTH as text
%        str2double(get(hObject,'String')) returns contents of edit48_SIRneededTH as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.Signal2InterferenceRatio = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit48_SIRneededTH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48_SIRneededTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox22_FwdMaskingFlag.
function checkbox22_FwdMaskingFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22_FwdMaskingFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22_FwdMaskingFlag
ChBx1 = get(hObject,'Value');
handles.BatSonarParams.MaskingFwdFlag = ChBx1;
guidata(hObject,handles); % save the changes to the structure 


function edit49_FwdMaskingTH_Callback(hObject, eventdata, handles)
% hObject    handle to edit49_FwdMaskingTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49_FwdMaskingTH as text
%        str2double(get(hObject,'String')) returns contents of edit49_FwdMaskingTH as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.MaskingFwdSIRth = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function edit49_FwdMaskingTH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49_FwdMaskingTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_FwdMakingTime_Callback(hObject, eventdata, handles)
% hObject    handle to edit50_FwdMakingTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50_FwdMakingTime as text
%        str2double(get(hObject,'String')) returns contents of edit50_FwdMakingTime as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.MaskingFwdTime = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function edit50_FwdMakingTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50_FwdMakingTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit51_MinPulsePower_Callback(hObject, eventdata, handles)
% hObject    handle to edit51_MinPulsePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit51_MinPulsePower as text
%        str2double(get(hObject,'String')) returns contents of edit51_MinPulsePower as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.PulseMinPower = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function edit51_MinPulsePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit51_MinPulsePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit52_Dist2Buzz_Callback(hObject, eventdata, handles)
% hObject    handle to edit52_Dist2Buzz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52_Dist2Buzz as text
%        str2double(get(hObject,'String')) returns contents of edit52_Dist2Buzz as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.DistanceToBuzz = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 



% --- Executes during object creation, after setting all properties.
function edit52_Dist2Buzz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52_Dist2Buzz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit53_Dist2Appraoch_Callback(hObject, eventdata, handles)
% hObject    handle to edit53_Dist2Appraoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit53_Dist2Appraoch as text
%        str2double(get(hObject,'String')) returns contents of edit53_Dist2Appraoch as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.DistanceFromPreyToApproach = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 



% --- Executes during object creation, after setting all properties.
function edit53_Dist2Appraoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53_Dist2Appraoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textSummaryGroup2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textSummaryGroup2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit54_BatMinVelocity_Callback(hObject, eventdata, handles)
% hObject    handle to edit54_BatMinVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54_BatMinVelocity as text
%        str2double(get(hObject,'String')) returns contents of edit54_BatMinVelocity as a double
ed1 = get(hObject,'String');
handles.BatFlightParams.MinVelocity = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 


% --- Executes during object creation, after setting all properties.
function edit54_BatMinVelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54_BatMinVelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% 
% % --- Executes on button press in radiobutton21.
% function radiobutton21_Callback(hObject, eventdata, handles)
% % hObject    handle to radiobutton21 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton21




% --- Executes on button press in checkbox29_LocalizationErrFlag.
function checkbox29_LocalizationErrFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29_LocalizationErrFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ChBx1 = get(hObject,'Value');
    handles.BatSonarParams.LocalizationErrFlag = ChBx1;
    guidata(hObject,handles) % save the changes to the structure 
    
% Hint: get(hObject,'Value') returns toggle state of checkbox29_LocalizationErrFlag



function edit_JARDelta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_JARDelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ed1 = get(hObject,'String');
handles.BatSonarParams.JARDelta = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 
% Hints: get(hObject,'String') returns contents of edit_JARDelta as text
%        str2double(get(hObject,'String')) returns contents of edit_JARDelta as a double


% --- Executes during object creation, after setting all properties.
function edit_JARDelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_JARDelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_JARmem_Callback(hObject, eventdata, handles)
% hObject    handle to edit_JARmem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ed1 = get(hObject,'String');
handles.BatSonarParams.JARMem = str2double(ed1);
guidata(hObject,handles) % save the changes to the structure 
% Hints: get(hObject,'String') returns contents of edit_JARmem as text
%        str2double(get(hObject,'String')) returns contents of edit_JARmem as a double


% --- Executes during object creation, after setting all properties.
function edit_JARmem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_JARmem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttonJAR_OnOff.
function radiobuttonJAR_OnOff_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonJAR_OnOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatSonarParams.JARMode = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobuttonJAR_OnOff


% --- Executes on key press with focus on radiobuttonSonarTimes and none of its controls.
function radiobuttonSonarTimes_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to radiobuttonSonarTimes (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on radiobuttonJAR_OnOff and none of its controls.
function radiobuttonJAR_OnOff_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to radiobuttonJAR_OnOff (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxJARMode.
handles.BatSonarParams.JARMode = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 

% hObject    handle to checkboxJARMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxJARMode


% --- Executes on button press in checkbox_JARMode.
function checkbox_JARMode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_JARMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatSonarParams.JARMode = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkbox_JARMode


% % --- Executes on button press in radiobutton23_BiStatMode.
% function radiobutton23_BiStatMode_Callback(hObject, eventdata, handles)
% % hObject    handle to radiobutton23_BiStatMode (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% handles.BatSonarParams.BiStatMode = get(hObject,'Value');	
% guidata(hObject,handles); % save the changes to the structure 
% % Hint: get(hObject,'Value') returns toggle state of radiobutton23_BiStatMode


% --- Executes on button press in checkbox32_BiStatMode.
function checkbox32_BiStatMode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32_BiStatMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatSonarParams.BiStatMode = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of checkbox32_BiStatMode



function edit58_PulseDetectionTH_Callback(hObject, eventdata, handles)
% hObject    handle to edit58_PulseDetectionTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit58_PulseDetectionTH as text
%        str2double(get(hObject,'String')) returns contents of edit58_PulseDetectionTH as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.PulseDetectionTH = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 



% --- Executes during object creation, after setting all properties.
function edit58_PulseDetectionTH_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit58_PulseDetectionTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox33_ReactToBatsFlag.
function checkbox33_ReactToBatsFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox33_ReactToBatsFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatFlightParams.ReactToBatsFlag = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkbox33_ReactToBatsFlag


% --- Executes on button press in checkbox34_PhantomEchoesFromConsFlag.
function checkbox34_PhantomEchoesFromConsFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34_PhantomEchoesFromConsFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.BatSonarParams.PhantomEchoesFromConsFlag = get(hObject,'Value');	
guidata(hObject,handles); % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkbox34_PhantomEchoesFromConsFlag



function edit59_NoiseLeveldB_Callback(hObject, eventdata, handles)
% hObject    handle to edit59_NoiseLeveldB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit59_NoiseLeveldB as text
%        str2double(get(hObject,'String')) returns contents of edit59_NoiseLeveldB as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.NoiseLeveldB = str2double(ed1);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit59_NoiseLeveldB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit59_NoiseLeveldB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Hint: get(hObject,'Value') returns toggle state of radiobutton24_PhantomCheckSim


% --- Executes on button press in checkbox35_PhantomCheckSim.
function checkbox35_PhantomCheckSim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox35_PhantomCheckSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox35_PhantomCheckSim
handles.BatSonarParams.PhantomEchoesCheckSimiliarity = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in radiobutton25_Directivity_Mode.
function radiobutton25_Directivity_Mode_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton25_Directivity_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton25_Directivity_Mode
handles.BatSonarParams.DirectivityMode = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 




function edit60_ClutterGainDB_Callback(hObject, eventdata, handles)
% hObject    handle to edit60_ClutterGainDB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit60_ClutterGainDB as text
%        str2double(get(hObject,'String')) returns contents of edit60_ClutterGainDB as a double
ed1 = get(hObject,'String');
handles.BatSonarParams.ClutterGainDB = str2double(ed1);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit60_ClutterGainDB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit60_ClutterGainDB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_TestMode.
function popupmenu_TestMode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_TestMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_TestMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_TestMode
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
handles.SimParams.TestMode = popupValue;
guidata(hObject,handles);% save the struct
 

% --- Executes during object creation, after setting all properties.
function popupmenu_TestMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_TestMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_swarm_start_pos_std_Callback(hObject, eventdata, handles)
% hObject    handle to edit_swarm_start_pos_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_swarm_start_pos_std as text
%        str2double(get(hObject,'String')) returns contents of edit_swarm_start_pos_std as a double
ed1 = get(hObject,'String');
handles.SimParams.swarm_initial_ybat_std = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit_swarm_start_pos_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_swarm_start_pos_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function TimePlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: place code in OpeningFcn to populate TimePlot


% --- Executes on selection change in popupmenu10_DetetctBehavior.
function popupmenu10_DetetctBehavior_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10_DetetctBehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
switch popupValue
    case 'Prey'
        handles.SimParams.DetectPrey = 1;
        handles.SimParams.DetectConsps = 0;
        handles.SimParams.DetectObs = 0; 
    case 'Conspecifics'
        handles.SimParams.DetectPrey = 0;
        handles.SimParams.DetectConsps = 1;
        handles.SimParams.DetectObs = 0; 
    case 'Obstacles'
        handles.SimParams.DetectPrey = 0;
        handles.SimParams.DetectConsps = 0;
        handles.SimParams.DetectObs = 1; 
        
end % switch

guidata(hObject,handles);% save the struct

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10_DetetctBehavior contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10_DetetctBehavior


% --- Executes during object creation, after setting all properties.
function popupmenu10_DetetctBehavior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10_DetetctBehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox10.
function listbox10_Callback(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox10


% --- Executes during object creation, after setting all properties.
function listbox10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton26_ConspJamming.
function radiobutton26_ConspJamming_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton26_ConspJamming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.MaskingByConsps = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobutton26_ConspJamming


% --- Executes on button press in radiobutton27_ClutterMode.
function radiobutton27_ClutterMode_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton27_ClutterMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.MaskingByClutter = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of radiobutton27_ClutterMode


% --- Executes on button press in radiobutton28_Acoustics.
function radiobutton28_Acoustics_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton28_Acoustics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.AcousticsCalcFlag = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of radiobutton28_Acoustics


% --- Executes on button press in checkbox37_MaskingByConsps.
function checkbox37_MaskingByConsps_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox37_MaskingByConsps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.MaskingByConsps = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of checkbox37_MaskingByConsps


% --- Executes on button press in checkbox38_MaskingByClutter.
function checkbox38_MaskingByClutter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox38_MaskingByClutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.MaskingByClutter = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkbox38_MaskingByClutter


% --- Executes on button press in checkbox_JamMothFlg.
function checkbox_JamMothFlg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_JamMothFlg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PreyFlightParams.React2BatsFlag = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkbox_JamMothFlg


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PreyFlightParams.jamPreySignal = generateAcousicPreyJamming(handles.PreyFlightParams, handles.SimParams.FsAcoustic);
handles.PreyFlightParams.JamLoadFlg  = 1;
set(handles.checkbox_LoadJamSig,'Value', handles.PreyFlightParams.JamLoadFlg)
guidata(hObject,handles)

% --- Executes on button press in radiobutton_loadJamSig.
function radiobutton_loadJamSig_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_loadJamSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_loadJamSig


% --- Executes on button press in checkbox_LoadJamSig.
function checkbox_LoadJamSig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LoadJamSig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PreyFlightParams.JamLoadFlg = get(hObject,'Value');
guidata(hObject,handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_LoadJamSig


% --- Executes on button press in radiobutton_PreyJam.
function radiobutton_PreyJam_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_PreyJam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataToPlot.JamPreyPosFlag =  get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure  
% Hint: get(hObject,'Value') returns toggle state of radiobutton_PreyJam



function edit_jamMothTxLvl_Callback(hObject, eventdata, handles)
% hObject    handle to edit_jamMothTxLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_jamMothTxLvl as text
%        str2double(get(hObject,'String')) returns contents of edit_jamMothTxLvl as a double
ed1 = get(hObject,'String');
handles.PreyFlightParams.JamTxLvl = str2double(ed1);
guidata(hObject,handles); % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit_jamMothTxLvl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_jamMothTxLvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_PlotAcoustics.
function pushbutton_PlotAcoustics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlotAcoustics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
if handles.DataToPlot.PlotAllBatsFlag
    BatNum = 'all';
else % if handles.DataToPlot.PlotAllBatsFlag
    BatNum = handles.DataToPlot.BatNumToPlot;
end % if handles.DataToPlot.PlotAllBatsFlag

if handles.checkbox41_NewFig.Value == 0
    axHandle = handles.TimePlot;
else
    figure();
    axHandle = gca;
end % if handles.checkbox41_NewFig.Value == 1

My_AcousticSig_VecTest (handles.DataToAnalyze, BatNum, axHandle);

if handles.checkbox42_LegendTime.Value == 1
   axHandle.Legend.Visible = 'on';
else
   axHandle.Legend.Visible = 'off';    
end % if handles.checkbox41_NewFig.Value == 1

guidata(hObject,handles); % save the changes to the structure 


% --- Executes on button press in checkbox41_NewFig.
function checkbox41_NewFig_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox41_NewFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox41_NewFig

handles.checkbox41_NewFig.Value =  get(hObject,'Value');
guidata(hObject,handles); % save the changes to the structure 


% --- Executes on button press in checkbox42_LegendTime.
function checkbox42_LegendTime_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox42_LegendTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox42_LegendTime
if get(hObject,'Value') == 1
    handles.TimePlot.Legend.Visible = 'on';
else
    handles.TimePlot.Legend.Visible = 'off';    
end 



function edit64_pulseToAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to edit64_pulseToAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit64_pulseToAnalyze as text
%        str2double(get(hObject,'String')) returns contents of edit64_pulseToAnalyze as a double


% --- Executes during object creation, after setting all properties.
function edit64_pulseToAnalyze_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit64_pulseToAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton26_analyzePulse.
function pushbutton26_analyzePulse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26_analyzePulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PulseNum = str2double(get(handles.edit64_pulseToAnalyze, 'String'));
BatNum = handles.DataToPlot.BatNumToPlot;
Analayze_Sig(handles.DataToAnalyze, PulseNum, BatNum)


% --- Executes on button press in checkbox43_ClassifyFlag.
function checkbox43_ClassifyFlag_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox43_ClassifyFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox43_ClassifyFlag
handles.BatSonarParams.PreyClassifierFlag = get(hObject,'Value');

guidata(hObject,handles) % save the changes to the structure 


% --- Executes on selection change in popupmenu_receiverType.
function popupmenu_receiverType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_receiverType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_receiverType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_receiverType

PopupContents = cellstr(get(hObject,'String'));
popupValue= PopupContents{get(hObject,'Value')};
% handles.BatSonarParams.ReceiverType  = popupValue;
%%%% Somthing Wrong with CorrealtionDetector .. so that's the solution
switch popupValue
    case 'MaxPowerTH'
        handles.BatSonarParams.ReceiverType = 'MaxPowerTH';
    
    case 'FilterBank'
        handles.BatSonarParams.ReceiverType = 'FilterBank';

    case 'CorrelationDetector'
        handles.BatSonarParams.ReceiverType = 'CorrelationDetector';
    
    case 'TimeSequenceTH'
        handles.BatSonarParams.ReceiverType = 'TimeSequenceTH';
        
    case 'NoMasking'
        handles.BatSonarParams.ReceiverType = 'NoMasking';
 
end % switch popupValue
guidata(hObject,handles) % save the changes to the structure 

% popupmenu_receiverType
% popupmenuReceiverType

% --- Executes during object creation, after setting all properties.
function popupmenu_receiverType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_receiverType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AnimationTimeRate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AnimationTimeRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AnimationTimeRate as text
%        str2double(get(hObject,'String')) returns contents of edit_AnimationTimeRate as a double

handles.AnimationTimeRate = str2num(get(hObject,'String')) ;

guidata(hObject,handles) % save the changes to the structure 

% --- Executes during object creation, after setting all properties.
function edit_AnimationTimeRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AnimationTimeRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxDetectPrey.
function checkboxDetectPrey_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDetectPrey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.DetectPrey = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of checkboxDetectPrey


% --- Executes on button press in checkboxDetectClutter.
function checkboxDetectClutter_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDetectClutter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.DetectObs = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkboxDetectClutter


% --- Executes on button press in checkboxDetectConsps.
function checkboxDetectConsps_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDetectConsps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SimParams.DetectConsps = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 
% Hint: get(hObject,'Value') returns toggle state of checkboxDetectConsps


% --- Executes on button press in radiobuttonAnumatOnlyCurrent.
function radiobuttonAnumatOnlyCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonAnumatOnlyCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.AnimmateCurrentOnly = get(hObject,'Value');
guidata(hObject,handles) % save the changes to the structure 

% Hint: get(hObject,'Value') returns toggle state of radiobuttonAnumatOnlyCurrent



function edit66_Callback(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit66 as text
%        str2double(get(hObject,'String')) returns contents of edit66 as a double
handles.DataToPlot.AnimateStartTime = str2double(get(hObject,'String'));
guidata(hObject,handles) % save the changes to the stru

% --- Executes during object creation, after setting all properties.
function edit66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double
handles.DataToPlot.AnimateEndTime = str2double(get(hObject,'String'));
guidata(hObject,handles) % save the changes to the stru

% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonAnalyzeMan.
function pushbuttonAnalyzeMan_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAnalyzeMan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PulseNum = str2double(get(handles.edit64_pulseToAnalyze, 'String'));
BatNum = handles.DataToPlot.BatNumToPlot;
if handles.checkbox41_NewFig.Value == 0
    axHandle = handles.xyPlot;
else
    figure(); hold on;
    axHandle = gca;
end % if handles.checkbox41_NewFig.Value == 1
myManuverAnlysisPlot(handles.DataToAnalyze, BatNum, PulseNum, axHandle);
