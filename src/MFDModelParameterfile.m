function varargout = MFDModelParameterfile(varargin)
%MFDMODELPARAMETERFILE M-file for MFDModelParameterfile.fig
%      MFDMODELPARAMETERFILE, by itself, creates a new MFDMODELPARAMETERFILE or raises the existing
%      singleton*.
%
%      H = MFDMODELPARAMETERFILE returns the handle to a new MFDMODELPARAMETERFILE or the handle to
%      the existing singleton*.
%
%      MFDMODELPARAMETERFILE('Property','Value',...) creates a new MFDMODELPARAMETERFILE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to MFDModelParameterfile_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      MFDMODELPARAMETERFILE('CALLBACK') and MFDMODELPARAMETERFILE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in MFDMODELPARAMETERFILE.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MFDModelParameterfile

% Last Modified by GUIDE v2.5 07-May-2023 15:40:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MFDModelParameterfile_OpeningFcn, ...
                   'gui_OutputFcn',  @MFDModelParameterfile_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before MFDModelParameterfile is made visible.
function MFDModelParameterfile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for MFDModelParameterfile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MFDModelParameterfile wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MFDModelParameterfile_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
