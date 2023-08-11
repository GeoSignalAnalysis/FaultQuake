function varargout = WeightFileFormat(varargin)
%WEIGHTFILEFORMAT M-file for WeightFileFormat.fig
%      WEIGHTFILEFORMAT, by itself, creates a new WEIGHTFILEFORMAT or raises the existing
%      singleton*.
%
%      H = WEIGHTFILEFORMAT returns the handle to a new WEIGHTFILEFORMAT or the handle to
%      the existing singleton*.
%
%      WEIGHTFILEFORMAT('Property','Value',...) creates a new WEIGHTFILEFORMAT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to WeightFileFormat_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      WEIGHTFILEFORMAT('CALLBACK') and WEIGHTFILEFORMAT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in WEIGHTFILEFORMAT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WeightFileFormat

% Last Modified by GUIDE v2.5 23-Sep-2015 15:25:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @WeightFileFormat_OpeningFcn, ...
                   'gui_OutputFcn',  @WeightFileFormat_OutputFcn, ...
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


% --- Executes just before WeightFileFormat is made visible.
function WeightFileFormat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for WeightFileFormat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes WeightFileFormat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = WeightFileFormat_OutputFcn(hObject, eventdata, handles)
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
