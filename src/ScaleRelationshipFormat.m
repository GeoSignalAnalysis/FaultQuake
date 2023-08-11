function varargout = ScaleRelationshipFormat(varargin)
%SCALERELATIONSHIPFORMAT M-file for ScaleRelationshipFormat.fig
%      SCALERELATIONSHIPFORMAT, by itself, creates a new SCALERELATIONSHIPFORMAT or raises the existing
%      singleton*.
%
%      H = SCALERELATIONSHIPFORMAT returns the handle to a new SCALERELATIONSHIPFORMAT or the handle to
%      the existing singleton*.
%
%      SCALERELATIONSHIPFORMAT('Property','Value',...) creates a new SCALERELATIONSHIPFORMAT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to ScaleRelationshipFormat_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SCALERELATIONSHIPFORMAT('CALLBACK') and SCALERELATIONSHIPFORMAT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SCALERELATIONSHIPFORMAT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ScaleRelationshipFormat

% Last Modified by GUIDE v2.5 03-Feb-2015 16:46:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ScaleRelationshipFormat_OpeningFcn, ...
                   'gui_OutputFcn',  @ScaleRelationshipFormat_OutputFcn, ...
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


% --- Executes just before ScaleRelationshipFormat is made visible.
function ScaleRelationshipFormat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for ScaleRelationshipFormat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ScaleRelationshipFormat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ScaleRelationshipFormat_OutputFcn(hObject, eventdata, handles)
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
