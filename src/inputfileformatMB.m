function varargout = inputfileformatMB(varargin)
% INPUTFILEFORMATMB MATLAB code for inputfileformatMB.fig
%      INPUTFILEFORMATMB, by itself, creates a new INPUTFILEFORMATMB or raises the existing
%      singleton*.
%
%      H = INPUTFILEFORMATMB returns the handle to a new INPUTFILEFORMATMB or the handle to
%      the existing singleton*.
%
%      INPUTFILEFORMATMB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTFILEFORMATMB.M with the given input arguments.
%
%      INPUTFILEFORMATMB('Property','Value',...) creates a new INPUTFILEFORMATMB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputfileformatMB_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputfileformatMB_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputfileformatMB

% Last Modified by GUIDE v2.5 25-Jun-2015 12:52:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputfileformatMB_OpeningFcn, ...
                   'gui_OutputFcn',  @inputfileformatMB_OutputFcn, ...
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


% --- Executes just before inputfileformatMB is made visible.
function inputfileformatMB_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inputfileformatMB (see VARARGIN)

% Choose default command line output for inputfileformatMB
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inputfileformatMB wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputfileformatMB_OutputFcn(hObject, eventdata, handles) 
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
