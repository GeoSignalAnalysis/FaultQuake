function varargout = outputfileformatMB(varargin)
% OUTPUTFILEFORMATMB MATLAB code for outputfileformatMB.fig
%      OUTPUTFILEFORMATMB, by itself, creates a new OUTPUTFILEFORMATMB or raises the existing
%      singleton*.
%
%      H = OUTPUTFILEFORMATMB returns the handle to a new OUTPUTFILEFORMATMB or the handle to
%      the existing singleton*.
%
%      OUTPUTFILEFORMATMB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OUTPUTFILEFORMATMB.M with the given input arguments.
%
%      OUTPUTFILEFORMATMB('Property','Value',...) creates a new OUTPUTFILEFORMATMB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before outputfileformatMB_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to outputfileformatMB_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help outputfileformatMB

% Last Modified by GUIDE v2.5 25-Jun-2015 12:52:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @outputfileformatMB_OpeningFcn, ...
                   'gui_OutputFcn',  @outputfileformatMB_OutputFcn, ...
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


% --- Executes just before outputfileformatMB is made visible.
function outputfileformatMB_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to outputfileformatMB (see VARARGIN)

% Choose default command line output for outputfileformatMB
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes outputfileformatMB wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = outputfileformatMB_OutputFcn(hObject, eventdata, handles) 
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
