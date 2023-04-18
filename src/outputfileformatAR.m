function varargout = outputfileformatAR(varargin)
% OUTPUTFILEFORMATAR MATLAB code for outputfileformatAR.fig
%      OUTPUTFILEFORMATAR, by itself, creates a new OUTPUTFILEFORMATAR or raises the existing
%      singleton*.
%
%      H = OUTPUTFILEFORMATAR returns the handle to a new OUTPUTFILEFORMATAR or the handle to
%      the existing singleton*.
%
%      OUTPUTFILEFORMATAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OUTPUTFILEFORMATAR.M with the given input arguments.
%
%      OUTPUTFILEFORMATAR('Property','Value',...) creates a new OUTPUTFILEFORMATAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before outputfileformatAR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to outputfileformatAR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help outputfileformatAR

% Last Modified by GUIDE v2.5 28-Jan-2015 11:35:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @outputfileformatAR_OpeningFcn, ...
                   'gui_OutputFcn',  @outputfileformatAR_OutputFcn, ...
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


% --- Executes just before outputfileformatAR is made visible.
function outputfileformatAR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to outputfileformatAR (see VARARGIN)

% Choose default command line output for outputfileformatAR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes outputfileformatAR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = outputfileformatAR_OutputFcn(hObject, eventdata, handles) 
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
