function varargout = infoRP(varargin)
% INFORP MATLAB code for infoRP.fig
%      INFORP, by itself, creates a new INFORP or raises the existing
%      singleton*.
%
%      H = INFORP returns the handle to a new INFORP or the handle to
%      the existing singleton*.
%
%      INFORP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INFORP.M with the given input arguments.
%
%      INFORP('Property','Value',...) creates a new INFORP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before infoRP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to infoRP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help infoRP

% Last Modified by GUIDE v2.5 25-Jun-2015 12:51:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @infoRP_OpeningFcn, ...
                   'gui_OutputFcn',  @infoRP_OutputFcn, ...
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


% --- Executes just before infoRP is made visible.
function infoRP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to infoRP (see VARARGIN)

% Choose default command line output for infoRP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes infoRP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = infoRP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
