% Copyright (c) 2014, Swallowing
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of AHI nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The data is based on LS-120 monitoring with the recorded id: 17812 
% The time frame is between 0:20 - 0:50 (30 min frame)
% Last Modified (record list below):
% ---- 11-Nov-2014 by Adam Lin
% ---- record required if edit by others
%

function varargout = Enge2(varargin)
% ENGE2 M-file for Enge2.fig
%      ENGE2, by itself, creates a new ENGE2 or raises the existing
%      singleton*.
%
%      H = ENGE2 returns the handle to a new ENGE2 or the handle to
%      the existing singleton*.
%
%      ENGE2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENGE2.M with the given input arguments.
%
%      ENGE2('Property','Value',...) creates a new ENGE2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Enge2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Enge2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Enge2

% Last Modified by GUIDE v2.5 17-Nov-2014 10:53:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Enge2_OpeningFcn, ...
                   'gui_OutputFcn',  @Enge2_OutputFcn, ...
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


% --- Executes just before Enge2 is made visible.
function Enge2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Enge2 (see VARARGIN)

% Choose default command line output for Enge2
handles.data = 'no_path';
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Enge2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Enge2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function max_value_Callback(hObject, eventdata, handles)
% hObject    handle to max_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
max_val = get(hObject, 'value');
max_val = max_val/10;
handles.max_value = max_val;
set(handles.txt_max, 'string', handles.max_value);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function min_value_Callback(hObject, eventdata, handles)
% hObject    handle to min_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
min_val = get(hObject, 'value');
min_val = min_val/10;
handles.min_value = min_val;
set(handles.txt_min, 'string', handles.min_value);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in diagnosis.
function diagnosis_Callback(hObject, eventdata, handles)
% hObject    handle to diagnosis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%fprintf('Data max: %.2f Data min %.2f \n', handles.max_value, handles.min_value);
%showPlot(data, handles.max_value, handles.min_value);
if strcmpi(handles.data, 'no_path')
    return;
else
    set(handles.plot_val_time,'Value',1);
    axes(handles.show_plot);
    output = showPlot(handles.data, handles.max_value, handles.min_value, 'diagnosis', 1);
    
    axes(handles.show_org_plot);
    showPlot(handles.data, 1, 1, 'original', 1);

    % display AHI message
    set(handles.ahi, 'string', sprintfc('%d T/H', output(1)));
    set(handles.max_ahi, 'string', sprintfc('%0.1f Sec', output(2)));
end
%axes(handles.org_plot);
%showPlot(handles.data, handles.max_value, handles.min_value, 'cancellation', 1);

% --- Executes during object creation, after setting all properties.
function txt_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function txt_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in plot_val_time.
function plot_val_time_Callback(hObject, eventdata, handles)
% hObject    handle to plot_val_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_val_time contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_val_time
times = get(hObject, 'value');
axes(handles.org_plot);
showPlot(handles.data, 1, 1, 'original', times);


% --- Executes during object creation, after setting all properties.
function plot_val_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_val_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browerFileName.
function browerFileName_Callback(hObject, eventdata, handles)
% hObject    handle to browerFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*csv'},'File Selector');
if pathname ~= 0
    fullpathname = strcat(pathname, filename);
    set(handles.txt_pathname, 'string', fullpathname);
    handles.data = fullpathname;
end
guidata(hObject, handles);


% --- Executes on button press in btn_import.
function btn_import_Callback(hObject, eventdata, handles)
% hObject    handle to btn_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(handles.data, 'no_path')
    return;
else
    axes(handles.org_plot);
    showPlot(handles.data, 0.01, 0.01, 'original', 1);
end


% --- Executes on selection change in plot_val_time2.
function plot_val_time2_Callback(hObject, eventdata, handles)
% hObject    handle to plot_val_time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_val_time2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_val_time2
times = get(hObject, 'value');
axes(handles.show_org_plot);
showPlot(handles.data, 1, 1, 'original', times);

% --- Executes during object creation, after setting all properties.
function plot_val_time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_val_time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check.
function check_Callback(hObject, eventdata, handles)
% hObject    handle to check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check
checked = get(hObject, 'Value');

if checked == 1
    axes(handles.low_show_high_plot);
    showPlot(handles.data, 1, 1, 'cancellation', 1);
else
    return;
end