function varargout = TDSviewer_m4(varargin)
% TDSVIEWER_M4 MATLAB code for TDSviewer_m4.fig
%      TDSVIEWER_M4, by itself, creates a new TDSVIEWER_M4 or raises the existing
%      singleton*.
%
%      H = TDSVIEWER_M4 returns the handle to a new TDSVIEWER_M4 or the handle to
%      the existing singleton*.
%
%      TDSVIEWER_M4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TDSVIEWER_M4.M with the given input arguments.
%
%      TDSVIEWER_M4('Property','Value',...) creates a new TDSVIEWER_M4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TDSviewer_m4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TDSviewer_m4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TDSviewer_m4

% Last Modified by GUIDE v2.5 04-May-2022 16:53:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TDSviewer_m4_OpeningFcn, ...
                   'gui_OutputFcn',  @TDSviewer_m4_OutputFcn, ...
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


% --- Executes just before TDSviewer_m4 is made visible.
function TDSviewer_m4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TDSviewer_m4 (see VARARGIN)

% Choose default command line output for TDSviewer_m4
handles.output = hObject;

set(gcf,'toolbar','figure')

set(handles.uibuttongroup1, 'Visible', 'on')

box(handles.axes1, 'on');
box(handles.axes2, 'on');

% --- Initialize sliders of filters

set(handles.slider2, 'Value', 5);
set(handles.slider2, 'Min', 5);
set(handles.slider2, 'Max', 100);

set(handles.slider4, 'Value', 1);
set(handles.slider4, 'Min', 1);
set(handles.slider4, 'Max', 20);

set(handles.slider1, 'Value', 0);
set(handles.slider1, 'Min', -5);
set(handles.slider1, 'Max', 5);

set(handles.slider3, 'Value', 0);
set(handles.slider3, 'Min', -2);
set(handles.slider3, 'Max', 2);

set(handles.slider5, 'Value', 1);
set(handles.slider5, 'Min', 0);
set(handles.slider5, 'Max', 5);

set(handles.slider6, 'Value', 5);
set(handles.slider6, 'Min', 5);
set(handles.slider6, 'Max', 100);


% --- Disable all sliders

set(handles.slider1, 'Enable', 'off')
set(handles.slider2, 'Enable', 'off')
set(handles.slider3, 'Enable', 'off')
set(handles.slider4, 'Enable', 'off')
set(handles.slider5, 'Enable', 'off')
set(handles.slider6, 'Enable', 'off')


% --- Load file list

selector = varargin{1};
path = varargin{2};

if selector
    selector_1 = 'STA';
    selector_2 = 'sta';
else
    selector_1 = 'STB';
    selector_2 = 'stb';
end

path = [path, '\', selector_1, '_TDS'];

files = dir(path);
files = files(3:end);
names_tds = cell(size(files));
names_mag = cell(size(files));
for i = 1:length(names_tds)
    names_tds{i} = [path, '\', files(i).name, '\', ...
        selector_2, '_l2_tds_lang_', files(i).name(10:24), '_v01.cdf'];  
end

handles.names_tds = names_tds;
handles.names_mag = names_mag;

handles.eventindex = 1;

% --- Plot first event
data = spdfcdfread(handles.names_tds{1});

[f, I, pE] = findmainfreq(data{6}, mean(diff(data{2})));
set(handles.slider2, 'Value', f(I));

plot(handles.axes1, data{2}, data{6});

plot(handles.axes2, f, pE);hold(handles.axes2, 'on');
plot(handles.axes2, (f(I)-0.5)*ones(1, 128), linspace(0, 1, 128), ...
    'LineWidth', 1, 'LineStyle', '-.', 'Color', 'black');
plot(handles.axes2, (f(I)+0.5)*ones(1, 128), linspace(0, 1, 128), ...
    'LineWidth', 1, 'LineStyle', '-.', 'Color', 'black');
plot(handles.axes2, ones(1, 128), linspace(0, 1, 128), ...
    'LineWidth', 1, 'LineStyle', '-.', 'Color', 'r');
plot(handles.axes2, 5*ones(1, 128), linspace(0, 1, 128), ...
    'LineWidth', 1, 'LineStyle', '-.', 'Color', 'g');
hold(handles.axes2, 'off');

set(handles.edit1, 'String', '1');
set(handles.edit2, 'String', num2str(length(names_tds)));
set(handles.edit2, 'Enable', 'off');

handles.data = data;
handles.mode = 'ori';
handles.comp = 'all';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TDSviewer_m4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function Ef = filterwaveform(E, dt, mode, varargin)
% Description:
%
% Apply Butterworth BPF on electric field time series
%
% Input arguments:
%
% E: time series of electric field
% t: time of observed electric field sequence
% fcen: central frequency
% bw: bandwidth
%
% Output arguments:
%
% Ef: bandpassed electric field
%

if dt < 0.005 % Determine sample rate
    Fs = 250;
else
    Fs = 125;
end

Rp = 1;  % Maximum attenuation in passband
Rs = 30; % Minimum attenuation in stopband

switch mode
    case 'bpf'
        fcen = varargin{1};
        bw = varargin{2};
        wp = [fcen - bw/2 fcen + bw/2]/(Fs/2);
        ws = [fcen - bw/2-1 fcen + bw/2+1]/(Fs/2);
        
        type = 'bandpass';
    case 'hpf'
        fp = varargin{1};
        fs = varargin{2};
        
        wp = fp/(Fs/2);
        ws = fs/(Fs/2);
        
        type = 'high';
    case 'lpf'
        fp = varargin{1};
        fs = varargin{2};
        
        wp = fp/(Fs/2);
        ws = fs/(Fs/2);
        
        type = 'low';
end

[N, Wn] = buttord(wp, ws, Rp, Rs);
[bb, ab] = butter(N, Wn, type);

Ef(:, 1) = filter(bb, ab, double(E(:, 1)));
Ef(:, 2) = filter(bb, ab, double(E(:, 2)));
Ef(:, 3) = filter(bb, ab, double(E(:, 3)));


function [f, I, P1] = findmainfreq(data, dt)
% Description:
%
% Apply Butterworth BPF on electric field time series
%
% Input arguments:
%
% E: time series of electric field
% t: time of observed electric field sequence
% fcen: central frequency
% bw: bandwidth
%
% Output arguments:
%
% Ef: bandpassed electric field
%

L = length(data);

Fs = 1/dt;

f = Fs*(0:(L/2))/L;

P2 = abs(fft(data)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

[~, I] = max(P1(f > 5));
I = I + length(P1(f <= 5));


function deleteplot(axes, color)
% Description:
%
% Apply Butterworth BPF on electric field time series
%
% Input arguments:
%
% E: time series of electric field
% t: time of observed electric field sequence
% fcen: central frequency
% bw: bandwidth
%
% Output arguments:
%
% Ef: bandpassed electric field
%

lines = get(axes, 'Children');

for i = 1:length(lines)
    if norm(get(lines(i), 'Color') - color) == 0
        delete(lines(i))
    end
end


function drawfields(data, axes, f1, f2, mode, comp)

    switch mode
        case 'ori'
            Ef = data{6};
        case 'bpf'
            deleteplot(axes(2), [0 0 0]);
            Ef = filterwaveform(data{6}, mean(diff(data{2})),...
                'bpf', f1, f2); % f1: central f2: bandwidth

            hold(axes(2), 'on');
            plot(axes(2), (f1 - f2/2)*ones(1, 128), ...
                linspace(0, 1, 128), 'LineWidth', 1, ...
                'LineStyle', '-.', 'Color', [0 0 0]);
            plot(axes(2), (f1 + f2/2)*ones(1, 128), ...
                linspace(0, 1, 128), 'LineWidth', 1, ...
                'LineStyle', '-.', 'Color', [0 0 0]);
            hold(axes(2), 'off');
        case 'hpf'
            deleteplot(axes(2), [1 0 0]);
            Ef = filterwaveform(data{6}, mean(diff(data{2})),...
                'hpf', f1, f2); % f1: cutoff f2: pass

            hold(axes(2), 'on');
            plot(axes(2), f2*ones(1, 128), linspace(0, 1, 128), ...
                'LineWidth', 1, 'LineStyle', '-.', 'Color', 'r');
            hold(axes(2), 'off');
        case 'lpf'
            deleteplot(axes(2), [0 1 0]);
            Ef = filterwaveform(data{6}, mean(diff(data{2})),...
                'lpf', f1, f2); % f1: pass f2: cutoff

            hold(axes(2), 'on');
            plot(axes(2), ones(1, 128), linspace(0, 1, 128), ...
                'LineWidth', 1, 'LineStyle', '-.', 'Color', 'g');
            hold(axes(2), 'off');
    end

    switch comp
        case 'all'
            plot(axes(1), data{2}, Ef(:, 1), 'r');hold(axes(1), 'on')
            plot(axes(1), data{2}, Ef(:, 2), 'g');
            plot(axes(1), data{2}, Ef(:, 3), 'b');hold(axes(1), 'off')
            legend(axes(1), 'E_B', 'E_{B\timesV}', 'E_{B\times(B\timesV)}')
        case 'x'
            plot(axes(1), data{2}, Ef(:, 1), 'r');
            legend(axes(1), 'E_B')
        case 'y'
            plot(axes(1), data{2}, Ef(:, 2), 'r');
            legend(axes(1), 'E_{B\timesV}')
        case 'z'
            plot(axes(1), data{2}, Ef(:, 3), 'b');
            legend(axes(1), 'E_{B\times(B\timesV)}')
    end
    
    xlim(axes(1), [data{2}(1) data{2}(end)])

    
function plotstftspec(data, axes, comp)

dt = double(mean(diff(data{2})));

if strcmp(comp, 'all') || strcmp(comp, 'x')
    I = 1;
elseif strcmp(comp, 'y')
    I = 2;
else
    I = 3;
end

[s,f,t] = spectrogram(data{6}(:,I),...
    chebwin(1024), 128, 1024, 1/dt);

surf(axes, t, f, log10(abs(s)));
shading(axes, 'interp');
view(axes, 2);

xlim(axes, [t(1) t(end)]);
ylim(axes, [f(1) f(end)])

% --- Outputs from this function are returned to the command line.
function varargout = TDSviewer_m4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fcen = get(handles.slider1, 'Value') + get(handles.slider2, 'Value');
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');

    drawfields(handles.data, [handles.axes1, handles.axes2],...
        fcen, bw, 'bpf', handles.comp);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fcen = get(handles.slider2, 'Value') + get(handles.slider1, 'Value');
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');

    drawfields(handles.data, [handles.axes1, handles.axes2],...
        fcen, bw, 'bpf', handles.comp);



% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fcen = get(handles.slider2, 'Value') + get(handles.slider1, 'Value');
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');

    drawfields(handles.data, [handles.axes1, handles.axes2],...
        fcen, bw, 'bpf', handles.comp);




% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fcen = get(handles.slider2, 'Value') + get(handles.slider1, 'Value');
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');

    drawfields(handles.data, [handles.axes1, handles.axes2],...
        fcen, bw, 'bpf', handles.comp);




% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject, 'tag')
    case 'radiobutton1'
        set(handles.slider1, 'Enable', 'off')
        set(handles.slider2, 'Enable', 'off')
        set(handles.slider3, 'Enable', 'off')
        set(handles.slider4, 'Enable', 'off')
        set(handles.slider5, 'Enable', 'off')
        set(handles.slider6, 'Enable', 'off')
        
        drawfields(handles.data, [handles.axes1, handles.axes2],...
            0, 0, 'ori', handles.comp);
        handles.mode = 'ori';
    case 'radiobutton2'
        set(handles.slider1, 'Enable', 'on')
        set(handles.slider2, 'Enable', 'on')
        set(handles.slider3, 'Enable', 'on')
        set(handles.slider4, 'Enable', 'on')
        set(handles.slider5, 'Enable', 'off')
        set(handles.slider6, 'Enable', 'off')
        
        fcen = get(handles.slider1, 'Value') + get(handles.slider2, 'Value');
        bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');

        drawfields(handles.data, [handles.axes1, handles.axes2],...
            fcen, bw, 'bpf', handles.comp);
        handles.mode = 'bpf';
    case 'radiobutton3'
        set(handles.slider1, 'Enable', 'off')
        set(handles.slider2, 'Enable', 'off')
        set(handles.slider3, 'Enable', 'off')
        set(handles.slider4, 'Enable', 'off')
        set(handles.slider5, 'Enable', 'off')
        set(handles.slider6, 'Enable', 'on')
        
        stb = get(handles.slider6, 'Value');
        
        drawfields(handles.data, [handles.axes1, handles.axes2],...
            stb, stb + 1, 'hpf', handles.comp);
        handles.mode = 'hpf';
    case 'radiobutton4'
        set(handles.slider1, 'Enable', 'off')
        set(handles.slider2, 'Enable', 'off')
        set(handles.slider3, 'Enable', 'off')
        set(handles.slider4, 'Enable', 'off')
        set(handles.slider5, 'Enable', 'on')
        set(handles.slider6, 'Enable', 'off')
        
        stb = get(handles.slider6, 'Value');
        
        drawfields(handles.data, [handles.axes1, handles.axes2],...
            stb, stb + 1, 'lpf', handles.comp);
        handles.mode = 'lpf';
end

guidata(hObject, handles)

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

name = get(hObject, 'String');
name = name(12:end);

switch name
    case 'All'
        handles.comp = 'x';
        set(hObject, 'String', 'Component: x');
    case 'x'
        handles.comp = 'y';
        set(hObject, 'String', 'Component: y');
    case 'y'
        handles.comp = 'z';
        set(hObject, 'String', 'Component: z');
    case 'z'
        handles.comp = 'all';
        set(hObject, 'String', 'Component: All');
end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if handles.eventindex > 1
        handles.eventindex = handles.eventindex - 1;
        handles.data = spdfcdfread(handles.names_tds{handles.eventindex});
    else
        return;
    end
    
    [f, I, pE] = findmainfreq(handles.data{6}, mean(diff(handles.data{2})));
    plot(handles.axes2, f, pE);
    
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');
    
    set(handles.slider2, 'Value', f(I));
    set(handles.slider1, 'Value', 0);
    
    switch handles.mode
        case 'ori'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                0, 0, 'ori', handles.comp);
        case 'bpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                f(I), bw, 'bpf', handles.comp);
        case 'hpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'hpf', handles.comp);
        case 'lpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'lpf', handles.comp);
    end
    
    plotstftspec(handles.data, handles.axes3, handles.comp);
    
    set(handles.edit1, 'String', num2str(handles.eventindex));
    
    guidata(hObject, handles)

    
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if handles.eventindex < length(handles.names_tds)
        handles.eventindex = handles.eventindex + 1;
        handles.data = spdfcdfread(handles.names_tds{handles.eventindex});
    else
        return;
    end
    
    [f, I, pE] = findmainfreq(handles.data{6}, mean(diff(handles.data{2})));
    plot(handles.axes2, f, pE);
    
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');
    
    set(handles.slider2, 'Value', f(I));
    set(handles.slider1, 'Value', 0);
    
    switch handles.mode
        case 'ori'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                0, 0, 'ori', handles.comp);
        case 'bpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                f(I), bw, 'bpf', handles.comp);
        case 'hpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'hpf', handles.comp);
        case 'lpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'lpf', handles.comp);
    end
    
    plotstftspec(handles.data, handles.axes3, handles.comp);
    
    set(handles.edit1, 'String', num2str(handles.eventindex));
    
    guidata(hObject, handles)

    
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    eventindex = str2double(get(handles.edit1, 'String'));
    
    if isnan(eventindex)
        return
    elseif eventindex >= 1 && eventindex <= length(handles.names_tds)
        handles.eventindex = eventindex;
        handles.data = spdfcdfread(handles.names_tds{handles.eventindex});
    else
        return
    end
    
    [f, I, pE] = findmainfreq(handles.data{6}, mean(diff(handles.data{2})));
    plot(handles.axes2, f, pE);
    
    bw = get(handles.slider4, 'Value') + get(handles.slider3, 'Value');
    
    set(handles.slider2, 'Value', f(I));
    set(handles.slider1, 'Value', 0);
    
    switch handles.mode
        case 'ori'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                0, 0, 'ori', handles.comp);
        case 'bpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                f(I), bw, 'bpf', handles.comp);
        case 'hpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'hpf', handles.comp);
        case 'lpf'
            drawfields(handles.data, [handles.axes1, handles.axes2],...
                stb, stb + 1, 'lpf', handles.comp);
    end
    
    plotstftspec(handles.data, handles.axes3, handles.comp);
    
    guidata(hObject, handles)
    
    
