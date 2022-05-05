function varargout = LFanalyzer(varargin)
% LFANALYZER MATLAB code for LFanalyzer.fig
%      LFANALYZER, by itself, creates a new LFANALYZER or raises the existing
%      singleton*.
%
%      H = LFANALYZER returns the handle to a new LFANALYZER or the handle to
%      the existing singleton*.
%
%      LFANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFANALYZER.M with the given input arguments.
%
%      LFANALYZER('Property','Value',...) creates a new LFANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LFanalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LFanalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LFanalyzer

% Last Modified by GUIDE v2.5 12-Feb-2022 16:46:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LFanalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @LFanalyzer_OutputFcn, ...
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


% --- Initialization of GUI
function LFanalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LFanalyzer (see VARARGIN)

% Choose default command line output for LFanalyzer

handles.output = hObject;

set(gcf,'toolbar','figure')

box(handles.axes1, 'on')
box(handles.axes2, 'on')
box(handles.axes3, 'on')

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
names = cell(size(files));
for i = 1:length(names)
    names{i} = [path, '\', files(i).name, '\', ...
        selector_2, '_l2_tds_lang_', files(i).name(10:24), '_v01.cdf'];  
end

handles.names = names;

%=>Initialize data
handles.eventindex = 1;
handles.compoindex = 1;
handles.data = spdfcdfread(handles.names{handles.eventindex});

handles.Fc = 3;
handles.Fb = 5;

[coefs, t, f] = calculate_cpxcwt(handles.data, ...
    handles.Fb, handles.Fc);

draw_figures(handles, t, f, coefs);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LFanalyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function Ef = filterwaveform(E, t)
% Description:
%
% Apply Butterworth HPF on electric field time series
%
% Input arguments:
%
% E: time series of electric field
% t: time of observed electric field sequence
%
% Output arguments:
%
% Ef: highpassed electric field
%

if mean(diff(t)) < 0.005 % Determine sample rate
    fs = 250;
else
    fs = 125;
end

wp = 4/(fs/2);
ws = 6/(fs/2);

Rp = 3;
Rs = 20;

[N, Wn] = buttord(wp, ws, Rp, Rs);
[bb, ab] = butter(N, Wn, 'high');

Ef(:, 1) = filter(bb, ab, double(E(:, 1)));
Ef(:, 2) = filter(bb, ab, double(E(:, 2)));
Ef(:, 3) = filter(bb, ab, double(E(:, 3)));


function [coefs, t, f] = calculate_cwt(data)
% Description:
%
% Calculate real wavelet transform
%
% Input arguments:
%
% data: time series of electric field
%
% Output arguments:
%
% coefs: real wavelet transform
% t: time of observed electric field sequence
% f: frequency calculated from wavelet scale
%

freq = 0.1:0.1:50; %Frequency in kHz

if mean(diff(data{2})) < 0.005
    Ts = 0.004;
    C = 203.125;
else
    Ts = 0.008;
    C = 101.5625;
end

t = data{2};

[~, I] = max(max(data{6}, [], 1));

[coefs, f] = cwt(double(data{6}(:, I)), C./freq, 'morl', Ts);


function coefs_env = getenv(coefs, f)
% Description:
%
% Calculate envelope of real wavelet transform by Hilbert transform
%
% Input arguments:
%
% coefs: real wavelet transform
% f: frequency of wavelet spectrum
%
% Output arguments:
%
% coefs_env: envelope of real wavelet transform
%

coefs_env = zeros(size(coefs));
for i = 1:length(f)
    coefs_env(i, :) = abs(hilbert(coefs(i, :)));
end


function  [coefs, t, f] = calculate_cpxcwt(data, Fb, Fc)
% Description:
%
% Calculate complex wavelet transform
%
% Input arguments:
%
% data: time series of electric field
% Fb: width of cpx wavelet basis
% Fc: central frequency of cpx wavelet basis
%

freq = [0.05:0.05:1, 1.1:0.1:50]; % Frequency in kHz

if mean(diff(data{2})) < 0.005 % Generate scale from frequency
    Ts = 0.004;
else
    Ts = 0.008;
end
Fs = 1/Ts;

t = data{2};

[~, I] = max(max(data{6}, [], 1));

nwlet = ['cmor', num2str(Fb), '-', num2str(Fc)];

[coefs, f] = cwt(double(data{6}(:, I)), (Fs*Fc)./freq, nwlet, Ts);


function [bicohe, f] = calculate_wtbicohe(coefs, t, f, gpuflag)
% Description:
%
% Calculate normalized bi-spectrum (bi-coherence) through
% wavelet spectrum
%
% Input arguments:
%
% coefs: wavelet spectrum
% t: time of observed electric field sequence
% f: frequency calculated from wavelet scale
% gpuflag: 1 for gpu computation, 0 for cpu computation
%

[F1, F2] = meshgrid(1:length(f), 1:length(f));

[rows, ~] = size(coefs);

f1f2 = F1 + F2;
f1f2(f1f2 > rows) = 1;

if gpuflag == 1
    slic1 = gpuArray(zeros(rows, rows, 128));
    slic2 = gpuArray(zeros(rows, rows, 128));
    slic3 = gpuArray(zeros(rows, rows, 128));
    grou1 = gpuArray(zeros(rows, rows, 128));
    grou2 = gpuArray(zeros(rows, rows, 128));
    grou3 = gpuArray(zeros(rows, rows, 128));
    t = gpuArray(t);
    coefs = gpuArray(coefs);
elseif gpuflag == 0
    slic1 = zeros(rows, rows, 128);
    slic2 = zeros(rows, rows, 128);
    slic3 = zeros(rows, rows, 128);
    grou1 = zeros(rows, rows, 128);
    grou2 = zeros(rows, rows, 128);
    grou3 = zeros(rows, rows, 128);
end

for i = 1:128
    for j = 1:128
        index = (i-1)*128 + j;
        coe_col = coefs(:, index);
        coe_col = triu(flip(coe_col(f1f2), 2));
        coe_col(((1:rows)-1)*rows + (1:rows)) = 0;
        coe_col = flip(coe_col, 2);
        slic1(:, :, j) = (coefs(:, index)*coefs(:, index).').*conj(coe_col);
        slic2(:, :, j) = abs(coefs(:, index)*coefs(:, index).').^2;
        slic3(:, :, j) = abs(coe_col).^2;
    end
    grou1(:, :, i) = trapz(t((i-1)*128+1:i*128), slic1, 3);
    grou2(:, :, i) = trapz(t((i-1)*128+1:i*128), slic2, 3);
    grou3(:, :, i) = trapz(t((i-1)*128+1:i*128), slic3, 3);
end

bispec = sum(grou1, 3);
bicohe = abs(bispec).^2./(sum(grou2, 3).*sum(grou3, 3));


function draw_figures(handles, t, f, coefs)
% Description:
%
% Draw wavelet spectrum and filtered waveform of electric field
%
% Input arguments:
%
% handles: handles of GUI
% t: time of observed electric field sequence
% f: frequency calculated from wavelet scale
% coefs: wavelet spectrum
%

    Ef = filterwaveform(handles.data{6}, handles.data{2});
    plot(handles.axes3, handles.data{2}, Ef);
    xlim(handles.axes3, [t(1) t(end)])
    
    xlabel(handles.axes3, 'Time (ms)')
    ylabel(handles.axes3, 'E-field (mV/m)')
    
    surf(handles.axes1, t, f, abs(coefs));
    shading(handles.axes1, 'interp');
    view(handles.axes1, 2);
    
    if strcmp(get(handles.pushbutton4, 'String'), 'Scale mode: Log')
        set(handles.axes1, 'YScale', 'log')
    end
    
    xlim(handles.axes1, [t(1) t(end)])
    ylim(handles.axes1, [f(1) f(end)])
    
    xlabel(handles.axes1, 'Time (ms)')
    ylabel(handles.axes1, 'Frequency (kHz)')
    
    set(handles.edit1, 'String', num2str(handles.eventindex))


% --- Outputs from this function are returned to the command line.
function varargout = LFanalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Edit: Event no.
% --- Function: Provide event number for Go-to button
function edit1_Callback(hObject, eventdata, handles)
% Description:
%
% This callback is used to check if the input is a valid numerical value
%

if isnan(str2double(get(hObject, 'String')))
    
end

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


% --- Pushbutton: Previous
% --- Function: Go to previous event
function pushbutton1_Callback(hObject, eventdata, handles)


    if handles.eventindex > 1

        handles.eventindex = handles.eventindex - 1;
        if isfield(handles, 'indlist')
            handles.data = spdfcdfread(handles.names{handles.indlist(handles.eventindex)});
        else
            handles.data = spdfcdfread(handles.names{handles.eventindex});
        end

    else
        return;
    end

    switch get(handles.radiobutton1, 'Value')
        case true
            [coefs, t, f] = calculate_cwt(handles.data);
            coefs = getenv(coefs, f);
        case false
            [coefs, t, f] = calculate_cpxcwt(handles.data, ...
                handles.Fb, handles.Fc);
    end
    
    draw_figures(handles, t, f, coefs);
    guidata(hObject, handles);

    
% --- Pushbutton: Next
% --- Function: Go to next event
function pushbutton2_Callback(hObject, eventdata, handles)


    if isfield(handles, 'indlist')
        if handles.eventindex < length(indlist)
            handles.eventindex = handles.eventindex + 1;
        else
            return;
        end
        evind = handles.indlist(handles.eventindex);
    else
        if handles.eventindex < length(handles.names)
            handles.eventindex = handles.eventindex + 1;
        else
            return;
        end
        evind = handles.eventindex;
    end
    handles.data = spdfcdfread(handles.names{evind});
    
    switch get(handles.radiobutton1, 'Value')
        case true
            [coefs, t, f] = calculate_cwt(handles.data);
            coefs = getenv(coefs, f);
        case false
            [coefs, t, f] = calculate_cpxcwt(handles.data, ...
                handles.Fb, handles.Fc);
    end
    
    draw_figures(handles, t, f, coefs);

    guidata(hObject, handles);


% --- Pushbutton: Go to
% --- Function: Go to event specified by event number
function pushbutton3_Callback(hObject, eventdata, handles)


    eventindex = str2double(get(handles.edit1, 'String'));

    if isfield(handles, 'indlist')
        if eventindex > length(handles.indlist) || eventindex < 1
            return;
        else
            handles.eventindex = eventindex;
        end
        evind = handles.indlist(handles.eventindex);
    else
        if eventindex > length(handles.names) || eventindex < 1
            return;
        else
            handles.eventindex = eventindex;
        end
        evind = handles.eventindex;
    end
    handles.data = spdfcdfread(handles.names{evind});

    switch get(handles.radiobutton1, 'Value')
        case true
            [coefs, t, f] = calculate_cwt(handles.data);
            coefs = getenv(coefs, f);
        case false
            [coefs, t, f] = calculate_cpxcwt(handles.data, ...
                handles.Fb, handles.Fc);
    end

    draw_figures(handles, t, f, coefs);

    guidata(hObject, handles)


% --- Pushbutton: Scale mode
% --- Function: Switches frequency scale mode of wavelet spectrum
function pushbutton4_Callback(hObject, eventdata, handles)

if strcmp(get(hObject, 'String'), 'Scale mode: Linear')
    set(handles.axes1, 'YScale', 'log')
    set(hObject, 'String', 'Scale mode: Log')
else
    set(handles.axes1, 'YScale', 'linear')
    set(hObject, 'String', 'Scale mode: Linear')
end


% --- Pushbutton: Bicoherence
% --- Function: Calculates normalized bi-spectrum by wavelet transform
function pushbutton5_Callback(hObject, eventdata, handles)

switch get(handles.radiobutton1, 'Value')
    case true
        [coefs, t, f] = calculate_cwt(handles.data);
    case false
        [coefs, t, f] = calculate_cpxcwt(handles.data, ...
            handles.Fb, handles.Fc);
end

switch get(handles.radiobutton5, 'Value')
    case true
        [bicohe, f] = calculate_wtbicohe(coefs, t, f, 1);
    case false
        [bicohe, f] = calculate_wtbicohe(coefs, t, f, 0);
end

surf(handles.axes2, f, f, bicohe);

shading(handles.axes2, 'interp');
view(handles.axes2, 2);

xlim(handles.axes2, [f(1) f(end)])
ylim(handles.axes2, [f(1) f(end)])

xlabel(handles.axes2, 'Frequency (kHz)')
ylabel(handles.axes2, 'Frequency (kHz)')

set(handles.axes2, 'YScale', 'log')
set(handles.axes2, 'XScale', 'log')


% --- Edit: Frequency bandwidth
% --- Function: Sets the bandwidth of bandpass filter
function edit2_Callback(hObject, eventdata, handles)
% Description:
%
% Fb: width of cpx wavelet basis
%

handles.Fb = str2double(get(hObject, 'String'));
guidata(hObject, handles)



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


% --- Edit: Central frequency
% --- Function: Sets the central frequency of bandpass filter
function edit3_Callback(hObject, eventdata, handles)
% Description:
%
% Fc: Central frequency of cpx wavelet basis
%

handles.Fc = str2double(get(hObject, 'String'));
guidata(hObject, handles)




% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Pushbutton: Load list
% --- Function: Load an beat-event list to analyze specified events
function pushbutton6_Callback(hObject, eventdata, handles)
% Description:
%
% Variable below must be included in file:
% indlist: indexs of each event
%

[file,path] = uigetfile;

load([path, file])

handles.indlist = indlist;

guidata(hObject, handles)


% --- Pushbutton: Automatic
% --- Function: Automatically calculates wavelet spectrum and bi-spectrum,
% --- and save them as pictures
function pushbutton7_Callback(hObject, eventdata, handles)
% Description:
%
% Picture will be saved in .png format. Folder contains
% pictures is specified by path
%

path = uigetdir;

for i = 1:length(handles.fsep)
    
    data = spdfcdfread(handles.names{handles.indlist(i)});
    
    %========Calculate wavelet transform and bi-spectrum==========
    switch get(handles.radiobutton1, 'Value')
        case true
            [coefs, t, f] = calculate_cwt(data);
        case false
            [coefs, t, f] = calculate_cpxcwt(data, ...
                handles.Fb, handles.Fc);
    end
    
    switch get(handles.radiobutton5, 'Value')
        case true
            [bicohe, f] = calculate_wtbicohe(coefs, t, f, 1);
        case false
            [bicohe, f] = calculate_wtbicohe(coefs, t, f, 0);
    end

    figure(11);

    %========Plot wavelet spectrum=========
    subplot(1, 2, 1);
    surf(t, f, abs(coefs)); shading interp; view(2);

    set(gca, 'YScale', 'log');

    xlim([t(1) t(end)])
    ylim([f(1) f(end)])
    
    xlabel('Time (ms)')
    ylabel('Frequency (kHz)')
    
    set(gca, 'YScale', 'log')
    set(gca, 'FontSize', 14)
    set(gca, 'FontWeight', 'bold')
    set(gca, 'LineWidth', 1.5)
    set(gca, 'FontName', 'Courier')
    box on

    %=======Plot wavelet bicoherence=========
    
    subplot(1, 2, 2);
    surf(f, f, bicohe); shading interp; view(2);

    xlim([f(1) f(end)])
    ylim([f(1) f(end)])
    
    xlabel('Frequency (kHz)')
    ylabel('Frequency (kHz)')

    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    set(gca, 'FontSize', 14)
    set(gca, 'FontWeight', 'bold')
    set(gca, 'LineWidth', 1.5)
    set(gca, 'FontName', 'Courier')
    box on
    
    %=======Save the figures as png pictures=======
    f = getframe(gcf);
    imwrite(f.cdata,[path, '\', num2str(i), '.png']);

end
