function varargout = TDStimeanalyzer_m2(varargin)
% TDSTIMEANALYZER_M2 MATLAB code for TDStimeanalyzer_m2.fig
%      TDSTIMEANALYZER_M2, by itself, creates a new TDSTIMEANALYZER_M2 or raises the existing
%      singleton*.
%
%      H = TDSTIMEANALYZER_M2 returns the handle to a new TDSTIMEANALYZER_M2 or the handle to
%      the existing singleton*.
%
%      TDSTIMEANALYZER_M2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TDSTIMEANALYZER_M2.M with the given input arguments.
%
%      TDSTIMEANALYZER_M2('Property','Value',...) creates a new TDSTIMEANALYZER_M2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TDStimeanalyzer_m2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TDStimeanalyzer_m2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TDStimeanalyzer_m2

% Last Modified by GUIDE v2.5 04-May-2022 15:26:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TDStimeanalyzer_m2_OpeningFcn, ...
                   'gui_OutputFcn',  @TDStimeanalyzer_m2_OutputFcn, ...
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


% --- Executes just before TDStimeanalyzer_m2 is made visible.
function TDStimeanalyzer_m2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TDStimeanalyzer_m2 (see VARARGIN)

% Choose default command line output for TDStimeanalyzer_m2
handles.output = hObject;

set(gcf,'toolbar','figure')

set(handles.axes1, 'box', 'on');
set(handles.axes2, 'box', 'on');
set(handles.axes3, 'box', 'on');
set(handles.axes4, 'box', 'on');
set(handles.edit2, 'enable', 'off');

selector = varargin{1};
path = varargin{2};

if selector
    selector_1 = 'STA';
    selector_2 = 'sta';
else
    selector_1 = 'STB';
    selector_2 = 'stb';
end


% --- Rotation matrix from SC to RTN

load([selector_2, '_sc2rtn.mat'])

handles.sc2rtn = zeros(3, 3, length(m00));

handles.sc2rtn(1, 1, :) = reshape(m00, [1 1 length(m00)]);
handles.sc2rtn(1, 2, :) = reshape(m01, [1 1 length(m00)]);
handles.sc2rtn(1, 3, :) = reshape(m02, [1 1 length(m00)]);
handles.sc2rtn(2, 1, :) = reshape(m10, [1 1 length(m00)]);
handles.sc2rtn(2, 2, :) = reshape(m11, [1 1 length(m00)]);
handles.sc2rtn(2, 3, :) = reshape(m12, [1 1 length(m00)]);
handles.sc2rtn(3, 1, :) = reshape(m20, [1 1 length(m00)]);
handles.sc2rtn(3, 2, :) = reshape(m21, [1 1 length(m00)]);
handles.sc2rtn(3, 3, :) = reshape(m22, [1 1 length(m00)]);

handles.time_in_min_rtn = time_in_min/60;


% --- Rotation matrix from SC to GSE (for wind)

load([selector_2, '_sc2gse.mat'])

handles.sc2gse = zeros(3, 3, length(m00));

handles.sc2gse(1, 1, :) = reshape(m00, [1 1 length(m00)]);
handles.sc2gse(1, 2, :) = reshape(m01, [1 1 length(m00)]);
handles.sc2gse(1, 3, :) = reshape(m02, [1 1 length(m00)]);
handles.sc2gse(2, 1, :) = reshape(m10, [1 1 length(m00)]);
handles.sc2gse(2, 2, :) = reshape(m11, [1 1 length(m00)]);
handles.sc2gse(2, 3, :) = reshape(m12, [1 1 length(m00)]);
handles.sc2gse(3, 1, :) = reshape(m20, [1 1 length(m00)]);
handles.sc2gse(3, 2, :) = reshape(m21, [1 1 length(m00)]);
handles.sc2gse(3, 3, :) = reshape(m22, [1 1 length(m00)]);

handles.time_in_min_gse = time_in_min/60;


% --- Direction cosine of solar wind speed

load([selector_2, '_vrtn.mat'])
load([selector_2, '_wind_vgse.mat'])

handles.tt_P_wind = tt_P_wind;
handles.vgse = cat(2, X, Y, Z);
handles.vrtn = cat(2, R, T, N);

handles.time_in_min_vdc = time_in_min;


% --- Equivalent antenna parameters

load antenna_90pf.mat;
M2i = zeros(3);

M2i(1, 1) = hx(1)*cos(hx(2));
M2i(2, 1) = hy(1)*cos(hy(2));
M2i(3, 1) = hz(1)*cos(hz(2));

M2i(1, 2) = hx(1)*sin(hx(2))*sin(hx(3));
M2i(2, 2) = hy(1)*sin(hy(2))*sin(hy(3));
M2i(3, 2) = hz(1)*sin(hz(2))*sin(hz(3));

M2i(1, 3) = -hx(1)*sin(hx(2))*cos(hx(3));
M2i(2, 3) = -hy(1)*sin(hy(2))*cos(hy(3));
M2i(3, 3) = -hz(1)*sin(hz(2))*cos(hz(3));

handles.M2i = M2i;


% --- Timestamps of waveforms

load([selector_2, '_moy.mat'])

handles.moy = moy;


% --- Waveform data

path = [path, '\', selector_1, '_TDS'];
files = dir(path);
files = files(3:end);

names_tds = cell(size(files));
names_mag = cell(size(files));
for i = 1:length(names_tds)
    names_tds{i} = [path, '\', files(i).name, '\', ...
        selector_2, '_l2_tds_lang_', files(i).name(10:24), '_v01.cdf'];
    names_mag{i} = [path, '\', files(i).name, '\', ...
        selector_2, '_l2_mag_lang_', files(i).name(10:24), '_v01.cdf'];    
end

set(handles.edit2, 'String', num2str(length(names_tds)));
set(handles.edit1, 'String', '1');

handles.names_tds = names_tds;
handles.names_mag = names_mag;


% --- Solar wind background parameters

load([selector_2, '_solarwindpara.mat'])
load([selector_2, '_wind_swpara.mat'])

handles.thermalene = zeros(size(moy));
moy_wind = moy(moy <= max(tt_P_wind));
moy_ster = moy(moy > max(tt_P_wind));

kb = 1.3806505e-23;

handles.thermalene(moy <= tt_P_wind(end)) = ...
    (kb*Tp_wind(moy_wind + 1)).*(Np_wind(moy_wind + 1)*1e6);
handles.thermalene(moy > tt_P_wind(end)) = ...
    (kb*Tp(moy_ster - tt_P(1) + 1)).*(Np(moy_ster - tt_P(1) + 1)*1e6);

handles.tt_P_wind_tene = tt_P_wind;

handles.globalindex = 1;

DrawFigure(handles, 1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TDStimeanalyzer_m2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function [Ef, E, t] = loadwaveform(filename)
% Description:
%
% Load field data and Apply Butterworth BPF on electric field time series
%
% Input arguments:
%
% filename: data file of waveform
%
% Output arguments:
%
% Ef: bandpassed electric field
% E: original electric field
% t: time of observed electric field sequence
%

data = spdfcdfread(filename);

if mean(diff(data{2})) < 0.005 % Determine sample rate
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

Ef(:, 1) = filter(bb, ab, double(data{6}(:, 1)));
Ef(:, 2) = filter(bb, ab, double(data{6}(:, 2)));
Ef(:, 3) = filter(bb, ab, double(data{6}(:, 3)));

E = data{6};

t = data{2};

function [B, t] = loadmagfield(filename)
% Description:
%
% Load magnetic field in SC
%
% Input arguments:
%
% filename: data file of waveform
%
% Output arguments:
%
% B: vector magnetic field
% t: time of observed magnetic field sequence
%

data = spdfcdfread(filename);

B = data{3};
t = data{2};

function swdc = sw_coord2sc(swdc, sc2coord)
% Description:
%
% Transform solar wind direction cosines from other coordinates to SC
%
% Input arguments:
%
% swdc: solar wind direction cosines
% sc2coord: rotation matrix from SC to other coordinates
%
% Output arguments:
%
% swdc: solar wind direction cosines in SC frame
%

coord2sc = sc2coord^-1;

if norm(size(swdc) - [1 3]) == 0
    swdc = coord2sc*swdc';
else
    swdc = coord2sc*swdc;
end

swdc = swdc/norm(swdc);


function efield_sc = efield_fa2sc(efield, B, swdc)
% Description:
%
% Transform electric field vector from FA frame to SC frame
%
% Input arguments:
%
% efield: electric field vector
% B: magnetic field vector
% swdc: solar wind direction cosines
%
% Output arguments:
%
% efield_sc: electric field vector in SC frame
%

e1 = B/norm(B);
e2 = cross(e1, swdc)';

e2 = e2/norm(e2);

e3 = cross(e1, e2)';

MR = [[e1(1) e2(1) e3(1)];[e1(2) e2(2) e3(2)];[e1(3) e2(3) e3(3)]];

efield_sc = MR*efield';

efield_sc = efield_sc';


function phi = efiled2potential(efield_sc, M2i)
% Description:
%
% Transform electric field in SC frame into floating potential
%
% Input arguments:
%
% efield_sc: electric field vector in SC frame
% M2i: equivalent antenna parameters
%
% Output arguments:
%
% phi: floating potential on antenna
%

phi = M2i*efield_sc';
phi = phi';


function data_low = extractlowfreq(data)
% Description:
%
% Apply Gaussian LPF on time series
%
% Input arguments:
%
% data: time series
%
% Output arguments:
%
% data_low: lowpassed time series
%

width = 30;
centre = 20;
filt_b = exp(-((1:16384) - (8192 - centre)).^2/width.^2) + ...
    exp(-((1:16384) - (8193 + centre)).^2/width.^2);

filt_b = filt_b'/max(filt_b);

data_fft = fftshift(fft(data) / length(data));
data_low = real(ifft(fftshift(data_fft.*filt_b))) * length(data);


function DrawFigure(handles, index)
% Description:
%
% Draw bandpassed waveform, envelope of field energy and density 
% fluctuation and their frequency spectra
%
% Input arguments:
%
% handles: GUI handles
% index: index of event
%
% Output arguments:
%
% none
%

[Ef, E, t] = loadwaveform(handles.names_tds{index});
[B, ~] = loadmagfield(handles.names_mag{index});

set(handles.edit1, 'String', num2str(index));
plot(handles.axes1, t, Ef);

ylabel(handles.axes1, 'E-field (mV/m)')

% --- Calculate solar wind direction cosines
if handles.moy(index) < handles.tt_P_wind(end)
    swdc = sw_coord2sc(handles.vgse(handles.moy(index) + 1, :), ...
        handles.sc2gse(:, :, handles.moy(index) + 1));
else
    swdc = sw_coord2sc(handles.vrtn(handles.moy(index) - handles.time_in_min_vdc(1) + 1, :), ...
        handles.sc2rtn(:, :, handles.moy(index) + 1));
end

% --- Calculate floating potential
efield_sc = efield_fa2sc(E, B(floor(length(B(:, 1))/2), :), swdc);
phi = efiled2potential(efield_sc, handles.M2i);

phi_low(:, 1) = extractlowfreq(phi(:, 1));
phi_low(:, 2) = extractlowfreq(phi(:, 2));
phi_low(:, 3) = extractlowfreq(phi(:, 3));

plot(handles.axes2, t, phi_low(:, 1), 'r', 'LineWidth', 1); hold(handles.axes2, 'on');
plot(handles.axes2, t, phi_low(:, 2), 'g', 'LineWidth', 1);
plot(handles.axes2, t, phi_low(:, 3), 'b', 'LineWidth', 1); hold(handles.axes2, 'off');

ylabel(handles.axes2, 'Potential (mV)')

% --- Compare field energy and density fluctuation
[~, I] = max(var(phi_low, 1));

W = (Ef(:, 1).^2 + Ef(:, 2).^2 + Ef(:, 3).^2)*1e-6*8.8541878e-12/2;
W = W/handles.thermalene(index);

dn = phi_low(:, I)*1e-3/3;

yyaxis(handles.axes3, 'right')
plot(handles.axes3, t, dn, 'r', 'LineWidth', 1.5);
ylim(handles.axes3, [-max(abs(dn)) max(abs(dn))])

ylabel(handles.axes3, 'Normalized energy')

yyaxis(handles.axes3, 'left')
plot(handles.axes3, t, W, 'b', 'Linewidth', 1);
ylim(handles.axes3, [-max(W) max(W)])

ylabel(handles.axes3, 'Normalized \deltan/n')

xlim(handles.axes1, [t(1) t(end)])
xlabel(handles.axes1, 'Time (ms)')
xlim(handles.axes2, [t(1) t(end)])
xlabel(handles.axes2, 'Time (ms)')
xlim(handles.axes3, [t(1) t(end)])
xlabel(handles.axes3, 'Time (ms)')

% --- Calculate spectrum of energy envelope and density fluctuation
pdn = abs(fft(dn/max(abs(dn)))/length(t));
pdn = pdn(1:length(t)/2+1);
pdn(2:end-1) = 2*pdn(2:end-1);

f1 = (1/mean(diff(t)))*(0:(length(t)/2))/length(t);

[pks, locs] = findpeaks(double(W/max(W)));

t_r = linspace(t(locs(1)), t(locs(end)), 16384);
pks = interp1(t(locs), pks, t_r, 'linear');

pen = abs(fft(pks)/length(t_r));
pen = pen(1:length(t)/2+1);
pen(2:end-1) = 2*pen(2:end-1);

f2 = (1/mean(diff(t_r)))*(0:(length(t_r)/2))/length(t_r);

plot(handles.axes4, f1, pdn);hold(handles.axes4, 'on')
plot(handles.axes4, f2, pen);hold(handles.axes4, 'off')

xlabel(handles.axes4, 'Frequency (kHz)')
ylabel(handles.axes4, 'Power density')

set(handles.axes4, 'XScale', 'log');



% --- Outputs from this function are returned to the command line.
function varargout = TDStimeanalyzer_m2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Pushbutton: Previous
% --- Function: Go to previous event
function pushbutton1_Callback(hObject, eventdata, handles)


if handles.globalindex > 1
    handles.globalindex = handles.globalindex - 1;
    index = handles.globalindex;
else
    return;
end

if isfield(handles, 'indlist')
    index = handles.indlist(handles.globalindex);
    set(handles.edit3, 'String', num2str(handles.globalindex));
end

DrawFigure(handles, index);

guidata(hObject, handles);


% --- Pushbutton: Next
% --- Function: Go to next event
function pushbutton2_Callback(hObject, eventdata, handles)

 if isfield(handles, 'indlist')
    if handles.globalindex < length(handles.indlist)
        handles.globalindex = handles.eventindex + 1;
    else
        return;
    end
    index = handles.indlist(handles.globalindex);
else
    if handles.globalindex < length(handles.names_tds)
        handles.globalindex = handles.globalindex + 1;
    else
        return;
    end
    index = handles.globalindex;
end

DrawFigure(handles, index);

guidata(hObject, handles);


% --- Pushbutton: Go to
% --- Function: Go to event specified by event number
function pushbutton3_Callback(hObject, eventdata, handles)


if isfield(handles, 'indlist')
    globalindex = str2double(get(handles.edit3, 'String')); 
   
    if globalindex > length(handles.indlist) || globalindex < 1
        return;
    else
        handles.globalindex = globalindex;
    end
    index = handles.indlist(handles.globalindex);
else
    globalindex = str2double(get(handles.edit1, 'String')); 
    if globalindex > length(handles.names_tds) || eventindex < 1
        return;
    else
        handles.globalindex = globalindex;
    end
    index = handles.globalindex;
end

DrawFigure(handles, index);

guidata(hObject, handles);


% --- Pushbutton: Load list
% --- Function: Load an beat-event list to analyze specified events
function pushbutton4_Callback(hObject, eventdata, handles)
% Description:
%
% Variable below must be included in file:
% indlist: indexs of each event
%

[file,path] = uigetfile;

load([path, file])

handles.indlist = indlist;

guidata(hObject, handles)



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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
% --- Function: reserved for further usage
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ps1 = zeros(315, 1);
% ps2 = zeros(315, 1);
% 
% for i = 1:315
%     index = handles.indlist(i);
%     
%     [~, E, t] = loadwaveform(handles.names_tds{index});
%     [B, ~] = loadmagfield(handles.names_mag{index});
% 
%     if handles.moy(index) < handles.tt_P_wind(end)
%         swdc = sw_coord2sc(handles.vgse(handles.moy(index) + 1, :), ...
%             handles.sc2gse(:, :, handles.moy(index) + 1));
%     else
%         swdc = sw_coord2sc(handles.vrtn(handles.moy(index) - handles.time_in_min_vdc(1) + 1, :), ...
%             handles.sc2rtn(:, :, handles.moy(index) + 1));
%     end
% 
%     efield_sc = efield_fa2sc(E, B(floor(length(B(:, 1))/2), :), swdc);
%     phi = efiled2potential(efield_sc, handles.M2i);
%     disp(index)
%     phi_low(:, 1) = extractlowfreq(phi(:, 1));
%     phi_low(:, 2) = extractlowfreq(phi(:, 2));
%     phi_low(:, 3) = extractlowfreq(phi(:, 3));
% 
%     [~, I] = max(var(phi_low, 1));
%     
%     ps1(i) = max(abs(phi_low(:, I)));
%     
% end
% 
% save denflu_10.mat ps1
