function varargout = RadioGUI(varargin) 
% RADIOGUI MATLAB code for RadioGUI.fig
%
%   Applied Digital Filter Project
%   AM superheterodyne receiver
%
%   This program was written in Matlab R2015b for OS X
%   and may have some GUI issues in Windows
%
%   Demodulation of radioA.mat data set.
%
%   jorgen.hoem@gmail.com
%

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RadioGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RadioGUI_OutputFcn, ...
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
% End initialization code



% --- Visual radiotuner
function tuneRadio(freq, hObject, handles)

set(handles.sliderTune, 'value', freq);
set(handles.edCarrier, 'string', num2str(freq));

bp1min = freq-str2double(get(handles.edBP1bw, 'string'))/2;
bp1max = freq+str2double(get(handles.edBP1bw, 'string'))/2;

[x,y] = FFTRadioPower(handles.signal, handles.Fs);
plot(handles.axesSignal, x, y);
hold(handles.axesSignal, 'on')

ylim = get(handles.axesSignal,'ylim');
plot(handles.axesSignal, [freq freq], ylim, 'r');

if get(handles.checkboxBPtuning, 'value') == 1
    plot(handles.axesSignal, [bp1min bp1min], ylim, 'g');
    plot(handles.axesSignal, [bp1max bp1max], ylim, 'g');
end

hold(handles.axesSignal, 'off')
ylabel(handles.axesSignal, 'Power [dB]')
xlabel(handles.axesSignal, 'Frequency [Hz]')
title(handles.axesSignal, 'Antenna signal FFT (Fs = 6 MHz)')

handles.Ftune = freq;

set(handles.axesSignal,'buttondownfcn',{@signalAxesButtonDown,handles}) %buggy

% Update handles structure
guidata(hObject, handles);



% --- Temporary error check for filter orders.
function order = checkIfOrderIsOK(type, orderHandle)

orderChk = str2double(get(orderHandle, 'string'));

if strcmp(type, 'bp')
    if mod(orderChk,2) > 0 || orderChk == 0    
        warningMessage = sprintf('Error: BP order must be even or above 0.\nSetting order to 2.');
        uiwait(msgbox(warningMessage));
        set(orderHandle, 'string', num2str(2));
        order = 2;
    else
        order = orderChk;
    end;
else
    if orderChk == 0 || orderChk < 0
        warningMessage = sprintf('Error: Filter order must be above 0.\nSetting order to 1.');
        uiwait(msgbox(warningMessage));
        set(orderHandle, 'string', num2str(1));
        order = 1;
    else
        order = orderChk;
    end;
end;

% --- normalize/amplify amplitude vector
function y = amplify(x, amp)
if ~exist('amp', 'var')
    amp = 1;
end
max_signal = max(abs(x));
factor = amp/max_signal;
y = x.*factor;


% --- Demodulation stages.
function demodulateSignal(hObject, handles)

tuneRadio(get(handles.sliderTune, 'value'), hObject, handles);

signal_tuned = handles.signal;

% 1. Preselection bandpass filter and amplification
if get(handles.checkboxBPtuning, 'value') == 1
    bp1order = checkIfOrderIsOK('bp', handles.edBP1order);
    bp1min = handles.Ftune-str2double(get(handles.edBP1bw, 'string'))/2;
    bp1max = handles.Ftune+str2double(get(handles.edBP1bw, 'string'))/2;

    bandpass1 = designfilt('bandpassiir','FilterOrder',bp1order, ...
        'HalfPowerFrequency1',bp1min,'HalfPowerFrequency2',bp1max, ...
        'SampleRate', handles.Fs);
    signal_tuned = filtfilt(bandpass1, handles.signal);
end

signal_tuned = amplify(signal_tuned, handles.amp);

% 2. Local ocillator mixer
Fif = str2double(get(handles.edIF, 'string'));

Fmixer = handles.Ftune-Fif;

t = (1:length(signal_tuned))/handles.Fs;

signal_mixed = signal_tuned.*cos(2*pi*Fmixer*t);

handles.Flo = Fmixer;
handles.Fif = Fif;

handles.signalMixer = signal_mixed;

% 3. IF bandpass filter
bp2order = checkIfOrderIsOK('bp', handles.edBP2order);
bp2min = Fif-str2double(get(handles.edBP2bw, 'string'))/2;
bp2max = Fif+str2double(get(handles.edBP2bw, 'string'))/2;

bandpass2 = designfilt('bandpassiir','FilterOrder',bp2order, ...
    'HalfPowerFrequency1',bp2min,'HalfPowerFrequency2',bp2max, ...
    'SampleRate', handles.Fs);

signal_IF = filtfilt(bandpass2, signal_mixed);

handles.signalBandpass = signal_IF;

% 4. IF amplifier
signal_AM = amplify(abs(signal_IF), handles.amp);

handles.signalIFamplifier = signal_AM;

% 5. Low pass demodulator
% todo: demodulate offset

lppass = str2double(get(handles.edLPpass, 'string'));
lpripple = str2double(get(handles.edLPripple, 'string'));
lporder = checkIfOrderIsOK('lp', handles.edLPorder);

lowpass1 = designfilt('lowpassiir','FilterOrder',lporder, ...
         'PassbandFrequency',lppass,'PassbandRipple',lpripple, ...
         'SampleRate', handles.Fs);

signal_LP = filtfilt(lowpass1, signal_AM);

% additional LP filtering
if get(handles.checkboxLP2active, 'value') == 1
    %disp('using LP2')
    lp2pass = str2double(get(handles.edLP2pass, 'string'));
    lp2ripple = str2double(get(handles.edLP2ripple, 'string'));
    lp2order = checkIfOrderIsOK('lp', handles.edLP2order);

    lowpass2 = designfilt('lowpassiir','FilterOrder',lp2order, ...
         'PassbandFrequency',lp2pass,'PassbandRipple',lp2ripple, ...
         'SampleRate', handles.Fs);

    signal_LP = filter(lowpass2, signal_LP);
end;

handles.signalDemodulated = signal_LP;

% 6. LF amplifier
handles.signalLFamplifier = amplify(signal_LP, handles.amp);

% Update handles structure
guidata(hObject, handles);

% Update selected view
updateView(hObject, handles);



% --- Update selected view.
function updateView(hObject, handles)

cView = handles.currentView;

if strcmp(cView, 'antenna')
    antennaBtn_Callback(hObject, 0, handles);
elseif strcmp(cView, 'mixer')
    mixerBtn_Callback(hObject, 0, handles);
elseif strcmp(cView, 'bandpass')
    bandpassBtn_Callback(hObject, 0, handles);
elseif strcmp(cView, 'ifamplifier')
    ifamplifierBtn_Callback(hObject, 0, handles);
elseif strcmp(cView, 'demodulator')
    demodulatorBtn_Callback(hObject, 0, handles);
elseif strcmp(cView, 'lfamplifier')
    lfamplifierBtn_Callback(hObject, 0, handles);
end;



% --- Single sided FFT power spectrum.
function [fftx,ffty] = FFTRadioPower(sig,fs)

N = length(sig);
N_2 = ceil(N/2); % single-sided spectrum

x_hertz = (0:N-1)*fs/N;

fft_sig = abs(fft(sig)); %fftshift

fftx = x_hertz(1:N_2);
ffty = 10*log10(fft_sig(1:N_2));



% --- Display frequency at mouse position on main signal. Buggy
function [] = signalAxesCursorPos(hObject, eventdata, handles)

XLM = get(handles.axesSignal, 'xlim');
AXP = get(handles.axesSignal, 'pos');
DFX = diff(XLM);

F = get(hObject, 'currentpoint');
tf1 = AXP(1) <= F(1) && F(1) <= AXP(1) + AXP(3);
tf2 = AXP(2) <= F(2) && F(2) <= AXP(2) + AXP(4);

if tf1 && tf2
    Cx =  XLM(1) + (F(1)-AXP(1)).*(DFX/AXP(3));
    set(handles.textFreqCursor, 'visible', 'on')
    set(handles.textFreqCursor, 'string', [num2str(round(Cx)/1e3), ' kHz'])
else
    set(handles.textFreqCursor, 'visible', 'off')
end



% --- Tune radio by clicking main signal spectrum.
function signalAxesButtonDown(hObject, eventdata, handles)
%disp('click')
seltype = get(handles.figure1, 'selectiontype'); % Right-or-left click?

if strmatch(seltype,'normal')
    p = get(handles.axesSignal, 'currentpoint');
    tuneRadio(p(1), hObject, handles)
end

set(handles.axesSignal,'buttondownfcn',{@signalAxesButtonDown,handles}) % buggy2k
%guidata(hObject, handles);



% --- Load data set.
% todo: open file
function sig = loadDataSet()

if exist('radioA.mat', 'file')
	sig = importdata('radioA.mat');
else
	warningMessage = sprintf('Warning: Could not find mat-file: radioA.mat\n\nLoading dummydata');
	uiwait(msgbox(warningMessage));
    % loading dummydata
    Fs = 6e6;
    t = (1:1000) / Fs;
    Fm = cos(2*pi*4e3*t); % Message signal
    Fc1 = cos(2*pi*5000e3*t); % Carrier signal 1
    Fc2 = cos(2*pi*1700e3*t); % Carrier signal 2
    sig = ((1+Fm).*Fc1) + ((1+Fm).*Fc2); % AM signal
end



% --- Init default values before RadioGUI is made visible
function RadioGUI_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

handles.Fs = 6e6; % todo: selectable

handles.signal = loadDataSet();

[x, y] = FFTRadioPower(handles.signal, handles.Fs);
[~, loc] = findpeaks(y, 'MinPeakHeight', mean(y)*1.5, 'MinPeakDistance', 30);
peaks = round(x(loc)/1e3)*1e3;
handles.signalFFTpeaks = peaks;

set(handles.checkboxBPtuning, 'value', 1);

handles.Ftune = peaks(1);
tuneRadio(handles.Ftune, hObject, handles);

% Tuner step slider
set(handles.sliderTuneStep, 'min', 0)
set(handles.sliderTuneStep, 'max', length(peaks)-1)
set(handles.sliderTuneStep, 'sliderstep', [1/(length(peaks)-1) 1/(length(peaks)-1)]) 

% Tuner slider
Fss = handles.Fs/2;
set(handles.sliderTune, 'min', 0);
set(handles.sliderTune, 'max', Fss);
set(handles.sliderTune, 'sliderstep', [1000/(Fss-1) 1000/(Fss-1)]);

% Set defaults
set(handles.edBP1bw, 'string', '100e3');
set(handles.edBP1order, 'string', '2');

set(handles.edIF, 'string', '450e3');

set(handles.edBP2bw, 'string', '5e3');
set(handles.edBP2order, 'string', '2');

set(handles.edLPpass, 'string', '5e3');
set(handles.edLPripple, 'string', '0.2');
set(handles.edLPorder, 'string', '4');

set(handles.edLP2pass, 'string', '2e3');
set(handles.edLP2ripple, 'string', '0.4');
set(handles.edLP2order, 'string', '1');

handles.Flo = 0;
handles.Fif = 0;

handles.currentView = 'antenna';

handles.signalMixer = 0;
handles.signalBandpass = 0;
handles.signalIFamplifier = 0;
handles.signalDemodulated = 0;
handles.signalLFamplifier = 0;

handles.amp = 3.3;

demodulateSignal(hObject, handles);

updateView(hObject, handles);

set(hObject, 'windowbuttonmotionfcn', {@signalAxesCursorPos,handles})
set(handles.axesSignal,'buttondownfcn', {@signalAxesButtonDown,handles})

% Update handles structure
guidata(hObject, handles);


    
% --- Outputs from this function are returned to the command line.
function varargout = RadioGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



% I/O
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on slider movement.
function sliderTune_Callback(hObject, eventdata, handles)

tuneRadio(get(hObject,'value'), hObject, handles);
%demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sliderTune_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edCarrier_Callback(hObject, eventdata, handles)

L = get(handles.sliderTune, {'min','max','value'}); 
E = str2double(get(hObject, 'string'));
if E >= L{1} && E <= L{2}
    tuneRadio(E, hObject, handles);
    %demodulateSignal(hObject, handles)
else
    set(hObject, 'string', L{3});
end



% --- Executes during object creation, after setting all properties.
function edCarrier_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edBP1bw_Callback(hObject, eventdata, handles)

%tune = get(handles.sliderTune, 'value');
%tuneRadio(tune, hObject, handles);

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edBP1bw_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edIF_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);

%Fif = str2double(get(hObject, 'string'));
%radioTune(1000e3, hObject, handles);


% --- Executes during object creation, after setting all properties.
function edIF_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edBP2bw_Callback(hObject, eventdata, handles)

%tune = get(handles.sliderTune, 'value');
%tuneRadio(tune, hObject, handles);
demodulateSignal(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edBP2bw_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edLPpass_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edLPpass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edLPripple_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edLPripple_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edBP1order_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edBP1order_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edBP2order_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edBP2order_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushBtnDemodulate.
function pushBtnDemodulate_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);


function edLPorder_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edLPorder_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in antennaBtn.
function antennaBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signal)
title(handles.axesResult, 'Antenna signal')

[x,y] = FFTRadioPower(handles.signal, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')
ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'antenna';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: Antenna')
guidata(hObject, handles);



% --- Executes on button press in mixerBtn.
function mixerBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signalMixer)
title(handles.axesResult, ['Mixer (LO = ', num2str(handles.Flo/1e3), ' kHz)'])
ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')

bp2min = handles.Fif-str2double(get(handles.edBP2bw, 'string'))/2;
bp2max = handles.Fif+str2double(get(handles.edBP2bw, 'string'))/2;

[x,y] = FFTRadioPower(handles.signalMixer, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylim = get(handles.axesResultFFT,'ylim');

hold(handles.axesResultFFT, 'on')
plot(handles.axesResultFFT, [bp2min bp2min], ylim, 'g');
plot(handles.axesResultFFT, [handles.Fif handles.Fif], ylim, 'r');
plot(handles.axesResultFFT, [bp2max bp2max], ylim, 'g');
hold(handles.axesResultFFT, 'off')

ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'mixer';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: Mixer')
guidata(hObject, handles);



% --- Executes on button press in bandpassBtn.
function bandpassBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signalBandpass)
title(handles.axesResult, 'Bandpass filtering')

[x,y] = FFTRadioPower(handles.signalBandpass, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')
ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'bandpass';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: Bandpass filter')
guidata(hObject, handles);



% --- Executes on button press in ifamplifierBtn.
function ifamplifierBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signalIFamplifier)
title(handles.axesResult, 'IF amplified')

[x,y] = FFTRadioPower(handles.signalIFamplifier, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')
ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'ifamplifier';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: IF amplifier')
guidata(hObject, handles);



% --- Executes on button press in demodulatorBtn.
function demodulatorBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signalDemodulated)
title(handles.axesResult, 'Demodulation')

[x,y] = FFTRadioPower(handles.signalDemodulated, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')
ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'demodulator';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: Demodulator')
guidata(hObject, handles);



% --- Executes on button press in lfamplifierBtn.
function lfamplifierBtn_Callback(hObject, eventdata, handles)

plot(handles.axesResult, handles.signalLFamplifier)
title(handles.axesResult, 'LF amplified')

[x,y] = FFTRadioPower(handles.signalLFamplifier, handles.Fs);
plot(handles.axesResultFFT, x, y);
title(handles.axesResultFFT, 'FFT')

ylabel(handles.axesResult, 'Amplitude [V]')
xlabel(handles.axesResult, 'Time [s]')
ylabel(handles.axesResultFFT, 'Power [dB]')
xlabel(handles.axesResultFFT, 'Frequency [Hz]')

handles.currentView = 'lfamplifier';
set(handles.uibuttongroup1, 'title', 'Superheterodyne: LF amplifier')
guidata(hObject, handles);



% --- Executes on slider movement.
function sliderTuneStep_Callback(hObject, eventdata, handles)

ind = round(get(hObject, 'value'))+1;
freq = handles.signalFFTpeaks(ind);
tuneRadio(freq, hObject, handles);
%demodulateSignal(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sliderTuneStep_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in checkboxLP2active.
function checkboxLP2active_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);



function edLP2pass_Callback(hObject, eventdata, handles)

set(handles.checkboxLP2active, 'value', 1);
demodulateSignal(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edLP2pass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edLP2ripple_Callback(hObject, eventdata, handles)

set(handles.checkboxLP2active, 'value', 1);
demodulateSignal(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edLP2ripple_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function edLP2order_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxBPtuning.
function checkboxBPtuning_Callback(hObject, eventdata, handles)

demodulateSignal(hObject, handles);

%[EOF]
