function varargout = Enge1(varargin)
% ENGE1 M-file for Enge1.fig
%      ENGE1, by itself, creates a new ENGE1 or raises the existing
%      singleton*.
%
%      H = ENGE1 returns the handle to a new ENGE1 or the handle to
%      the existing singleton*.
%
%      ENGE1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENGE1.M with the given input arguments.
%
%      ENGE1('Property','Value',...) creates a new ENGE1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Enge1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Enge1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Enge1

% Last Modified by GUIDE v2.5 23-Oct-2014 19:33:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Enge1_OpeningFcn, ...
                   'gui_OutputFcn',  @Enge1_OutputFcn, ...
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

% --- Executes just before Enge1 is made visible.
function Enge1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Enge1 (see VARARGIN)

% Choose default command line output for Enge1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Enge1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Enge1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*S.wav');
cd(pathname);
if ~isequal(filename, 0)
    handles.nSwallow = 0;
    
    handles.fname = filename(1:7);
    
    % .wavƒtƒ@ƒCƒ‹‚©‚ç?A“ª?iš‹‰º?j‰¹ƒf?[ƒ^‚ð“Ç‚Ý?ž‚Þ?B
    % [Y,handles.Fs,NBITS] = wavread(['S' handles.fname]);
     [Y,Fs,NBITS] = wavread([handles.fname 'S']);
    handles.SW1 = Y(:,1) * 10;      % ?A“ª?iš‹‰º?j‰¹
    % handles.SW1 = Y(3000000:4000000,1) * 10;      % ?A“ª?iš‹‰º?j‰¹
    handles.Fs_sound = Fs;

    % .wavƒtƒ@ƒCƒ‹‚©‚ç?ã?œƒf?[ƒ^‚ð“Ç‚Ý?ž‚Þ?B
    % [Y,handles.Fs,NBITS] = wavread(['D' handles.fname]);
     [Y,Fs,NBITS] = wavread([handles.fname 'D']);
    % handles.SW2 = - Y(:,1) * 10;     
    handles.SW2 = Y(:,1) * 10;      % debug 140506   
    % handles.SW2 = Y(300000:400000,1) * 10;      % debug 140506   
    handles.Fs_hyoid = Fs;

    % .wavƒtƒ@ƒCƒ‹‚©‚çŒÄ‹zƒf?[ƒ^‚ð“Ç‚Ý?ž‚Þ?B
    % [Y,handles.Fs,NBITS] = wavread(['B' handles.fname]);
     [Y,Fs,NBITS] = wavread([handles.fname 'B']);
    handles.flow = Y(:,1) * 10;   
    % handles.flow = Y(300000:400000,1) * 10;   
    handles.Fs_breath = Fs;

    trigger_sel_index = get(handles.trigger_signal, 'Value');
    if trigger_sel_index == 1
        % .wavƒtƒ@ƒCƒ‹‚©‚çƒgƒŠƒK?[ƒf?[ƒ^‚ð“Ç‚Ý?ž‚Þ?B
        [Y,Fs,NBITS] = wavread([handles.fname 'A']);
        handles.trigger = Y(:,1) * 10 + 10;
        handles.Fs_trigger = Fs;
    end

    handles.Fs =handles.Fs_sound;
    % handles.Fs = 1 * 10^4; % sampling frequency is 10kHz
    % handles.Fs_breath = 1 * 10^3;       % for murata: sampling frequency of breath = 1kHz 
    % handles.Fs_hyoid = 1 * 10^3;        % for murata: sampling frequency of hyoid = 1kHz 

    handles.SW1 = detrend(handles.SW1, 'constant');     % fft‚È‚Ç‚ð?s‚¤‚½‚ß?A•½‹Ï‚ð0‚É‚·‚é?B

    highpass_sel_index = get(handles.highpass_filter, 'Value');
    if highpass_sel_index == 1
        % highpass filter; % debug 140505
        Hd = HP1;
        handles.SW1 = filter(Hd, handles.SW1) * 10;
    end

    % 1kHz sampling‚É•ÏŠ·
    handles.SW11 = zeros(length(handles.SW2), 1);
    for iter = 1 : length(handles.SW2)-1    
        handles.SW11(iter) = mean(handles.SW1(10*(iter-1)+1:10*iter));
    end
    handles.SW11(length(handles.SW2)) = mean(handles.SW1(10*(length(handles.SW2)-1)+1:end));

    handles.SW2 = detrend(handles.SW2, 'constant');
    % ?ã?œ•ÏˆÊISW‚ðŒvŽZ?B
    % ?ã?œ•ÏˆÊ?i‘º“c?»ƒZƒ“ƒT‚Í?A•ÎˆÊ‘¬“x‚Å?o—Í‚³‚ê‚é?B?Ï•ª‚µ‚½‚Æ‚«‚É?A?ã•ûˆÚ“®‚ª?³’l‚ð‚Æ‚é‚æ‚¤‚É‚µ‚Ä‚¢‚é?B
    handles.ISW2 = zeros(length(handles.SW2), 1);
    % handles.ISW2(1) = -handles.SW2(1);
    handles.ISW2(1) = handles.SW2(1);       % debug 140506
    for iter = 1 : length(handles.SW2) - 1
        % handles.ISW2(iter + 1) = handles.ISW2(iter) - handles.SW2(iter + 1);
        handles.ISW2(iter + 1) = handles.ISW2(iter) + handles.SW2(iter + 1);        % debug 140506
    end

    % handles.flow = detrend(handles.flow, 'constant');
    % handles.flow = handles.flow - 2.5; % handles.flow=0 is 2.5V in Nagano-Keiki KL17 differential transmitter.

    % ‘S‘Ì‚ÌŽžŠÔŒo‰ß•\Ž¦ƒEƒBƒ“ƒhƒE‚Ì?ì?¬
     handles.t = 0 : 1 / handles.Fs : 1 / handles.Fs * (length(handles.SW1)-1);
     scrsz = get(0,'ScreenSize');
     handles.h = figure('Position',[100 50 scrsz(3)-200 scrsz(4)/2]);

     handles.h2 = figure('Position',[scrsz(3)/2+50 scrsz(4)/2 scrsz(3)/2-100 scrsz(4)/2-100]);
     set(handles.h2, 'Visible', 'off');

     handles.h3 = figure('Position',[50 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2-100]);
     set(handles.h3, 'Visible', 'off');

     figure(handles.h); subplot(4,1,1); plot(handles.t', handles.SW1);
     xlabel('Time (sec)');
     title('raw sound signal');

     handles.t_hyoid = 0 : 1 / handles.Fs_hyoid : 1 / handles.Fs_hyoid * (length(handles.SW2)-1);
     figure(handles.h); subplot(4,1,3); plot(handles.t_hyoid', handles.SW2);
     xlabel('Time (sec)');
     title('hyoid movement signal');

    % ?¶‘Ì‰¹‹æŠÔ‚ÌŒŸ?o
    % ‘S”g?®—¬‚µ‚½Œã‚ÉƒŠ?[ƒN?Ï•ª‚µ?A‰¹‚Ì‹­“x‚ÌŽžŠÔ•Ï‰»‚ðŒv‘ª‚·‚é
    % full-wave rectification
    handles.step = 0.2; % (sec)
    handles.timeScale = 0 : handles.step : max(handles.t)-1;
    d = abs(handles.SW1);
    % 'leaky' integration with time constant = tc (sec)
    % tc = 0.1;
    % for iter = 2 : length(d)
    %     d(iter) = d(iter-1) * handles.exp(-handles.Fs/tc) + d(iter);
    % end
    % compressing analog data ....
    bin_width = handles.Fs * handles.step;
    % handles.xf = zeros(1, floor(length(d) / bin_width));
    handles.xf = zeros(1, length(handles.timeScale));       % debug 140506
    for iter = 1 : length(handles.xf)
        handles.xf(iter) = mean(d(1+bin_width*(iter-1):bin_width*(iter)));       % debug 140506
    end

    % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
    figure(handles.h); subplot(4,1,1); 
    xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
    %ylim = [-2 2];
    axis([xlim ylim]);

    % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
    figure(handles.h); subplot(4,1,3); 
    xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
    %ylim = [-1 1];
    axis([xlim ylim]);

    % ŒÄ‹zƒtƒ??[ƒVƒOƒiƒ‹‚Ætrigger signal‚ð?¶‘Ì‰¹ƒŠ?[ƒN?Ï•ª‚Ìbin_width‚É‚ ‚í‚¹‚Äˆ³?k
    % % compressing analog data ....
    % bin_width = handles.Fs * handles.step;
    bin_width = handles.Fs_breath * handles.step;     % for murata: sampling frequency of breath = 1kHz 

    handles.data = zeros(1, length(handles.timeScale));
    for iter = 1 : length(handles.data)
        handles.data(iter) = mean(handles.flow(1+bin_width*(iter-1):bin_width*(iter)));      % for murata: sampling frequency of breath = 1kHz 
    end

    if trigger_sel_index == 1
        bin_width2 = handles.Fs_trigger * handles.step;     % for murata: sampling frequency of breath = 1kHz 
        handles.Trig = zeros(1, length(handles.timeScale));
        for iter = 1 : length(handles.Trig)
            handles.Trig(iter) = max(handles.trigger(1+bin_width2*(iter-1):bin_width2*(iter)));  
        end
    end

    % To detect zero flow automatically
    zeroflow_sel_index = get(handles.zero_flow, 'Value');
    switch zeroflow_sel_index
        case 1
            % fine tune handles.flow zero-level
            min_fluc = (max(handles.data) - min(handles.data)) * 0.02;
            minimum_pause = 0.6;
            for iter = 1 : length(handles.data) - minimum_pause/handles.step
                fluctuation = max(handles.data(iter:round(iter+minimum_pause/handles.step))) ...
                    - min(handles.data(iter:round(iter+minimum_pause/handles.step)));
            %     if (fluctuation < min_fluc) && (fluctuation > 0.01)
                  if (fluctuation < min_fluc)
                    min_fluc = fluctuation;
                    shift = mean(handles.data(iter:iter+minimum_pause/handles.step));
                  end
            end
            set(handles.zero_level, 'String', num2str(shift));
        case 2
            shift = str2double(get(handles.zero_level, 'String'));
            % shift = -1.86;
    end
    handles.data = handles.data - shift;
    handles.flow = handles.flow - shift;

    % ŒÄ‹zƒtƒ??[ƒVƒOƒiƒ‹‚Ì•\Ž¦
    figure(handles.h); subplot(4,1,2);
    hold off; 
    handles.tt = 0 : handles.step : (length(handles.timeScale)-1) * handles.step;
    plot(handles.tt, handles.data); 
    xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
    %ylim = [-2 2];
    axis([xlim ylim]);
    xlabel('Time (sec)');
    title('respiratory handles.flow signal');
    hold all;

    % This routine detects inspiration and expiration
    handles.insp = zeros(500, 1);
    handles.exp = zeros(500, 1);
    handles.pause = zeros(500, 2);

    handles.Fs2 = 1 / handles.step;

    idx = 1;
    ibreath = 0;
    ipause = 0;

    resp_detect_sel_index = get(handles.resp_detect, 'Value');
    switch resp_detect_sel_index
        case 1
            handles.ith = str2double(get(handles.edit7, 'String'));    % -0.15
            handles.eth = str2double(get(handles.edit8, 'String'));      % 0.15
            handles.isth = str2double(get(handles.edit9, 'String'));  % -0.08
            handles.esth = str2double(get(handles.edit10, 'String'));  % 0.05; 
        case 2
            handles.ith = min(handles.data) * 0.1;
            handles.eth = max(handles.data) * 0.1;  % 0.2
            handles.isth = min(handles.data) * 0.05;
            handles.esth = max(handles.data) * 0.05;    % 0.1
    end

    handles.ti_min = 0.3; % minimum Ti in sec
    handles.te_min = 0.5; % minimum Te in sec
    handles.ts_min = str2double(get(handles.edit4, 'String')); % minimum Ts (pause duration for a possible swallowing)

    figure(handles.h); subplot(4,1,2);
    while idx < length(handles.data)
        if ibreath==0
            while handles.data(idx) > handles.eth
                idx = idx + 1;
            end
        else
            while (handles.data(idx) > handles.eth) && ((idx - iend) / handles.Fs2 > handles.te_min)
                idx = idx + 1;
            end
        end

        % detect E-I transition
        while handles.data(idx) > handles.ith
            if idx >= length(handles.data)
                break;
            end
            idx = idx + 1;
        end
        if idx < length(handles.data)
            while (handles.data(idx) < handles.isth) && (idx > 1)
                idx = idx - 1;
            end
            idx = idx + 1;
        else
            break;
        end       
        istart = idx - 1; % modified 13/06/26

        % detect I-E transition
        while handles.data(idx) <= handles.isth    % This could be data(idx) <= 0
            idx = idx + 1;
            if idx > length(handles.data)
                    break;
            end
        end
    %     if idx > length(handles.data)
    %         break;
    %     else
            iend = idx - 1; % modified 13/06/26
            ti = (iend - istart) / handles.Fs2;
            ipeak = min(handles.data(istart:iend));
    %     end

        if (ti > handles.ti_min) && (ipeak < handles.ith)
            ibreath = ibreath + 1;
            handles.insp(ibreath) = istart;
            handles.exp(ibreath) = iend;
            y = [min(handles.data) max(handles.data)];

            % detect pause within a breath
            if ibreath > 1
                in_pause = false;
                ts = 0;
                for iter = handles.insp(ibreath - 1) + 1 : handles.insp(ibreath)
                   if (handles.data(iter) > handles.isth) && (handles.data(iter) < handles.esth)    % pause
                       if in_pause
                           ts = ts + 1 / handles.Fs2;

                           %----- for murata
                           if iter == handles.insp(ibreath)
                               pend = iter - 1; % modified 13/06/26
                               if ts >= handles.ts_min
                                   ipause = ipause + 1;
                                   handles.pause(ipause, 1) = pstart;
                                   handles.pause(ipause, 2) = pend;
                                   for pos = pstart : pend
                                       x = [pos pos] * handles.step;
                                       plot(x, y, 'y');
                                   end
                               end
                           end 
                           %----

                       else
                           in_pause = true;
                           pstart = iter - 1; % modified 13/06/26
                           ts = 1 / handles.Fs2;
                       end
                   else
                       if in_pause
                           in_pause = false;
                           pend = iter - 1; % modified 13/06/26
                           if ts >= handles.ts_min
                               ipause = ipause + 1;
                               handles.pause(ipause, 1) = pstart;
                               handles.pause(ipause, 2) = pend;
                               for pos = pstart : pend
                                   x = [pos pos] * handles.step;
                                   plot(x, y, 'y');
                               end
                           end
                           ts = 0;
                       end
                   end % if handles.pause ...
                end % for iter = handles.insp(ibreath - 1) + 1 : handles.insp(ibreath)
            end % if ibreath > 1 ...

            x = [istart istart] * handles.step;
            plot(x,y,'r');
            x = [iend iend] * handles.step;
            plot(x,y,'g');

        end % if (ti > handles.ti_min) && (ipeak < handles.ith) ...
    end % while idx <= length(handles.data) ...
    handles.nbreath = ibreath;
    handles.npause = ipause;

    % •½‹Ï‚ÌˆêŒÄ‹z‚É—v‚·‚éŽžŠÔ‚ð‹?‚ß‚é
    handles.ts_mTT = (handles.insp(handles.nbreath) - handles.insp(1)) / (handles.nbreath - 1) * handles.step;
    % This does not consider the time spent during swallowing.

    % ˆêŒÄ‹z‚ð‚P‚Æ‚µ‚Ä?A•½‹Ï‚Ì‹z‹C?\ŒÄ‹C‚Ì?Ø‚è‘Ö‚í‚è‚Ìƒ^ƒCƒ~ƒ“ƒO‚ð‹?‚ß‚é
    sTI = 0;
    for ibreath = 1 : handles.nbreath - 1
        sTI = sTI + handles.exp(ibreath) - handles.insp(ibreath);
    end
    % mTI = sTI * handles.step / (handles.nbreath - 1);
    % IE_transition = mTI / handles.ts_mTT;

    plot(handles.tt, handles.data, 'b');
    x = [0 (length(handles.timeScale)-1) * handles.step];
    y = [0 0];
    plot(x, y, 'k');

    handles.tt = 0 : handles.step : (length(handles.timeScale)-1) * handles.step;
    figure(handles.h); subplot(4,1,4);
    plot(handles.tt, handles.xf); xlim = [0 (length(handles.timeScale)-1) * handles.step]; axis([xlim ylim]);

    hold all;
    if trigger_sel_index == 1
        for iter = 1 : length(handles.Trig)
            if handles.Trig(iter) > 2.5  
                    p = [handles.tt(iter) handles.tt(iter)];
                    plot(p, ylim, 'b'); % Trigger‚ð?Â‚Å•\Ž¦       
            end
        end
    end

    % Update the handles on the interface
    guidata(handles.output,handles);

end % if ~isequal(filename, 0)

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.h)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

close all;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.h2, 'Visible', 'on');
set(handles.h3, 'Visible', 'on');

% ‚±‚±‚©‚çš‹‰º‰¹‚ðŒŸ?o‚·‚éroutine
% anal_start = 2;
% anal_end = length(handles.timeScale) - 1;

% ------------------- debug by Naomi.Yagi 140506 ----------------
% handles.xth = 0.08;     % 20140429
handles.xth = 0.03;     % 20140429

% š‹‰º‰¹‚Å‚È‚¢‚Æ—\‘z‚³‚ê‚é‹æŠÔ‚ð?í?œ‚µ‚Äè‡’lxth‚ðŒˆ’è   20140505
% pause_max = zeros(handles.npause, 1);
% pause_sum = zeros(handles.npause, 1);
% for ipause = 1 : handles.npause
%     pstart = handles.pause(ipause, 1);
%     pend = handles.pause(ipause, 2);
%     pause_max(ipause) = max(handles.xf(pstart:pend));
%     pause_sum(ipause) = sum(handles.xf(pstart:pend));
% end
% for ipause = handles.npause : (-1) : 1 
%     if pause_max(ipause) > 1 || pause_sum(ipause) > 1.2
%         pause_max(ipause) = [];
%     end
% end
% handles.xth = mean(pause_max);  
% handles.xth = mean(pause_max) - std (pause_max) * 0.4;
% handles.xth = median(pause_max);
% ------------------- debug by Naomi.Yagi 140506 ----------------

figure(handles.h); subplot(4,1,4); hold off;
plot(handles.tt, handles.xf); xlim = [0 (length(handles.timeScale)-1) * handles.step]; axis([xlim ylim]);

hold all;
trigger_sel_index = get(handles.trigger_signal, 'Value');
if trigger_sel_index == 1
    for iter = 1 : length(handles.Trig)
        if handles.Trig(iter) > 2.5  
                p = [handles.tt(iter) handles.tt(iter)];
                plot(p, ylim, 'b'); % Trigger‚ð?Â‚Å•\Ž¦       
        end
    end
end

y = ones(1,length(handles.xf)) * handles.xth;
plot(handles.tt, y);

iSound = 0;
SW_flag = zeros(1, round(length(handles.timeScale)/2));
S_start =  zeros(1, round(length(handles.timeScale)/2));
SW_param = zeros(round(length(handles.timeScale)/2), 10);
handles.SW_time = zeros(1, 100);

% melFT = zeros(melFilterNum, round(length(handles.timeScale)/2));
P_count = zeros([500 1]);
P_pos = zeros([500 100]);
P_width = zeros([500 100]);
freqVec = 50:50:2500;
PP = zeros([500 100 length(freqVec)]);

CWTpowerTh = str2double(get(handles.edit1, 'String'));       % mel-scale spectrogram %power 500-2300Hz = 20
HMampTh = str2double(get(handles.edit2, 'String'));          % set Hyoid displacement threshold HMampTh = 5;
P_count_max = str2double(get(handles.edit3, 'String'));      % set maximum pulse count P_count_max = 20;
minimal_pause = str2double(get(handles.edit4, 'String'));    % set minimum pause (sec)  minimal_pause = 0.4;
maximal_pwidth = 42;    % set maximal pulse width (ms) 66 for murata_140408
HM_offset_max = 650;   % set maximum HM elevation delay time (ms) ?]‚±‚Ì—pŒê‚Í‚ ‚Ü‚è“K?Ø‚Å‚Í‚È‚¢?B
HM_offset_min = -100;   % set maximum HM elevation delay time (ms)

fid = fopen([handles.fname '_ana.txt'],'w');

skip = false;
for ipause = 1 : handles.npause
    pstart = handles.pause(ipause, 1);
    pend = handles.pause(ipause, 2);
    
    % 0.4•bˆÈ‰º‚Ì’Z‚¢ŠÔŠu‚ÅŒÄ‹z’âŽ~‚ª‹N‚±‚Á‚½?ê?‡‚Í?A‚Ð‚Æ‚Â‚ÌŒÄ‹z’âŽ~‹æŠÔ‚Æ‚Ý‚È‚·?B140416•Ï?X
    if skip == true
        pause_duration = 0;
        skip = false;
    else
        if ipause < handles.npause
            if handles.pause(ipause+1, 1) <= round(pend + 0.4/handles.step)
                pend = handles.pause(ipause+1, 2);
                skip = true;
            end
        end
        pause_duration = handles.step * (pend - pstart); % sec 
    end
    
    % ŒÄ‹z’âŽ~‹æŠÔ’·‚ªminimal_pauseˆÈ?ã‚Å‚ ‚è?AŒÄ‹z’âŽ~‹æŠÔ“à‚Å‰¹‚ª”­?¶‚µ‚Ä‚¢‚é?ê?‡‚Ì‚ÝˆÈ‰º‚Ì?ˆ—?‚ð?s‚¤
    sound_cutoff_sel_index = get(handles.sound_cutoff, 'Value');
    switch sound_cutoff_sel_index
        case 1
            sound_criteria = (pause_duration >= minimal_pause) && ...
                (max(handles.xf(pstart:pend)) > handles.xth) && ...
                (max(handles.xf(pstart:pend)) < 1.2) && ...
                (sum(handles.xf(pstart:pend)) < 1.8 );   
            % 20140505 by Naomi.Yagi
%             sound_criteria =  (pause_duration >= minimal_pause) && ...
%                 (max(handles.xf(pstart:pend)) > handles.xth) && ...
%                 (max(handles.xf(pstart:pend)) < handles.xth*40 )  
            % modified on 20140429 by Naomi.Yagi
        case 2
            sound_criteria = (pause_duration >= minimal_pause) && (max(handles.xf(pstart:pend)) > handles.xth);
    end
    
    if sound_criteria
        iSound = iSound + 1;
        S_start(iSound) = pstart * handles.step;
%         S_duration(iSound) = pause_duration;
    
        % continuous wavelet transform using Morlet mother function.
        [TFR,timeVec,freqVec] = morletcwt(handles.SW1(pstart*handles.step*handles.Fs : pend*handles.step*handles.Fs)', freqVec, handles.Fs, 7); % debug 140409
        sTFR = sum(TFR);
%         TFRth = mean(sTFR) + 2 * std(sTFR);
        TFRth = max(sTFR) * 0.1;

%         a = sum(TFR, 2);
%         a = a ./ sum(a) .* 100;

        % ˜A‘±ƒEƒF?[ƒuƒŒƒbƒg•ÏŠ·‚Å‚Í?AŽã‚¢‰¹ˆ³‚Ì‰¹‚ª•\Ž¦‚³‚ê‚È‚¢‚½‚ß?A“r’†‚©‚ç‰¹‚Ì•\Ž¦‚É’ZŽžŠÔƒt?[ƒŠƒG•ÏŠ·‚ðŽg—p?B
        fftsize = 256;
        D_epoch2 = 0.1; % (sec)
        step2 = 0.01; % (sec)
        w2 = D_epoch2 * handles.Fs;
        noverlap2 = w2 - step2 * handles.Fs;
        [S, Fscale2, handles.timeScale2, P2] = ...
            spectrogram(handles.SW1(round(pstart*handles.step*handles.Fs) : round(pend*handles.step*handles.Fs)), hamming(w2), noverlap2, fftsize, handles.Fs);  % debug 140409

         % ƒ?ƒ‹ƒXƒyƒNƒgƒ‹‚É•ÏŠ·
         melFilterNum = 8;                      % ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚Ì•ªŠ„?”
         [mBank, bandpassFreq] = melFilterBank( Fscale2, handles.Fs, fftsize, melFilterNum );

         % —?‘z‰ž“š‚Ìƒ?ƒ‹ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚É?U•?ƒXƒyƒNƒgƒ‹‚ðŠ|‚¯?‡‚í‚¹‚Ä?AŠe‘Ñˆæ‚Ì
         % ƒXƒyƒNƒgƒ‹‚Ì˜a‚ð‹?‚ß?A?U•?ƒXƒyƒNƒgƒ‹‚ðmelFilterNumŽŸŒ³‚Éˆ³?k‚·‚é
         AdftSum = zeros(melFilterNum, length(handles.timeScale2));
         for i_epoch = 1 : length(handles.timeScale2);
             for count = 1 : melFilterNum
                 % ƒtƒBƒ‹ƒ^‚ð‚©‚¯‚é
                 AdftFilterBank = P2(:, i_epoch) .* mBank(count,:)';

                 % ˜a‚ð‚Æ‚é
                AdftSum(count, i_epoch) = sum(AdftFilterBank);
             end
         end
         a = sum(abs(AdftSum),2);
         a = a ./ sum(a) * 100;

        % ŒÄ‹z’âŽ~‹æŠÔ‚ÉŠÜ‚Ü‚ê‚é‰¹ƒpƒ‹ƒX‚ÌˆÊ’u?iP_pos?j?Aƒpƒ‹ƒX•??iP_width?j?AŽü”g?”?¬•ª?iPP?j?Aƒpƒ‹ƒX‘??”?inPulse?j‚ð‹?‚ß‚é
        iter2 = 0;
        iPulse = 0;
        first_pulse = 1;
        pulse_pos = 1;
        while iter2 < length(sTFR)
                iter2 = iter2 + 1;
                if sTFR(iter2) > TFRth
                    pulse_start = iter2;
                    while (sTFR(iter2) > TFRth) && (iter2 < length(sTFR))
                        iter2 = iter2 + 1;
                    end
                    pulse_end = iter2;
                    iPulse = iPulse + 1;
                    pulse_pos = iter2;
                    if iPulse == 1
                        first_pulse = pulse_pos;
                    end
                    P_pos(iSound, iPulse) = pulse_start;
                    P_width(iSound, iPulse) = (pulse_end - pulse_start) / handles.Fs * 1000; %(ms)
                    PP(iSound, iPulse, :) = sum(TFR(:, pulse_start:pulse_end-1), 2);
                end
        end
        P_count(iSound) = iPulse;
        last_pulse = pulse_pos;

%         a = squeeze(PP(iSound, 1, 5:end)); % debug 140425
%         for iPulse = 2 : P_count(iSound)
%             a = a + squeeze(PP(iSound, iPulse, 5:end));
%         end
%         a = a ./ sum(a) .* 100;
        
        istart = pstart - round(0.8 / handles.step); % ŒÄ‹z’âŽ~‚æ‚è0.8sec‘O
        if (istart <= 0) 
            istart = 1;
        end
        iend = pend + round(0.8 / handles.step); % ŒÄ‹z’âŽ~‚æ‚è0.8secŒã
        
        t2 = handles.t(round(istart*handles.step*handles.Fs) : round(iend*handles.step*handles.Fs))-handles.t(round(istart*handles.step*handles.Fs));  % debug 140409
        t2_breath = t2(1) : 1 / handles.Fs_breath : t2(end) + 1 / handles.Fs;    % for murata
        t2_hyoid = t2(1) : 1 / handles.Fs_hyoid : t2(end) + 1 / handles.Fs;     % for murata

        SW1_epoch = handles.SW11(round(istart*handles.step*handles.Fs_hyoid) : round(iend*handles.step*handles.Fs_hyoid));  % debug 140409
        SW1_epoch_hyoid = handles.SW2(round(istart*handles.step*handles.Fs_hyoid) : round(iend*handles.step*handles.Fs_hyoid));  % debug 140409
        SW_epoch = smooth(handles.SW2(round(istart*handles.step*handles.Fs_hyoid) : round(iend*handles.step*handles.Fs_hyoid))', 11);   % debug 140409
        flow_epoch = handles.flow(round(istart*handles.step*handles.Fs_breath) : round(iend*handles.step*handles.Fs_breath));  % debug 140409
        ISW2_epoch = detrend(handles.ISW2(round(istart*handles.step*handles.Fs_hyoid) : round(iend*handles.step*handles.Fs_hyoid)));   % debug 140409
        ISW2_epoch = ISW2_epoch / max(ISW2_epoch); 
        pause_start = handles.Fs_hyoid * 0.8; % debug 140425
        pause_end = length(SW_epoch) - handles.Fs_hyoid * 0.8; % debug 140425

        SWmax = max(SW_epoch(pause_start:pause_end));       % find local maxima within pause;
        SWmin = min(SW_epoch(pause_start:pause_end));       % find local minima within pause.
        HM_scale = max(handles.SW2) - min(handles.SW2);
        HM_amplitude = (SWmax - SWmin) / HM_scale * 100;    % (% max hyoid displacement)

        popup_sel_index = get(handles.popupmenu2, 'Value');
        switch popup_sel_index
            case 1
                criteria = (sum(a(3:5)) >= CWTpowerTh) && (HM_amplitude >= HMampTh) && ...
                (P_count(iSound) <= P_count_max) && (max(P_width(iSound,:)) <= maximal_pwidth);
                % ŒÄ‹z’âŽ~‹æŠÔ‚É‚¨‚¯‚é750HzˆÈ?ã‚Ì‰¹‚Ì?¬•ªƒp?[ƒZƒ“ƒg?isum(a(15:end))?j‚ªCWTpowerThˆÈ?ã?A‚©‚Â?A
                % ‰¹ƒpƒ‹ƒX?”‚ªP_count_maxˆÈ‰º?A‚©‚Â?A?Å‘åƒpƒ‹ƒX•?‚ªmaximal_pwidthˆÈ‰º?A‚©‚Â
                % ?ã?œ‚Ì•ÎˆÊ?U•??iHM_amplitude?j‚ªHMampThˆÈ?ã?A
                % ‚Å‚ ‚é?ê?‡‚É?A‚»‚ÌŒÄ‹z’âŽ~‹æŠÔ“à‚Åš‹‰º”½ŽË‚ª‹N‚«‚½‚Æ”»’è‚·‚é
            case 2
                criteria = true;
        end
        
        if criteria
            figure(handles.h); subplot(4,1,4);
            p = [handles.tt(istart) handles.tt(istart)];
            plot(p, ylim, 'r'); % š‹‰º?„’è‹æŠÔ‚ÌŽn‚Ü‚è‚ð?Ô‚Å•\Ž¦
            p = [handles.tt(iend) handles.tt(iend)];
            plot(p, ylim, 'g'); % š‹‰º?„’è‹æŠÔ‚Ì?I‚í‚è‚ð—Î‚Å•\Ž¦

            figure(handles.h2); subplot(2, 2, 1);
            hold off;
            plot(t2_hyoid, SW1_epoch, 'r');   % š‹‰º‰¹?i?¶ƒf?[ƒ^?j‚ð?Ô‚Å•\Ž¦
            hold on;
            plot(t2_hyoid, ISW2_epoch, 'g');  % ?ã?œ•ÏˆÊ?i?Ï•ª’l?j‚ð—Î‚Å•\Ž¦
            plot(t2_breath, flow_epoch, 'b');  % ŒÄ‹zƒtƒ??[‚ð?Â‚Å•\Ž¦
            xlim = [0 max(t2_hyoid)];
            axis([xlim ylim]);
            x = [t2_breath(pause_start) t2_breath(pause_start)];
            plot(x, ylim, 'y');         % ŒÄ‹z’âŽ~‹æŠÔ‚ÌŽn‚Ü‚è‚ð‰©‚Å•\Ž¦
            x = [t2_breath(pause_end) t2_breath(pause_end)];
            plot(x, ylim, 'y');         % ŒÄ‹z’âŽ~‹æŠÔ‚Ì?I‚í‚è‚ð‰©‚Å•\Ž¦
            
            [SWmax2, Indmax2] = max(sTFR);   % find sound local maxima within pause.
            Sound_II = pause_start+round(Indmax2/10)-1;

            %         ioffset = min(handles.Fs * 0.1, pause_start - 1);
            ioffset = min(handles.Fs_hyoid * 0.2, pause_start - 1);       % for murata: sampling frequency of breath = 1kHz
            % find HM local maxima within pause before Sound II. debug 140506
            [SWmax, Indmax] = findpeaks(SW_epoch(pause_start - ioffset : Sound_II + ioffset), 'SORTSTR', 'descend');   
%             [SWminX, IndminX] = findpeaks(-SW_epoch(pause_start - ioffset : pause_end));   % find HM local minima within pause.
            % ‚±‚±‚Å‚Í?A?ã?œ‚Ì“®‚«Žn‚ß‚©‚ç?Å‘åˆÊ‚É’B‚·‚é‚Ü‚Å‚ÌŽžŠÔ‚ðLEDT(laryngeal elevation delay time?j‚Æ’è‹`?A
            % ?¶‘Ì‰¹‹æŠÔ‚Ì’†‚Å?ã?œ‚ÌˆÚ“®‘¬“x‚ª?Å‘å’l‚ð‚Æ‚éˆÊ’u‚ðfindpeaksŠÖ?”‚Å‹?‚ß?A‘OŒü‚«?AŒã‚ëŒü‚«‚É?A
            % ‘¬“x‚ªƒ[ƒ?‚É‚È‚éˆÊ’u‚ð?A‚»‚ê‚¼‚êHM_end?AHM_start‚Æ‚µ‚Ä‹?‚ß?A‚»‚Ì?·‚ª45msˆÈ?ã‚Å‚ ‚ê‚ÎLEDT‚Æ‚·‚é?B
            LEDT = 0;
            ipeak = 0;
            HM_end = pause_end;
            HM_return = pause_end + handles.Fs_hyoid * 0.1;
            while ((HM_end >= pause_end) || (HM_end < pause_start) || (ISW2_epoch(HM_end) <= 0) || (LEDT < 45)) && ...
                    (ipeak < length(Indmax))      % debug 140927
                ipeak = ipeak + 1;
                Indmax1 = Indmax(ipeak);
                HM_slope_max = pause_start-ioffset+Indmax1-1;
                ipos = HM_slope_max;
                HM_Initial_Level = ISW2_epoch(HM_slope_max);
                while SW_epoch(ipos) > 0     % debug 140506
                    ipos = ipos - 1;
                    if ipos < 1
                        ipos = 1;
                        break;
                    end
                end
                HM_start = ipos;
                ipos = pause_start - ioffset + Indmax1 - 1;
    %           while SW_epoch(ipos) > 0
                while (SW_epoch(ipos) > 0) || (ISW2_epoch(ipos) <= -0.1)     % debug 140506
                    ipos = ipos + 1;
                    if ipos > length(SW_epoch)
                        ipos = length(SW_epoch);
                        break;
                    end
                end
                HM_end = ipos;
    %             LEDT = (HM_end - HM_start) / handles.Fs * 1000;        % (ms) 
                LEDT = (HM_end - HM_start) / handles.Fs_hyoid * 1000;     % (ms)       % for murata: sampling frequency of hyoid = 1kHz
                while (ISW2_epoch(ipos) > HM_Initial_Level)
%                 while (ISW2_epoch(ipos) > HM_Initial_Level) || (ipos <= Sound_II)    % debug 140512
                    ipos = ipos + 1;
                    if ipos > length(ISW2_epoch)
                        ipos = length(ISW2_epoch);
                        break;
                    end
                end
                HM_return = ipos;
            end % ((HM_end >= pause_end) || (ISW2_epoch(HM_end) <= 0) || (LEDT < 45) || ...
            
            x = [t2_hyoid(HM_start) t2_hyoid(HM_start)];
            plot(x, ylim, 'r');     % ?ã?œ?i?b?ó“î?œ?j‚Ì“®‚«Žn‚ß‚ð?Ô‚Å•\Ž¦
            x = [t2_hyoid(HM_slope_max) t2_hyoid(HM_slope_max)];
            h_slope_max = plot(x, ylim, 'g');     % ?ã?œ‚Ì“®‚«‚ª‹É‘å‚É‚È‚Á‚½Žž“_(?ã?œ‹}‘¬ˆÚ“®ŠJŽnŽž)‚ð—Î‚Å•\Ž¦ 140512‹L?q?C?³
            x = [t2_hyoid(HM_end) t2_hyoid(HM_end)];
            plot(x, ylim, 'b');     % ?ã?œ‚Ì?Å‘åˆÊ‚É’B‚·‚é?i?„’è?j“_‚ð?Â‚Å•\Ž¦
            x = [t2_hyoid(HM_return) t2_hyoid(HM_return)];
            h_return = plot(x, ylim, 'g');     % ?ã?œˆÊ’u‚ª?ã?œ‹}‘¬ˆÚ“®ŠJŽnŽž‚É–ß‚é?i?„’è?j“_‚ð—Î‚Å•\Ž¦ 140512‹L?q?C?³
            HM_offset = (Indmax2 / handles.Fs - (Indmax1 - ioffset) / handles.Fs_hyoid) * 1000; % debug 140409   
%             if (ipeak == length(Indmax)) || (ipos == length(ISW2_epoch)) || (HM_end < pause_start) || ...
%                     (HM_return >= pause_end + handles.Fs_hyoid * 0.4) || (ISW2_epoch(Sound_II) <= 0) % debug 140512
            if ((ipeak ~= 1) && (ipeak == length(Indmax))) || (ipos == length(ISW2_epoch))
                     
                title('unlikely swallow');
            elseif (HM_offset > HM_offset_max)
                title('delayed Sound II');
            elseif (HM_return >= pause_end + handles.Fs_hyoid * 0.4)
                title('delayed laryngeal relaxation')
            else
                title(['LEDT: ' num2str(LEDT) ' ms']);
            end

            figure(handles.h2); subplot(2,2,2);
%             surf(handles.timeScale2, Fscale2, 10*log10(abs(P2)),'EdgeColor','none');

            % ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚É‚æ‚Á‚Ä melFilterNum ŽŸŒ³‚Éˆ³?k‚³‚ê‚½‘Î?”?U•?ƒXƒyƒNƒgƒ‹‚ð‹?‚ß‚é
            % matlab note ‰¹?º‚Ì‰ð?Í?ihttp://shower.human.waseda.ac.jp/~m-kouki/pukiwiki_public/73.html?j 
            bandpassMedianFreq = median(bandpassFreq,2);    % ƒoƒ“ƒhƒpƒXƒtƒBƒ‹ƒ^‚Ì’†?SŽü”g?”
            surf(handles.timeScale2, bandpassMedianFreq, 10*log10(abs(AdftSum)),'EdgeColor','none');
            axis xy; axis tight; colormap(jet); view(0,90);
            xlabel('Time (sec)');
            ylabel('Frequency (Hz)');
            title(['mel-scale spectrogram %power 500-2300Hz: ' num2str(sum(a(3:5)), '%0.3g')]); % debug 140505

            figure(handles.h2); subplot(2, 2, 3); 
            hold off;
            plot(t2_hyoid, SW1_epoch_hyoid);       % for murata: sampling frequency of hyoid = 1kHz (20140313’²?®’†)  
            xlim = [0 max(t2)];
            axis([xlim ylim]);
            hold on;
            
            p = [t2_hyoid(Sound_II) t2_hyoid(Sound_II)]; 
            % sound local maxima within pause: normally squeezing sound during stage II transport.
            plot(p, ylim, 'y');     % ŒÄ‹z’âŽ~‹æŠÔ’†‚Ì?Å‘å‰¹”­?¶ˆÊ’u?i’Ê?í‚Ístage II transport’†‚Ìsqueezing‰¹?j‚ð‰©‚Å•\Ž¦ 
            
            p = [t2_hyoid(HM_slope_max) t2_hyoid(HM_slope_max)]; % HM local minima within pause        
            % (ms)  % for murata: sampling frequency of hyoid = 1kHz
            if (HM_offset > HM_offset_min) && (HM_offset < HM_offset_max)
            % Hyoid movement - sound offset is within HM_offset_max (ms).
            % ?ã?œ‚Ì“®‚«‚ª‹É‘å‚É‚È‚Á‚½Žž“_‚©‚ç?Å‘å‰¹”­?¶‚Ü‚Å‚ÌŽžŠÔ‚ªHM_offset_min(-100ms)?`HM_offset_max(650ms)‚Å‚ ‚ê‚Î?A
            % ?ã?œ‚Ì“®‚«‚ª‹É‘å‚É‚È‚Á‚½Žž“_‚ð?Ô‚Å?A‚»‚¤‚Å‚È‚¯‚ê‚Î?Â‚Å•\Ž¦‚·‚é?B’Ê?í?Aš‹‰º‚Ì?ê?‡?A?Ô‚É‚È‚é?B
                plot(p, ylim, 'r');
            else
                plot(p, ylim, 'b');
            end
            title(['Hyoid movement: ' num2str(HM_amplitude, '%0.3g') ' %max  offset: ' num2str(HM_offset) ' ms']);        
%             title(['Hyoid movement: ' num2str(HM_amplitude, '%0.3g') ' %max']); 
                       
            figure(handles.h2); subplot(2,2,4);
            hold off;
            plot(sTFR, 'b');    % ‰¹—Ê‚ð?Â‚Å•\Ž¦
            hold on;
            y = ones(1,length(sTFR)) * TFRth;
            plot(y, 'g');   % ‰¹ƒpƒ‹ƒX”»•Êè‡’l?iTFRth?j‚ð—Î‚Å•\Ž¦
            xlim = [0 length(sTFR)];
            axis([xlim ylim]);
            for iPulse = 1 : P_count(iSound)
                ipos = P_pos(iSound, iPulse);
                p = [ipos*handles.Fs ipos*handles.Fs];
                plot(p, ylim, 'r');     % ‰¹ƒpƒ‹ƒX‚ð?Ô‚Å•\Ž¦
            end                   
            title(['pulse count: ' num2str(P_count(iSound)) ' maximal pulse width: ' num2str(max(P_width(iSound,:)))]);

            % Figure 3 ‰æ–Ê‚Ì•\Ž¦
            % inspiration immediately before pause
            ibreath = 1;
            while (ibreath < handles.nbreath) && (handles.insp(ibreath) < pstart)
                ibreath = ibreath + 1;
            end
            
%             SW_param(iSound,6) = HM_start;
%             SW_param(iSound,7) = HM_slope_max;
%             SW_param(iSound,8) = HM_end;
%             SW_param(iSound,9) = Sound_II;
%             SW_param(iSound,10) = HM_return;
            
            oldphase = 0;
            insp_restart = 0;
            SWtype = 'unknown';
            SWtype2 = 'unknown';
            if (ibreath > 2)
               
                 % for murata: sampling frequency of breath = 1kHz
                if ibreath >= handles.nbreath - 1
                    last = min(length(handles.flow), length(handles.t_hyoid)); % debug 14/04/25
                    flow_epoch2 = handles.flow(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : last);
                    % SW_epoch2 = handles.SW1(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : int64(length(handles.SW1)/10) );
                    SW_epoch2 = handles.SW11(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : last);
                    ISW2_epoch2 = detrend(handles.ISW2(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : last));
                    t3 = handles.t_hyoid(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : last) ...
                        -handles.t_hyoid(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath)); % debug 14/04/30         
                else
                    flow_epoch2 = handles.flow(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : round(handles.insp(ibreath+1)*handles.step*handles.Fs_breath));
                    % SW_epoch2 = handles.SW1(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : round(handles.insp(ibreath+1)*handles.step*handles.Fs_breath));
                    SW_epoch2 = handles.SW11(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : round(handles.insp(ibreath+1)*handles.step*handles.Fs_breath));
                    ISW2_epoch2 = detrend(handles.ISW2(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : round(handles.insp(ibreath+1)*handles.step*handles.Fs_breath)));
                    t3 = handles.t_hyoid(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath) : ...
                        round(handles.insp(ibreath+1)*handles.step*handles.Fs_breath))-handles.t_hyoid(round(handles.exp(ibreath-2)*handles.step*handles.Fs_breath)); % debug 14/04/30                    
                end
                ISW2_epoch2 = ISW2_epoch2 / max(ISW2_epoch2); 
                
                figure(handles.h3); hold off; plot(t3, flow_epoch2, 'b');   % ŒÄ‹zƒtƒ??[‚ð?Â‚Å•\Ž¦
                hold all; plot(t3, SW_epoch2*5, 'r');                 % ?A“ª‰¹‚ð?Ô‚Å•\Ž¦
                plot(t3, ISW2_epoch2, 'g');   % ?A“ª•ÏˆÊ‚Ì?Ï•ª’l‚ð—Î‚Å•\Ž¦
                ind2 = round((handles.insp(ibreath-1)-handles.exp(ibreath-2))*handles.step*handles.Fs_breath);      % for murata: sampling frequency of breath = 1kHz
                while (flow_epoch2(ind2) <= 0) && (ind2 > 1)
                    ind2 = ind2 -1 ;
                end
                istart2 = ind2;
                ind2 = round((pstart-handles.exp(ibreath-2))*handles.step*handles.Fs_breath);                  % for murata: sampling frequency of breath = 1kHz
                                
%                 while (flow_epoch2(ind2) >= handles.isth) && (flow_epoch2(ind2) <= handles.esth)
%                     ind2 = ind2 - 1;
%                 end
                x = [t3(istart2) t3(istart2)];
                plot(x, ylim, 'r');         % ‹z‘§ŠJŽnˆÊ’u‚ð?Ô‚Å•\Ž¦
                
                iswallow_start = round((istart - handles.exp(ibreath-2)) * handles.step * handles.Fs_breath + HM_slope_max);

                % ŒÄ‹z’âŽ~‹æŠÔ‚ð“¯’è‚·‚é
                pause_duration = 0;
                while (ind2 < length(flow_epoch2)) && (pause_duration < minimal_pause) 
                    while (flow_epoch2(ind2) <= handles.isth) || (flow_epoch2(ind2) >= handles.esth)
                        ind2 = ind2 + 1;
                        if ind2 >= length(flow_epoch2)
                            break;
                        end
                    end
                    ipause_start = ind2;
                    ind2 = max(ind2, iswallow_start); % debug 140430
                    while (flow_epoch2(ind2) < handles.eth) && (flow_epoch2(ind2) > handles.ith)
                        ind2 = ind2 + 1;
                        if ind2 >= length(flow_epoch2)
                            break;
                        end
                    end
                    if ind2 < length(flow_epoch2)
                        while (flow_epoch2(ind2) <= handles.isth) || (flow_epoch2(ind2) >= handles.esth) 
                            ind2 = ind2 - 1;
                        end
                        ind2 = ind2 + 1;
                    end
                    ipause_end = min(ind2, length(flow_epoch2)); % debug 140428
                    pause_duration = (ipause_end - ipause_start) / handles.Fs_breath;      % for murata: sampling frequency of breath = 1kHz
                end
                
                % for murata: sampling frequency of breath = 1kHz
                while ((ind2 <= length(flow_epoch2) - handles.Fs_breath * handles.ti_min) && ...
                         ((flow_epoch2(ind2) > handles.isth) || (flow_epoch2(ind2 + handles.Fs_breath * handles.ti_min) > handles.ith)))          
                    ind2 = ind2 + 1;
                end
                irestart = min(ind2, length(flow_epoch2)); % debug 140428
                
                x = [t3(ipause_start) t3(ipause_start)];
                plot(x, ylim, 'y');     % ŒÄ‹z’âŽ~‹æŠÔ‚ÌŽn‚Ü‚è‚ð‰©‚Å•\Ž¦
                x = [t3(ipause_end) t3(ipause_end)];
                plot(x, ylim, 'y');     % ŒÄ‹z’âŽ~‹æŠÔ‚Ì?I‚í‚è‚ð‰©‚Å•\Ž¦
                
                x = [t3(iswallow_start) t3(iswallow_start)];
                plot(x, ylim, 'm');     % š‹‰º”½ŽËŠJŽn‚ðƒ}ƒ[ƒ“ƒ^‚Å•\Ž¦

                x = [t3(irestart) t3(irestart)];
                plot(x, ylim, 'g');     % š‹‰ºŒã‚Ì‹z‘§‚ÌŠJŽn‚ð—Î‚Å•\Ž¦

                % ŒÄ‹z’âŽ~‹æŠÔ’¼‘O100ms‚ÌŒÄ‹zƒtƒ??[‚É‚æ‚Á‚Äš‹‰º‘O‚ÌŒÄ‹z‚ª‹z‘§‚©ŒÄ‘§‚©‚ð”»’è‚·‚é
                if ipause_start - handles.Fs_breath * 0.2 < 0          % for murata: sampling frequency of breath = 1kHz
                    mflow = mean(flow_epoch2(1 : ipause_start));
                else
                    mflow = mean(flow_epoch2(ipause_start - handles.Fs_breath * 0.2 : ipause_start));     % for murata: sampling frequency of breath = 1kHz
                end
                if mflow > handles.esth
                    SWtype = 'E-SW';
                    SW_param(iSound,1) = 1;
                elseif mflow < handles.isth
                    SWtype = 'I-SW';
                    SW_param(iSound,1) = 2;
                else
                    SWtype = 'unknown';
                    SW_param(iSound,1) = 3;
                end
                
                % ŒÄ‹z’âŽ~‹æŠÔ’¼ŒãTi_min(300ms)‚ÌŒÄ‹zƒtƒ??[‚É‚æ‚Á‚Äš‹‰ºŒã‚ÌŒÄ‹z‚ª‹z‘§‚©ŒÄ‘§‚©‚ð”»’è‚·‚é
                if ipause_end == length(flow_epoch2)
                    SWtype2 = 'unknown';
                else
%                     if ipause_end + handles.Fs * handles.ti_min > length(flow_epoch2)
                    if ipause_end + handles.Fs_breath * handles.ti_min > length(flow_epoch2)       
                        mflow2 = mean(flow_epoch2(ipause_end : end));
                    else
%                         mflow2 = mean(flow_epoch2(ipause_end : ipause_end + handles.Fs * handles.ti_min));
                        mflow2 = mean(flow_epoch2(ipause_end : ipause_end + handles.Fs_breath * handles.ti_min));     % for murata: sampling frequency of breath = 1kHz
                        % taking negative pressure associated with bolus
                        % passage into account.
                    end
%                     if mflow2 > handles.esth
                    if mflow2 > 0
                        SWtype2 = 'SW-E';
                        SW_param(iSound,2) = 1;
                    elseif mflow2 < handles.isth
                        SWtype2 = 'SW-I';
                        SW_param(iSound,2) = 2;
                    else
                        SWtype2 = 'unknown';
                        SW_param(iSound,2) = 3;
                    end
                end
%                 oldphase = (ipause_start - istart2) / handles.Fs / handles.ts_mTT;
%                 insp_restart = (irestart - ipause_end) / handles.Fs;
                oldphase = (iswallow_start - istart2) / handles.Fs_breath / handles.ts_mTT;     % for murata: sampling frequency of breath = 1kHz
                insp_restart = (irestart - iswallow_start) / handles.Fs_breath;
                
                title(['SWtype: ' SWtype ' SWtype2: ' SWtype2 ' pause duration: ' num2str(pause_duration, '%0.3g') ...
                     ' oldphase: ' num2str(oldphase, '%0.3g') ' insp-restart: ' num2str(insp_restart)]);
                 
                SW_param(iSound,3) = pause_duration;
                SW_param(iSound,4) = oldphase;
                SW_param(iSound,5) = insp_restart;
           
            else
%             SWtype = 'unknown';
%             SWtype2 = 'unknown';
%             oldphase = 0;
%             insp_restart = 0;
        
            end % if ibreath > 2 ...

            button = 'Repeat';
            while strcmp(button, 'Repeat')
%                 sound(handles.SW1(istart*handles.step*handles.Fs : iend*handles.step*handles.Fs), handles.Fs);  % ?¶‘Ì‰¹‚ð?Ä?¶‚µ‚Äš‹‰º‰¹‚ðŠm”F debug 140409
                %sound(handles.SW1(pstart*handles.step*handles.Fs : pend*handles.step*handles.Fs), handles.Fs);  % ?¶‘Ì‰¹‚ð?Ä?¶‚µ‚Äš‹‰º‰¹‚ðŠm”F
%                 button = questdlg(['ipause = ', num2str(ipause), ' Is this swallowing?'],...
%                     ['Time = ', num2str(istart * handles.step), ' sec'],'Yes','No','Repeat','Yes');
                button = MFquestdlg([ 0.5 , 0.5 ], ['ipause = ', num2str(ipause), ' Is this swallowing?'],...
                    ['Time = ', num2str(istart * handles.step), ' sec'],'Yes','No','Repeat','Yes');
            end
            
            PRD_detection_sel_index = get(handles.PRD_detection, 'Value');
            switch button
                case 'Yes'
                    if PRD_detection_sel_index == 2
                        figure(handles.h2); subplot(2, 2, 1);
                        delete(h_slope_max);
                        delete(h_return);
                        w = 1;
                        while w == 1                                            
                            w = waitforbuttonpress; % wait until mouse button is pressed.
                        end
                        p = get(gca, 'CurrentPoint');
                        x = [p(1,1), p(1,1)];
                        plot(x, ylim, 'g');     % ?ã?œ‚Ì“®‚«‚ª‹É‘å‚É‚È‚Á‚½Žž“_(?ã?œ‹}‘¬ˆÚ“®ŠJŽnŽž)‚ð—Î‚Å•\Ž¦ 140512‹L?q?C?³
                        HM_slope_max = round((p(1,1) - t2(1)) * handles.Fs_hyoid);
                        w = 1;
                        while w == 1                                            
                            w = waitforbuttonpress; % wait until mouse button is pressed.
                        end
                        p = get(gca, 'CurrentPoint');
                        x = [p(1,1), p(1,1)];
                        plot(x, ylim, 'g');     % ?ã?œˆÊ’u‚ª?ã?œ‹}‘¬ˆÚ“®ŠJŽnŽž‚É–ß‚é?i?„’è?j“_‚ð—Î‚Å•\Ž¦ 140512‹L?q?C?³
                        HM_return = round((p(1,1) - t2(1)) * handles.Fs_hyoid);                    
                    end    
                    handles.nSwallow = handles.nSwallow + 1;
                    handles.SW_time(handles.nSwallow) = handles.tt(istart) + t2_hyoid(HM_slope_max);
                    SW_flag(iSound) = 1;
                    fprintf(fid,'%10d%5d%12.4f%12.4f%8s%8s%12.4f%12.4f%10d%10d%10d%10d%10d\n', ...
                          iSound, SW_flag(iSound), S_start(iSound), pause_duration,SWtype, SWtype2, oldphase, insp_restart, ...
                          HM_start, HM_slope_max, HM_end, Sound_II, HM_return);            
                    saveas(handles.h2, [handles.fname num2str(ipause) '.tif'],'tif') 
                case 'No'
                    SW_flag(iSound) = -1;
                case ''
                    break;
            end
      
%             fprintf(fid,'%10d%5d%12.4f%12.4f%12.4f%12.4f%10d%10d%10d%10d%10d\n', ...
%                 iSound, SW_flag(iSound), S_start(iSound), pause_duration, oldphase, insp_restart, ...
%                 HM_start, HM_slope_max, HM_end, Sound_II, HM_return);
            fprintf(fid,'%10d%5d%12.4f%12.4f%8s%8s%12.4f%12.4f%10d%10d%10d%10d%10d\n', ...
                  iSound, SW_flag(iSound), S_start(iSound), pause_duration,SWtype, SWtype2, oldphase, insp_restart, ...
                  HM_start, HM_slope_max, HM_end, Sound_II, HM_return);            
        end % if (sum(a(3:5)) >= CWTpowerTh) && ...
    end % if (pause_duration >= minimal_pause) && (max(handles.xf(pstart:pend)) > handles.xth) ...
end
fclose(fid);

% Œ‹‰Ê‚ðƒtƒ@ƒCƒ‹‚É•Û‘¶
% if ~strcmp(button, '')
%     save([fname '_ana.csv'], '-ASCII', '-TABS', 'iSound', 'SW_flag', 'S_start', 'SW_param');
%     save([fname '_ana2.csv'], '-ASCII', '-TABS', 'iSound', 'SW_flag', 'P_count', 'P_pos', 'P_width');
% end

% Update the handles on the interface
guidata(handles.output,handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spect_disp_sel_index = get(handles.spectrum_display, 'Value');
duration = 5; % (sec)
increment = 0.2;
p0_old = 0;
w = waitforbuttonpress;
while w == 0
    figure(handles.h); 
    p = get(gca, 'CurrentPoint');
    p0 = p(1,1);
    if p0 == p0_old
        break;
    else
       p0_old = p0;
    end
    x = [p0 p0]; 
    subplot(4,1,4);hold on;
    y1 = [1000 3000]; 
    z1 = [max(zlim), max(zlim)];
    
    switch spect_disp_sel_index
        case 1
            hl = plot3(x,y1,z1, 'r');
        case 2
            hl = plot(x,ylim, 'r');
    end
    
    ntimes = round(duration / increment);
    tic;
    sound(handles.SW1(round(p0*handles.Fs):round((p0+duration)*handles.Fs)), handles.Fs);
    for iter = 1 : ntimes - 1
        p0 = p0 + increment;
        x = [p0 p0];
        while toc < increment * iter
            % do nothing
        end
        delete(hl);
        switch spect_disp_sel_index
            case 1
                hl = plot3(x,y1,z1, 'r');
            case 2
                hl = plot(x,ylim, 'r');
        end
        drawnow;
    end
    delete(hl);        
    w = waitforbuttonpress;
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ŒÄ‹zƒtƒ??[ƒVƒOƒiƒ‹‚Ì•\Ž¦
figure(handles.h); subplot(4,1,2);
hold off; 
handles.tt = 0 : handles.step : (length(handles.timeScale)-1) * handles.step;
plot(handles.tt, handles.data); 
xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
%ylim = [-2 2];
axis([xlim ylim]);
xlabel('Time (sec)');
title('respiratory handles.flow signal');
hold all;

% ‚±‚±‚©‚çŒÄ‹zƒtƒ??[?M?†‚ð—p‚¢‚Ä?A‹z‹CŒÄ‹C‚Ìƒ^ƒCƒ~ƒ“ƒO‚ÆŒÄ‹zƒ|?[ƒY?i’âŽ~?j‚ÌŒŸ?o
handles.insp = zeros(500, 1);
handles.exp = zeros(500, 1);
handles.pause = zeros(500, 2);

handles.Fs2 = 1 / handles.step;

idx = 1;
ibreath = 0;
ipause = 0;

resp_detect_sel_index = get(handles.resp_detect, 'Value');
switch resp_detect_sel_index
    case 1
        handles.ith = str2double(get(handles.edit7, 'String'));    % -0.15
        handles.eth = str2double(get(handles.edit8, 'String'));      % 0.15
        handles.isth = str2double(get(handles.edit9, 'String'));  % -0.08
        handles.esth = str2double(get(handles.edit10, 'String'));  % 0.05; 
    case 2
        handles.ith = min(handles.data) * 0.1;
        handles.eth = max(handles.data) * 0.1;  % 0.2
        handles.isth = min(handles.data) * 0.05;
        handles.esth = max(handles.data) * 0.05;    % 0.1
end

handles.ti_min = 0.3; % minimum Ti in sec
handles.te_min = 0.5; % minimum Te in sec
handles.ts_min = str2double(get(handles.edit4, 'String')); % minimum Ts (pause duration for a possible swallowing)

figure(handles.h); subplot(4,1,2);
while idx < length(handles.data)
    if ibreath==0
        while handles.data(idx) > handles.eth
            idx = idx + 1;
        end
    else
        while (handles.data(idx) > handles.eth) && ((idx - iend) / handles.Fs2 > handles.te_min)
            idx = idx + 1;
        end
    end
    
    % detect E-I transition
    while handles.data(idx) > handles.ith
        if idx >= length(handles.data)
            break;
        end
        idx = idx + 1;
    end
    if idx < length(handles.data)
        while (handles.data(idx) < handles.isth) && (idx > 1)
            idx = idx - 1;
        end
        idx = idx + 1;
    else
        break;
    end       
    istart = idx - 1; % modified 13/06/26
        
    % detect I-E transition
    while handles.data(idx) <= handles.isth    % This could be data(idx) <= 0
        idx = idx + 1;
        if idx > length(handles.data)
                break;
        end
    end
%     if idx > length(handles.data)
%         break;
%     else
        iend = idx - 1; % modified 13/06/26
        ti = (iend - istart) / handles.Fs2;
        ipeak = min(handles.data(istart:iend));
%     end
    
    if (ti > handles.ti_min) && (ipeak < handles.ith)
        ibreath = ibreath + 1;
        handles.insp(ibreath) = istart;
        handles.exp(ibreath) = iend;
        y = [min(handles.data) max(handles.data)];

        % detect pause within a breath
        if ibreath > 1
            in_pause = false;
            ts = 0;
            for iter = handles.insp(ibreath - 1) + 1 : handles.insp(ibreath)
               if (handles.data(iter) > handles.isth) && (handles.data(iter) < handles.esth)    % pause
                   if in_pause
                       ts = ts + 1 / handles.Fs2;
                       
                       %----- for murata
                       if iter == handles.insp(ibreath)
                           pend = iter - 1; % modified 13/06/26
                           if ts >= handles.ts_min
                               ipause = ipause + 1;
                               handles.pause(ipause, 1) = pstart;
                               handles.pause(ipause, 2) = pend;
                               for pos = pstart : pend
                                   x = [pos pos] * handles.step;
                                   plot(x, y, 'y');
                               end
                           end
                       end 
                       %----
                       
                   else
                       in_pause = true;
                       pstart = iter - 1; % modified 13/06/26
                       ts = 1 / handles.Fs2;
                   end
               else
                   if in_pause
                       in_pause = false;
                       pend = iter - 1; % modified 13/06/26
                       if ts >= handles.ts_min
                           ipause = ipause + 1;
                           handles.pause(ipause, 1) = pstart;
                           handles.pause(ipause, 2) = pend;
                           for pos = pstart : pend
                               x = [pos pos] * handles.step;
                               plot(x, y, 'y');
                           end
                       end
                       ts = 0;
                   end
               end % if handles.pause ...
            end % for iter = handles.insp(ibreath - 1) + 1 : handles.insp(ibreath)
        end % if ibreath > 1 ...
        
        x = [istart istart] * handles.step;
        plot(x,y,'r');
        x = [iend iend] * handles.step;
        plot(x,y,'g');
        
    end % if (ti > handles.ti_min) && (ipeak < handles.ith) ...
end % while idx <= length(handles.data) ...
handles.nbreath = ibreath;
handles.npause = ipause;

% •½‹Ï‚ÌˆêŒÄ‹z‚É—v‚·‚éŽžŠÔ‚ð‹?‚ß‚é
handles.ts_mTT = (handles.insp(handles.nbreath) - handles.insp(1)) / (handles.nbreath - 1) * handles.step;
% This does not consider the time spent during swallowing.

% ˆêŒÄ‹z‚ð‚P‚Æ‚µ‚Ä?A•½‹Ï‚Ì‹z‹C?\ŒÄ‹C‚Ì?Ø‚è‘Ö‚í‚è‚Ìƒ^ƒCƒ~ƒ“ƒO‚ð‹?‚ß‚é
sTI = 0;
for ibreath = 1 : handles.nbreath - 1
    sTI = sTI + handles.exp(ibreath) - handles.insp(ibreath);
end
% mTI = sTI * handles.step / (handles.nbreath - 1);
% IE_transition = mTI / handles.ts_mTT;

plot(handles.tt, handles.data, 'b');
x = [0 (length(handles.timeScale)-1) * handles.step];
y = [0 0];
plot(x, y, 'k');

handles.tt = 0 : handles.step : (length(handles.timeScale)-1) * handles.step;
figure(handles.h); subplot(4,1,4);
plot(handles.tt, handles.xf); xlim = [0 (length(handles.timeScale)-1) * handles.step]; axis([xlim ylim]);

% Update the handles on the interface
guidata(handles.output,handles);


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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

Tb_sel_index = get(handles.togglebutton1, 'Value');
zoom_size = str2double(get(handles.edit5, 'String'));

figure(handles.h);

switch Tb_sel_index
    case 0
        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,1); 
        xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
        axis([xlim ylim]);

        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,2); 
        xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
        axis([xlim ylim]);

        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,3); 
        xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
        axis([xlim ylim]);

        subplot(4,1,4);
        hold off; 
        plot(handles.tt, handles.xf); 
        xlim = [0 (length(handles.timeScale)-1) * handles.step]; 
        axis([xlim ylim]);
        hold on;
        for iter = 1 : handles.nSwallow
            x = [handles.SW_time(iter) handles.SW_time(iter)];
            plot(x,ylim,'m');
        end
        
    case 1
        p = get(gca, 'CurrentPoint');
        p0 = round(p(1,1));
        if p0 <= 0
            p0 = 1;
        end
        
        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,1); 
        xlim = [p0 p0+zoom_size]; 
        axis([xlim ylim]);

        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,2); 
        xlim = [p0 p0+zoom_size]; 
        axis([xlim ylim]);

        % •\Ž¦‚Ì‚w?À•WŽ²‚ð‚ ‚í‚¹‚é
        subplot(4,1,3); 
        xlim = [p0 p0+zoom_size]; 
        axis([xlim ylim]);

        spect_disp_sel_index = get(handles.spectrum_display, 'Value');
        switch spect_disp_sel_index
            case 2
                subplot(4,1,4);
                hold off; 
                plot(handles.tt, handles.xf); 
                xlim = [p0 p0+zoom_size]; 
                axis([xlim ylim]);
                hold on;
                for iter = 1 : handles.nSwallow
                     if (handles.SW_time(iter) > min(xlim)) && (handles.SW_time(iter) < max(xlim))
                         x = [handles.SW_time(iter) handles.SW_time(iter)];
                         plot(x,ylim,'m');
                     end
                end
            case 1
                 % ’ZŽžŠÔƒt?[ƒŠƒG•ÏŠ·
                 step = 0.2; % (sec)
                 fftsize = 256;
                 D_epoch = 1.5; % (sec)
                 w = D_epoch * handles.Fs;
                 noverlap = w - step * handles.Fs;
                 if (p0+zoom_size)*handles.Fs >= length(handles.SW1)
                  [S, fscale, timeScale, P] = ...
                     spectrogram(handles.SW1(p0*handles.Fs : end), hamming(w), noverlap, fftsize, handles.Fs);
                 else   
                  [S, fscale, timeScale, P] = ...
                     spectrogram(handles.SW1(p0*handles.Fs : (p0+zoom_size)*handles.Fs), hamming(w), noverlap, fftsize, handles.Fs);
                 end

                 % ƒ?ƒ‹ƒXƒyƒNƒgƒ‹‚É•ÏŠ·
                 melFilterNum = 8;                      % ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚Ì•ªŠ„?”
                 [mBank, bandpassFreq] = melFilterBank( fscale, handles.Fs, fftsize, melFilterNum );

                 % —?‘z‰ž“š‚Ìƒ?ƒ‹ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚É?U•?ƒXƒyƒNƒgƒ‹‚ðŠ|‚¯?‡‚í‚¹‚Ä?AŠe‘Ñˆæ‚Ì
                 % ƒXƒyƒNƒgƒ‹‚Ì˜a‚ð‹?‚ß?A?U•?ƒXƒyƒNƒgƒ‹‚ðmelFilterNumŽŸŒ³‚Éˆ³?k‚·‚é
                 AdftSum = zeros(melFilterNum, length(timeScale));
                 for i_epoch = 1 : length(timeScale);
                     for count = 1 : melFilterNum
                         % ƒtƒBƒ‹ƒ^‚ð‚©‚¯‚é
                         AdftFilterBank = P(:, i_epoch) .* mBank(count,:)';

                         % ˜a‚ð‚Æ‚é
                        AdftSum(count, i_epoch) = sum(AdftFilterBank);
                     end
                 end

                % ƒtƒBƒ‹ƒ^ƒoƒ“ƒN‚É‚æ‚Á‚Ä melFilterNum ŽŸŒ³‚Éˆ³?k‚³‚ê‚½‘Î?”?U•?ƒXƒyƒNƒgƒ‹‚ð‹?‚ß‚é
                % matlab note ‰¹?º‚Ì‰ð?Í?ihttp://shower.human.waseda.ac.jp/~m-kouki/pukiwiki_public/73.html?j 
                 bandpassMedianFreq = median(bandpassFreq,2);    % ƒoƒ“ƒhƒpƒXƒtƒBƒ‹ƒ^‚Ì’†?SŽü”g?”
                %  h1 = figure;
                 figure(handles.h); subplot(4,1,4); 
                 hold off;
                 %xlim = [p0 p0+zoom_size]; 
                 surf(timeScale+p0, bandpassMedianFreq, 10*log10(abs(AdftSum)),'EdgeColor','none');
                 axis xy; axis tight; colormap(jet); view(0,90);
                 xlabel('Time (sec)');
                 ylabel('Frequency (Hz)');
                 title('mel-scale spectrogram');
                 hold on;
                 y1 = [1000 3000]; 
                 z1 = [max(zlim), max(zlim)];
                 for iter = 1 : handles.nSwallow
                     if (handles.SW_time(iter) > min(xlim)) && (handles.SW_time(iter) < max(xlim))
                         x = [handles.SW_time(iter) handles.SW_time(iter)];
                         plot3(x,y1,z1, 'm');
                     end
                 end
        end
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zero_level_Callback(hObject, eventdata, handles)
% hObject    handle to zero_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zero_level as text
%        str2double(get(hObject,'String')) returns contents of zero_level as a double


% --- Executes during object creation, after setting all properties.
function zero_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in sound_cutoff.
function sound_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to sound_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sound_cutoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sound_cutoff


% --- Executes during object creation, after setting all properties.
function sound_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sound_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in resp_detect.
function resp_detect_Callback(hObject, eventdata, handles)
% hObject    handle to resp_detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns resp_detect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from resp_detect


% --- Executes during object creation, after setting all properties.
function resp_detect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resp_detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in zero_flow.
function zero_flow_Callback(hObject, eventdata, handles)
% hObject    handle to zero_flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns zero_flow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zero_flow


% --- Executes during object creation, after setting all properties.
function zero_flow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_flow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in highpass_filter.
function highpass_filter_Callback(hObject, eventdata, handles)
% hObject    handle to highpass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns highpass_filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from highpass_filter


% --- Executes during object creation, after setting all properties.
function highpass_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highpass_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PRD_detection.
function PRD_detection_Callback(hObject, eventdata, handles)
% hObject    handle to PRD_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PRD_detection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PRD_detection


% --- Executes during object creation, after setting all properties.
function PRD_detection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PRD_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trigger_signal.
function trigger_signal_Callback(hObject, eventdata, handles)
% hObject    handle to trigger_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns trigger_signal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trigger_signal


% --- Executes during object creation, after setting all properties.
function trigger_signal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trigger_signal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in spectrum_display.
function spectrum_display_Callback(hObject, eventdata, handles)
% hObject    handle to spectrum_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns spectrum_display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from spectrum_display


% --- Executes during object creation, after setting all properties.
function spectrum_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spectrum_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


