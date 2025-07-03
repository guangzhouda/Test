%% realtime_voice_beautifier.m
% -------------------------------------------------------------
% 实时低复杂度人声美化 + 谱减降噪 + A/B 比较 Demo
% -------------------------------------------------------------
clear; clc;

%% 用户可调参数 -----------------------------------------------
fs          = 16000;   % 采样率
frameLen    = 256;     % 帧长 (≈16 ms)
tgtRMS      = 0.20;    % AGC 目标 RMS
lowGain     = 6;       % 低频搁架 EQ @120 Hz
midGain     = 3;       % 中频峰值 EQ @1 kHz
highGain    = -4;      % 高频搁架 EQ @6 kHz
delay_ms    = 20;      % 混响延迟
decay       = 0.30;    % 混响衰减
mode        = 'x';   % raw | beautified | ab
%% 降噪参数
alphaNoise  = 0.95;    % 噪声谱平滑因子
gainFloor   = 0.05;    % 最低增益
win         = hann(frameLen,'periodic');
fftLen      = frameLen;
noisePSD    = zeros(fftLen,1);
% -------------------------------------------------------------

%% 实时设备
mic = audioDeviceReader('SampleRate',fs,'SamplesPerFrame',frameLen);
spk = audioDeviceWriter('SampleRate',fs);

%% 预设计三段 EQ
[bL,aL] = shelfEQ('low',  fs, 120,   lowGain);
[bM,aM] = peakEQ(         fs,1000,1, midGain);
[bH,aH] = shelfEQ('high', fs,6000, highGain);

%% 简易混响
delaySamp = round(fs*delay_ms/1000);
delayBuf  = zeros(delaySamp,1);
idx       = 1;

fprintf('实时人声美化已启动，按 Ctrl+C 结束...\n');
%% 主循环
while true
    %% 读麦克风
    x = mic();
    x = mean(x,2);          % 立体声转单声
    
    %% 0) 谱减降噪
    X    = fft(win.*x,fftLen);
    Px   = abs(X).^2;
    noisePSD = alphaNoise*noisePSD + (1-alphaNoise)*Px;
    G    = max(gainFloor,1-noisePSD./(Px+eps));
    x_dn = real(ifft(X.*G));
    
    %% 1) AGC
    g  = tgtRMS/(rms(x_dn)+eps);
    x1 = x_dn * g;
    
    %% 2) 三段 EQ
    y = filter(bL,aL,x1);
    y = filter(bM,aM,y);
    y = filter(bH,aH,y);
    
    %% 3) 轻混响
    yOut = zeros(size(y));
    for n=1:frameLen
        yOut(n)       = y(n) + decay*delayBuf(idx);
        delayBuf(idx) = y(n);
        idx = idx + 1; if idx>delaySamp, idx=1; end
    end
    
    %% 4) 监听模式选择
    switch lower(mode)
        case 'raw'
            outFrame = x;                  % 仅原始
        case 'beautified'
            outFrame = yOut;               % 仅美化
        case 'ab'
            outFrame = [x yOut];           % 左原始右美化
        otherwise
            outFrame = yOut;               % 默认美化
    end
    
    %% 5) 输出
    spk(outFrame);
end
%% -------------------------------------------------------------
%% 辅助 EQ 设计函数
function [b,a] = shelfEQ(type,fs,fc,gain_dB)
    A   = 10^(gain_dB/40); w0=2*pi*fc/fs; cosw=cos(w0);
    alpha = sin(w0)/sqrt(2);
    if strcmpi(type,'low')
        b=[A*((A+1)-(A-1)*cosw+2*sqrt(A)*alpha),...
           2*A*((A-1)-(A+1)*cosw),...
           A*((A+1)-(A-1)*cosw-2*sqrt(A)*alpha)];
        a=[(A+1)+(A-1)*cosw+2*sqrt(A)*alpha,...
          -2*((A-1)+(A+1)*cosw),...
           (A+1)+(A-1)*cosw-2*sqrt(A)*alpha];
    else
        b=[A*((A+1)+(A-1)*cosw+2*sqrt(A)*alpha),...
          -2*A*((A-1)+(A+1)*cosw),...
           A*((A+1)+(A-1)*cosw-2*sqrt(A)*alpha)];
        a=[(A+1)-(A-1)*cosw+2*sqrt(A)*alpha,...
            2*((A-1)-(A+1)*cosw),...
           (A+1)-(A-1)*cosw-2*sqrt(A)*alpha];
    end
    b=b/a(1); a=a/a(1);
end

function [b,a] = peakEQ(fs,fc,Q,gain_dB)
    A=10^(gain_dB/40); w0=2*pi*fc/fs; alpha=sin(w0)/(2*Q); cosw=cos(w0);
    b=[1+alpha*A,-2*cosw,1-alpha*A];
    a=[1+alpha/A,-2*cosw,1-alpha/A];
    b=b/a(1); a=a/a(1);
end
