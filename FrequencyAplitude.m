%
% Read in and label data
%

[Raw_data,fs] = audioread('/path/to/wav'); % This function provides amplitude data from however many audio channels recorded in the file.
% fs = sampling frequency of the recorder using to record the audio file.


x = Raw_data(:,1); % Amplitude data for channel 1. 
Amp_ch2 = Raw_data(:,2); % Amplitude data for channel 2.
T = 1/fs; % sampling period.
L = numel(x);


%
% Create time vector
%

t = zeros(L,1);

for n = 2:numel(x)
    t(n) = t(n-1) + T;
end

%
% Plot time series data
%

figure
yyaxis left
plot(t,x,'b');
yyaxis right
plot(t,Amp_ch2,'r');
yyaxis left
title('Audio recording')
xlabel('Time / seconds')
ylabel('Channel 1 / arbitrary units')
yyaxis right
ylabel('Channel 2 / arbitrary units')
set(gcf,'color','w')
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
xlim([0 25]);

%
% Fast Fourier Transform
%

fft_1 = fft(x);

P2 = abs(fft_1/L); % Complex magnitude of fft normalised by L
P1 = P2(1:L/2+1); % Reed up on the sampling - i.e. mx f = fs/2?
P1(2:end-1) = 2*P1(2:end-1);

%
% Plot fft
%

figure
f = fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 7000]);


windowSize = 100; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y = filter(b,a,P1);

%B = 1/10*ones(100,1);
%outa = filter(B,1,P1);

hold on
plot(f,y);



winLength = round(0.05*fs);
overlapLength = round(0.045*fs);
[f0,idx] = pitch(x,fs,'Method','SRH','WindowLength',winLength,'OverlapLength',overlapLength);
tf0 = idx/fs;

%sound(x,fs) %plays the audio - can comment out 

figure
tiledlayout(2,1)

nexttile
plot(t,x)
ylabel('Amplitude')
title('Audio Signal')
axis tight

nexttile
plot(tf0,f0)
xlabel('Time (s)')
ylabel('Pitch (Hz)')
title('Pitch Estimations')
axis tight
