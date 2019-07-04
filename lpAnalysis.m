function [ lpc, lpres, pitch, formants ] = lpAnalysis( sig, p, fs )

% Estimate LP Coefficients of the given signal
% using Auto-Correlation Method.
% 
% INPUT:
% sig     =   signal upon which lp analysis is to be performed
% p   =   lp order
% fs  =   sampling frequency of the signal (preferably 8kHz)
% 
% OUTPUT:
% lpc     =   LP Coeffiecients
% lpres   =   LP Residual
% pitch   =   Pitch value in milli-seconds
% formants    =   Formants of the given segment
%
%
% Example:
% [s,fs]  =   wavread('Example.wav'); s=s(:,1);
% sig =   s(78501:80000);
% sig   =   resample(sig,8000,fs);
% fs    =   8000;
% p=12;
% [ lpc, lpres, pitch, formants ] = lpAnalysis( sig, p, fs )

sig =   sig./max(abs(sig)); % Normalize the signal
time    =   (1:length(sig))*1000/fs;

win =   window(@hamming,length(sig)); % Get hamming window equal to size of signal
winsig =   sig.*win; % Windowed signal
% winsig  =   sig;

acorrwinsig     =   xcorr(winsig); % Compute autocorrelation of windowed signal
acorrwinsig     =   acorrwinsig./max(abs(acorrwinsig));
acorrwinsig     =   acorrwinsig(length(sig):end);

r   =   acorrwinsig(2:p+1); % pth order autocorrelation sequence
R   =   toeplitz(acorrwinsig(1:p)); % Toeplitz matrix
lpc =   -inv(R)*r; % LP Coefficients

lpres   =   filter([1; lpc],1,winsig);  % LP Residual
lpres   =   lpres./max(abs(lpres));

acorrlpres  =   xcorr(lpres);   % Autocorrelation of LP Residual
acorrlpres  =   acorrlpres./max(abs(acorrlpres));
acorrlpres  =   acorrlpres(length(lpres):end);

% Pitch Extraction from LP Residual
[PKS,LOCS]  =   findpeaks(acorrlpres);
[~, I]  =   max(PKS);
pitch   =   LOCS(I)-1;
pitch   =   pitch*1000/fs;  % in milliseconds

% Formant Extraction from LP Residual
[mag,freq] = freqz(1,[1; lpc],256,fs);
mag     =   abs(mag);
[~, I]  =   findpeaks(mag);
nf  =   floor((p-2)/2); % no. of formants
if length(I)<nf
    formants     =   freq(I);
else
    formants    =   freq(I(1:nf)); % formants
end

error=[];
for p=1:100
    
    r1   =   acorrwinsig(2:p+1); % pth order autocorrelation sequence
    R1   =   toeplitz(acorrwinsig(1:p)); % Toeplitz matrix
    lpc1 =   -inv(R1)*r1; % LP Coefficients
    
    lpres1   =   filter([1; lpc1],1,winsig);  % LP Residual
    lpres1   =   lpres1./max(abs(lpres1));
    temp    =   xcorr(lpres1);
%     error(end+1)    =   acorrwinsig(1)+(lpc1'*acorrwinsig(1:p));
error(end+1)    =   temp(1)/acorrwinsig(1);
end

% Plots
h   =   figure;
ax(1)   =   subplot(311); plot(time, sig, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Signal'); grid on;
ax(2)   =   subplot(312); plot(time, winsig, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Windowed Signal'); grid on;
ax(3)   =   subplot(313); plot(time, acorrwinsig, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Auto-Correlation of the Windowed Signal'); grid on;
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
axis tight;

h   =   figure;
ax(1)   =   subplot(311); plot(time, sig, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Signal'); grid on;
ax(2)   =   subplot(312); plot(time, lpres, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('LP Residual of the given Signal'); grid on;
ax(3)   =   subplot(313); plot(time, acorrlpres, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Auto-Correlation of the LP Residual'); grid on;
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
axis tight; 

h   =   figure;
subplot(211); plot(time, sig, 'LineWidth',2.0); xlabel('time(ms)'); ylabel('Amplitue'); title('Signal'); grid on;
subplot(212); plot(freq,mag, 'LineWidth',2.0); xlabel('Freqeuncy(hz)'); ylabel('Magnitude(dB)'); title('Magnitude Spectrum of LP Filter'); grid on;
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);

h   =   figure;
plot(error,'LineWidth', 2.0); xlabel('LP Order'); ylabel('Normalized Error'); title('Normalized Error'); grid on;
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
end

