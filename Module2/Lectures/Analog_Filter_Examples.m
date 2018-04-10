clear; close('all');

%% Filter specifications (lowpass)
% Passband/stopband frequencies in Hz
fp = 4; fs = 4.5;
Wp = 2*pi*fp; Ws = 2*pi*fs;
% Passband/stopband attenuation
Gp = 0.95; Gs = 0.05;
Rp = -20*log10(Gp);
Rs = -20*log10(Gs);
ep = sqrt(1/(Gp*Gp) - 1); es = sqrt(1/(Gs*Gs) - 1); 

%% Butterworth
% Compute filter order and 3dB cutoff frequency
[N,Wc] = buttord(Wp,Ws,Rp,Rs,'s');
[B,A] = butter(N,Wc,'s');
W = linspace(0,2*pi*10,1000);
% Evaluate magnitude response
H = freqs(B,A,W);
% Approximate group delay by differencing unwrapped phase response
grpdelay = -diff(unwrap(angle(H)))./diff(W);
figure(1)
subplot(321)
plot(W/(2*pi),abs(H),'.')
xlabel('Frequency (Hz)'); ylabel('Magnitude Response');
title_string = ['Butterworth, ',sprintf('N=%2d',N)];
title(title_string);
axis([0,10,0,1.1]); modify_figure; grid on;
subplot(322)
plot(W(2:end)/(2*pi),grpdelay,'.')
xlabel('Frequency (Hz)'); ylabel('Group Delay (s)');
axis([0,10,0,2.5]); modify_figure; grid on;

%% Chebyshev Type I
% Compute filter order and revised passband frequency
[N,Wpnew] = cheb1ord(Wp,Ws,Rp,Rs,'s');
[B,A] = cheby1(N,Rp,Wpnew,'s');
H = freqs(B,A,W);
grpdelay = -diff(unwrap(angle(H)))./diff(W);
figure(1)
subplot(323)
plot(W/(2*pi),abs(H),'.')
xlabel('Frequency (Hz)'); ylabel('Magnitude Response');
title_string = ['Chebyshev Type I, ',sprintf('N=%2d',N)];
title(title_string);
axis([0,10,0,1.1]); modify_figure; grid on;
subplot(324)
plot(W(2:end)/(2*pi),grpdelay,'.')
xlabel('Frequency (Hz)'); ylabel('Group Delay (s)');
axis([0,10,0,2.5]); modify_figure; grid on;

%% Chebyshev Type II
% Compute filter order and revised stopband frequency
[N,Wsnew] = cheb2ord(Wp,Ws,Rp,Rs,'s');
[B,A] = cheby2(N,Rs,Wsnew,'s');
H = freqs(B,A,W);
grpdelay = -diff(unwrap(angle(H)))./diff(W);
figure(1)
subplot(325)
plot(W/(2*pi),abs(H),'.')
xlabel('Frequency (Hz)'); ylabel('Magnitude Response');
title_string = ['Chebyshev Type II, ',sprintf('N=%2d',N)];
title(title_string);
axis([0,10,0,1.1]); modify_figure; grid on;
subplot(326)
plot(W(2:end)/(2*pi),grpdelay,'.')
xlabel('Frequency (Hz)'); ylabel('Group Delay (s)');
axis([0,10,0,2.5]); modify_figure; grid on;

%% Elliptic
% Compute filter order and revised passband frequency
[N,Wpnew] = ellipord(Wp,Ws,Rp,Rs,'s');
[B,A] = ellip(N,Rp,Rs,Wpnew,'s');
H = freqs(B,A,W);
grpdelay = -diff(unwrap(angle(H)))./diff(W);
figure(2)
subplot(321)
plot(W/(2*pi),abs(H),'.')
xlabel('Frequency (Hz)'); ylabel('Magnitude Response');
title_string = ['Elliptic, ',sprintf('N=%2d',N)];
title(title_string);
axis([0,10,0,1.1]); modify_figure; grid on;
subplot(322)
plot(W(2:end)/(2*pi),grpdelay,'.')
xlabel('Frequency (Hz)'); ylabel('Group Delay (s)');
axis([0,10,0,2.5]); modify_figure; grid on;

%% Bessel
% Filter order chosen arbitrarily, no control over magnitude response
N = 25;
[B,A] = besself(N,Wp);
H = freqs(B,A,W);
grpdelay = -diff(unwrap(angle(H)))./diff(W);
figure(2)
subplot(323)
plot(W/(2*pi),abs(H),'.')
xlabel('Frequency (Hz)'); ylabel('Magnitude Response');
title_string = ['Bessel, ',sprintf('N=%2d',N)];
title(title_string);
axis([0,10,0,1.1]); modify_figure; grid on;
subplot(324)
plot(W(2:end)/(2*pi),grpdelay,'.')
xlabel('Frequency (Hz)'); ylabel('Group Delay (s)');
axis([0,10,0,2.5]); modify_figure; grid on;
