%% Problem 1
wp = 0.6*pi;
ws = 0.5*pi;
Rs = 40;    %dB

% Compute kaiser window filter order, wn, and beta parameters
[N,wn,beta,ftype] = kaiserord( [ws/pi wp/pi], [0 1], [10^(-Rs/20) 10^(-Rs/20)] );

% Create filter using kaiser window
b = fir1(N, wn, ftype, kaiser(N+1), 'noscale');

% Impulse Response
figure(1)
stem((0:N),b);
title('Filter Impulse Response'); xlabel('Samples (n)'); ylabel('Amplitude');

% Frequency Response
figure(2)
[H,w] = freqz(b);
plot(w/pi,20*log10(abs(H)));
title('Filter Frequency Response');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude (dB)');
line([0 wp/pi],[-Rs -Rs],'color','red','LineStyle','--');
line([wp/pi 1],[0.01 0.01],'color','red','LineStyle','--');
line([wp/pi wp/pi],[-Rs 0],'color','red','LineStyle','--');
legend('Kaiser Window LPF', 'Ideal LPF');

% Zoomed view of stopband
figure(3)
plot(w/pi,20*log10(abs(H)));
title('Zoomed View of Stopband');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude (dB)');
axis([0.4 0.55 -Rs-5 -Rs+5]);
line([0 wp/pi],[-Rs -Rs],'color','red','LineStyle','--');
line([wp/pi 1],[0.01 0.01],'color','red','LineStyle','--');
line([wp/pi wp/pi],[-Rs 0],'color','red','LineStyle','--');
legend('Kaiser Window LPF', 'Ideal LPF');

% Zoomed view of passband
figure(4)
plot(w/pi,20*log10(abs(H)));
title('Zoomed View of Passband');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude (dB)');
axis([0.55 0.7 -20 1]);
line([0 wp/pi],[-Rs -Rs],'color','red','LineStyle','--');
line([wp/pi 1],[0.01 0.01],'color','red','LineStyle','--');
line([wp/pi wp/pi],[-Rs 0],'color','red','LineStyle','--');
legend('Kaiser Window LPF', 'Ideal LPF');


%% Problem 2
clear;

wp1 = 0.4*pi;%0.09;
wp2 = 0.6*pi;%0.09;
ws1 = 0.3*pi;
ws2 = 0.7*pi;
Rp  =  4;     %dB
Rs  = 40;     %dB


% Create the Dolph_Chebyshev window using N+1 since the fir1() function
% expects an N+1 window. In this case the filter order is 30. 
N = 30;


% Adjusting initial parameters so filter meets spec @ minimum order
wp1 = wp1 - 0.05;           % Increasing filter bandwidth
wp2 = wp2 + 0.05;           % Increasing filter bandwidth
sLobeLevel  = Rs - 10;      % Adjusting window sidelobe level


% Create Dolph-Chebyshev window
w = chebwin(N+1,sLobeLevel);


% Plotting impulse response of window
figure(5)
stem((0:length(w)-1),w);
title('Dolph-Chebyshev Impulse Response');
xlabel('Samples(n)'); ylabel('Amplitude');


% Plotting frequency response of window
figure(6)
freqz(w);
title('Dolph-Chebyshev Frequency response');


% Using fir1() to create bandpass filter using a Dolph-Chebyshev window.
% The window size is N+1 since fir1() expects an N+1 window.
h = fir1(N,[wp1/pi wp2/pi],'bandpass',chebwin(N+1,sLobeLevel));


% Compute frequency response of the filter.
[H,wn] = freqz(h);
H_dB = 20*log10(abs(H));


% Plot Magnitude Response of BPF
figure(7)
hold on
plot(wn/pi,H_dB);
title('BP-Filter Frequency Response');
xlabel('Normalized Frequency (x pi rad/samp)'); ylabel('Magnitude (dB)');
grid('on');
line([0.4 0.4],[-40 0],'color','red','LineStyle','--');
line([0.6 0.6],[-40 0],'color','red','LineStyle','--');
line([0.3 0.3],[-40 -4],'color','red','LineStyle','--');
line([0.7 0.7],[-40 -4],'color','red','LineStyle','--');
line([0.3 0.7],[-4 -4],'color','red','LineStyle','--');
line([0 0.3],[-40 -40],'color','red','LineStyle','--');
line([0.7 1],[-40 -40],'color','red','LineStyle','--');
hold off


% Zoomed in view of passband
figure(8)
hold on
plot(wn/pi,H_dB);
plot(0.5,0,'x');
title('Zoomed View of BPF Passband');
xlabel('Normalized Frequency (x pi rad/samp)'); ylabel('Magnitude (dB)');
axis([0.39 0.61 -4.5 0.1]);
grid('on');
text(0.5,0.1,num2str(abs(H_dB(255))));
line([0.4 0.4],[-100 0],'color','red','LineStyle','--');
line([0.6 0.6],[-100 0],'color','red','LineStyle','--');
line([0.3 0.7],[-4 -4],'color','red','LineStyle','--');
hold off


% Zoomed in view of stopband
figure(9)
hold on
plot(wn/pi,H_dB);
plot(0.5,0,'x');
title('Zoomed View of BPF Stopband');
xlabel('Normalized Frequency (x pi rad/samp)'); ylabel('Magnitude (dB)');
axis([0.2 0.8 -45 -39.5]);
grid('on');
line([0.3 0.3],[-100 0],'color','red','LineStyle','--');
line([0.7 0.7],[-100 0],'color','red','LineStyle','--');
line([0 0.3],[-40 -40],'color','red','LineStyle','--');
line([0.7 1],[-40 -40],'color','red','LineStyle','--');
hold off


% Plot impulse response of BPF
figure(10)
stem((0:length(h)-1),h);
title('BPF Impulse Response')
xlabel('Samples (n)'); ylabel('Amplitude');

