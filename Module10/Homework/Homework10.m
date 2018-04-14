%% Problem 1

clear all;

% H0(z) specifications
Rp = 0.5;   %dB
Rs = 30;    %dB
d1 = 1 - 10^(-Rp/20);
d2 = 10^(-Rs/20);
k = d1/d2;

% wc = 0.5241*pi;
wc = 0.5219*pi;
wp = wc - 0.05*pi;
ws = wc + 0.05*pi;
dw = 0.1*pi;

% Estimate the filter order required to meet specs
Norder = (-10*log10(d1*d2) - 13)/(2.324*dw);
Norder = ceil(Norder);

% Make sure that filter order is odd
if ~mod(Norder,2)
    Norder = Norder+1;
end

% Increase filter order since the estimate was underestimated
Norder = Norder+4;

% Compute H0(z) impulse response
b0 = firpm(Norder,[0 wp ws pi]/pi,[1 1 0 0],[1/k 1]);

% Compute H1(z) impulse response by modulating H0(z) impulse response
b1 = b0.*exp(-1j*pi*(0:length(b0)-1));

% Compute frequency responses H0(z) and H1(z)
N_FFT = 2048;
H0 = fft(b0,N_FFT);
H1 = fft(b1,N_FFT);
w = (0:1023)/1024;

% Compute phase response of H0(z) and H1(z)
[P0,w0] = phasez(b0,2048);
[P1,w1] = phasez(b1,2048);

% Compute 2X the distortion transfer function
T0 = (1/2)*(H0.^2 - H1.^2);
T0 = T0*2;

% Compute impulse response of distortion transfer function
t0 = ifft(T0);

% Plot magnitude response of H0(z) and H1(z)
figure(1)
hold on
plot(w, 20*log10(abs(H0(1:1024))));
plot(w, 20*log10(abs(H1(1:1024))));
title('QMF Filter Bank');xlabel('Normalized Frequency (x pi (rad/sample))');
ylabel('Magnitude (dB)');
legend('H_0(z)','H_1(z)');
line([0 1],[-30 -30],'color','red','LineStyle','--');
line([0 1],[-0.5 -0.5],'color','red','LineStyle','--');
line([0 1],[0.5 0.5],'color','red','LineStyle','--');

% Plot phase response of H0(z) and H1(z)
figure(2)
hold on
plot(w0/pi,P0);
plot(w1/pi,P1);
title('QMF Filter Bank');xlabel('Normalized Frequency (x pi (rad/sample))');
ylabel('Phase (radians)');
legend('H_0(theta)','H_1(theta)');

% Plot magnitude/phase response of distortion transfer function
figure(3)
freqz(t0);
title('Distortion Transfer Function');
legend('T(z)');

%% Problem 2
M = 2;

h0 = b0;
h1 = b1;

% Create each polyphase impulse response by resampling h(n) from problem 1.
e0_0 = M*h0(1:M:end);
e0_1 = M*h0(2:M:end);

e1_0 = M*h1(1:M:end);
e1_1 = M*h1(2:M:end);


% Interleaving polyphase components so that I could plot the impulse
% response and make sure it matched the direct form impulse response.
e0_L = zeros(1, length(e0_0)+length(e0_1));
e0_L(1:M:end) = e0_0;
e0_L(2:M:end) = e0_1;

e1_L = zeros(1, length(e1_0)+length(e1_1));
e1_L(1:M:end) = e1_0;
e1_L(2:M:end) = e1_1;


% Compute magnitude response for each polyphase component
[E0_0,w0_0] = freqz(e0_0);
[E0_1,w0_1] = freqz(e0_1);

[E1_0,w1_0] = freqz(e1_0);
[E1_1,w1_1] = freqz(e1_1);


% Compute phase response for each polyphase component
[P0_0,wp0_0] = phasez(e0_0);
[P0_1,wp0_1] = phasez(e0_1);

[P1_0,wp1_0] = phasez(e1_0);
[P1_1,wp1_1] = phasez(e1_1);

% Plot phase response for e0 and e1
figure(7);
subplot(211)
hold on
plot(wp0_0/pi, P0_0);
plot(wp0_1/pi, P0_1);
title('Polyphase Components');
ylabel('Phase (radians)');
legend('E_0(z)','E_1(z)');
hold off


% Plot magnitude response for e0 and e1
subplot(212)
hold on
plot(w0_0/pi, 20*log10(abs(E0_0)));
plot(w0_1/pi, 20*log10(abs(E0_1)));
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
legend('E_0(z)','E_1(z)');
hold off

% Plot phase response for e0 and e1
figure(8);
subplot(211)
hold on
plot(wp1_0/pi, P1_0);
plot(wp1_1/pi, P1_1);
title('Polyphase Components');
ylabel('Phase (radians)');
legend('E_0(z)','E_1(z)');
hold off


% Plot magnitude response for e0 and e1
subplot(212)
hold on
plot(w1_0/pi, 20*log10(abs(E1_0)));
plot(w1_1/pi, 20*log10(abs(E1_1)));
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
legend('E_0(z)','E_1(z)');
hold off
