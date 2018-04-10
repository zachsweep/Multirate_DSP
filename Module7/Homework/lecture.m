%% Problem 1
% Design a 2-stage downsampler for M=16 that minimizes the overall
% computation rate. Use the Parks-McClellan design algorithm for each
% filter stage. The requirements for the anti-aliasing filter are as
% follows:

clear;
close all;


% Decimation factor
M = 50;

% Passband/stopband requirements
wp = (7/8)*(pi/M);
ws = pi/M;

% Ripple parameters
d1 = 0.01;
d2 = 0.001;

% Model filter
[n,f0,m0,w] = firpmord([wp/pi ws/pi],[1 0],[d1 d2]);
b = firpm(n,f0,m0,w);

% Plot model filter magnitude response
figure(1)
[H,w] = freqz(b,1,10000);
plot(w/pi,20*log10(abs(H)));
title('Model Filter Magnitude Response');xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');


% COMPUTE MULTISTAGE DECIMATOR FILTER PARAMETERS
% Find stretch factors that minimize total computational cost
M1 = 2:M-1;
M1 = M1(~(rem(M, M1)));
M2 = flip(M1);

D = (-10*log10(d1/2*d2) - 13)/2.324;
cost = zeros(1,length(M1));
for i=1:length(M1)
    Ng = D/(M1(i)^2*M2(i)*(ws-wp));
    Ni = D/((2*pi)-M1(i)*(ws+wp));
    cost(i) = Ng+Ni;
    fprintf('\n%f\n',cost(i));
end
[min_cost,i] = min(cost);
M1 = M1(i);
M2 = M2(i);
fprintf('\nM1 = %d\nM2 = %d\n',M1,M2);


% STAGE 1 I(z)
wp_i = wp;
ws_i = (2*pi/M1)-ws;
[n,f0,m0,w] = firpmord([wp_i/pi ws_i/pi],[1 0],[d1/2 d2]);
b_i = firpm(142,f0,m0,w);

% Compute frequency response I(z)
[H,w] = freqz(b_i,1,10000);
figure(2)
plot(w/pi,20*log10(abs(H)));
title('I(z) Magnitude Response');xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');

% % Downsample by M1
% b_i_ds = downsample(b_i,M1);


% STAGE 2 G(z)
wp_g = M1*wp;
ws_g = M1*ws;
[n,f0,m0,w] = firpmord([wp_g/pi ws_g/pi],[1 0],[d1/2 d2]);
b_g = firpm(90,f0,m0,w);

% Upsample G(z) impulse response
b_g_interp = upsample(b_g,M1);

% Compute frequency response I(z)
[H,w] = freqz(b_g_interp,1,10000);
figure(3)
plot(w/pi,20*log10(abs(H)));
title('G(z) Magnitude Response');xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');

% b_g_interp_ds = downsample(b_g_interp,M2);

h = conv(b_i,b_g_interp);
% h = downsample(h,M2);
% h = h/sum(h);
[H,w] = freqz(h,1,10000);
figure(4)
plot(w/pi,20*log10(abs(H)));
title('Convolved Magnitude Response');xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
line([0 1],[0.0864 0.0864],'color','red','LineStyle','--');
line([0 1],[-0.0873 -0.0873],'color','red','LineStyle','--');