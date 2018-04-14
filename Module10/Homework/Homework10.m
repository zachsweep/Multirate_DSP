clear all;

% H0(z) specifications
% Rp = 0.5;   %dB
% Rs = 30;    %dB
Rp = 0.1;
Rs = 50;

Rp_linear = min([(1-10^(-Rp/20)) (10^(Rp/20)-1)]);

% Transition band width
dw = 0.1*pi;

M = 2;
% wp = (pi/6);
% ws = (pi/4);
wp = 0.4385*pi;
ws = 0.6385*pi;

f = pi;
fs = 2*pi;

[n,fo,mo,w] = firpmord( [wp/pi ws/pi], [1 0], [Rp_linear 10^(-Rs/20)]);
b0 = firpm(n,fo,mo,w);
b1 = b0.*exp(1j*pi*(0:length(b0)-1));

% b0 = downsample(b0,M);
% b1 = downsample(b1,M);
H0 = freqz(b0,10000);
H1 = freqz(b1,10000);

He = abs(H0).^2 + abs(H1).^2;

figure(1)
hold on
freqz(b0,10000);
hold on
freqz(b1,10000);

figure(2)
plot(20*log10(abs(He)));
