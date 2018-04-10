%% Homework #3, Part 4

%........................... REPEAT OF PART 1
wp = 0.45*pi;
ws = 0.55*pi;
Rp = 0.1; %dB
Rs = 30;  %dB

[N,wpnew] = ellipord(wp/pi,ws/pi,Rp,Rs);
[z,p,k]   = ellip(N,Rp,Rs,wpnew);

% Calculate Z domain numerator polynomial
Bz = poly(z);
Bz = k*Bz;

% Calculate Z domain denominator polynomial
Az = poly(p);

% Calculating Magnitude Response
[H_orig,wh]   = freqz(Bz,Az);

%.......................... REPEAT OF PART 2
% 2nd Order All-Pass Filter
AMaz = poly([p(3) p(4)]);
AMbz = flip(AMaz);

% 3nd Order All-Pass Filter
ANaz = poly([p(1) p(2) p(5)]);
ANbz = flip(ANaz);

%.......................... START OF PART 4
K1 = tf2latc(ANaz);
K2 = tf2latc(AMaz);

disp('All-Pass filter 1 multiplier coefficients (3rd Order):')
disp(K1)

disp('All-Pass filter 2 multiplier coefficients (2nd Order):')
disp(K2)


%% Homework #3, Part 5

%...........16 bit Quantization

% Quantize direct form LPF coefficients and Compute
% Frequency response of quantized direct form LPF
[H_dir,wh]   = freqz(quantize(Bz,16),quantize(Az,16));

% Quantize Latice form LPF coefficients and
% Compute Latice realization transfer function for both all-pass sections
[num1,den1]=latc2tf(quantize(K1,16),'allpass');
[num2,den2]=latc2tf(quantize(K2,16),'allpass');

% Compute phase response for each all-pass section
[phAN,wn] = phasez(num1,den1);
[phAM,wm] = phasez(num2,den2);

% Compute frequency response of parallel all-pass sections
H_lat = (exp(1j*phAN)+exp(1j*phAM))/2;

figure(6)
hold on
plot(wh/pi,20*log10(abs(H_orig)));
plot(wh/pi,20*log10(abs(H_dir)));
plot(wn/pi,20*log10(abs(H_lat)));
title('LPF with 16 bit Quantization');
xlabel('Normalized Frequency'); ylabel('Magnitude Response');
legend('Direct form LPF (No Quantization)','Direct form LPF (I16)','Latice form LPF (I16)') 


%...........10 bit Quantization

% Quantize direct form LPF coefficients and Compute
% Frequency response of quantized direct form LPF
[H_dir,wh]   = freqz(quantize(Bz,10),quantize(Az,10));

% Quantize Latice form LPF coefficients and
% Compute Latice realization transfer function for both all-pass sections
[num1,den1]=latc2tf(quantize(K1,10),'allpass');
[num2,den2]=latc2tf(quantize(K2,10),'allpass');

% Compute phase response for each all-pass section
[phAN,wn] = phasez(num1,den1);
[phAM,wm] = phasez(num2,den2);

% Compute frequency response of parallel all-pass sections
H_lat = (exp(1j*phAN)+exp(1j*phAM))/2;

figure(7)
hold on
plot(wh/pi,20*log10(abs(H_orig)));
plot(wh/pi,20*log10(abs(H_dir)));
plot(wn/pi,20*log10(abs(H_lat)));
title('LPF with 16 bit Quantization');
xlabel('Normalized Frequency'); ylabel('Magnitude Response');
legend('Direct form LPF (No Quantization)','Direct form LPF (I16)','Latice form LPF (I16)') 