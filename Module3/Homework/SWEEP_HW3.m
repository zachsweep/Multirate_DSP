%% Homework #3, Part 1

wp = 0.45*pi;
ws = 0.55*pi;
Rp = 0.1; %dB
Rs = 30;  %dB

[N,wpnew] = ellipord(wp/pi,ws/pi,Rp,Rs);
[z,p,k]   = ellip(N,Rp,Rs,wpnew);

% Calculate Z domain numerator polynomial
Bz = poly(z);
Bz = k*Bz;
disp('Numerator Coefficients:\n');
disp(Bz);

% Calculate Z domain denominator polynomial
Az = poly(p);
disp('Denominator Coefficients:\n');
disp(Az);

% Calculating Magnitude Response
[H_orig,wh]   = freqz(Bz,Az);

% Calculating Group Delay
[Gpd,wg] = grpdelay(Bz,Az);

% Plotting Magnitude Response and Group Delay
figure(1) 
plot(wh/pi,20*log10(abs(H_orig)),'-');
title('LPF (direct form)');
xlabel('Normalized Frequency'); ylabel('Magnitude Response');

figure(2)
plot(wg/pi,Gpd,'-');
title('LPF (direct form)');
xlabel('Normalized Frequency'); ylabel('Group Delay (s)');



