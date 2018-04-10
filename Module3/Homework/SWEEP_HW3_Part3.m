%% Homework #3, Part 3
% I got the stopband to be within +-1dB of the LPF stopband by keeping the
% filter order at 5 and changing Rp to 0.00435;

wp = 0.45*pi;
ws = 0.55*pi;

% new Rp
Rp = 0.00435; %dB
Rs = 30;      %dB

% Recompute order and w
[N,wpnew] = ellipord(wp/pi,ws/pi,Rp,Rs);

% fix filter order to 5th order
[z,p,k]   = ellip(5,Rp,Rs,wpnew);


% Recompute allpass filters
AMaz = poly([p(3) p(4)]);
AMbz = flip(AMaz);
ANaz = poly([p(1) p(2) p(5)]);
ANbz = flip(ANaz);
[pAM,wm]  = phasez(AMbz,AMaz);
[pAN,wn]  = phasez(ANbz,ANaz);

% High pass filter response
Hhp = (exp(1j*pAM)-exp(1j*pAN))/2;

% Low pass filter response
Hlp = (exp(1j*pAM)+exp(1j*pAN))/2;

figure(5)
hold on
plot(wn/pi,20*log10(abs(Hhp)));
plot(wn/pi,20*log10(abs(Hlp)));
line([0 wn(end)/pi],[-30 -30],'color','black');
title('HPF (power complement to lowpass)');
xlabel('Normalized Frequency'); ylabel('Magnitude Response(dB)');
legend('HPF','LPF');