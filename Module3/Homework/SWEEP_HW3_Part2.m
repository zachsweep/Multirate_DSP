%% Homework #3, Part 2

%........................... REPEAT OF PART 1
wp = 0.45*pi;
ws = 0.55*pi;
Rp = 0.1; %dB
Rs = 30;  %dB

[N,wpnew] = ellipord(wp/pi,ws/pi,Rp,Rs);
[z,p,k]   = ellip(N,Rp,Rs,wpnew);

%............................ START OF PART 2
% 2nd Order All-Pass Filter
AMaz = poly([p(3) p(4)]);
AMbz = flip(AMaz);

% 3nd Order All-Pass Filter
ANaz = poly([p(1) p(2) p(5)]);
ANbz = flip(ANaz);

disp('2nd Order All-Pass Numerator:');
disp(AMbz);
disp('2nd Order All-Pass Denominator:');
disp(AMaz);

disp('3rd Order All-Pass Numerator:');
disp(ANbz);
disp('3rd Order All-Pass Denominator:');
disp(ANaz);

figure(3)
hold on
phasez(ANbz,ANaz);
phasez(AMbz,AMaz);
title('LPF (parallel all-pass)');
xlabel('Normalized Frequency'); ylabel('Phase Response');
legend('AN(z)','AM(z)');
hold off