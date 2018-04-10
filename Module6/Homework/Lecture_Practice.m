%% Lecture Practice
clear;

% Decimation Parameter
M = 3;

% Transition band
dw = pi/50;

% Passband/Stopband parameters
wp = pi/M - dw;
ws = pi/M + dw;

% Passband/Stopband ripple parameters
Rp =  1;    %dB
Rs = 40;    %dB

Rp_linear = min([(1-10^(-Rp/20)) (10^(Rp/20)-1)]);

R = [Rp_linear 10^(-Rs/20)];

% Compute estimated filter order and other design params using firpmord()
[N,fo,mo,w] = firpmord([wp/pi ws/pi], [1 0], R);

% N = 23;
% Increase filter order by 2 since original order was under estimated
N = 57;
b = firpm(N,fo,mo,w);


% Compute Frequency Response
[H,w] = freqz(b);


% Plot Impulse Response
figure(1)
stem((0:N),b,'MarkerFaceColor','b');
% axis([0 N -0.2 0.6]);
ylim([-0.2 0.6]);
title('Impulse Response');
xlabel('Samples (n)'); ylabel('Amplitude');


% Plot Magnitude Frequency Response
figure(2)
plot(w/pi,20*log10(abs(H)),'-o','MarkerFaceColor','b');
axis([0 1 -60 10]);
title('Magnitude Response');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude (dB)');
line([0 wp/pi],[-Rp -Rp],'color','red','LineStyle','--');
line([0 wp/pi],[Rp Rp],'color','red','LineStyle','--');
line([ws/pi 1],[-Rs -Rs],'color','red','LineStyle','--');


h = b;

% el = zeros(
clear e0 e1 e2 e3 e4;

% 
% el = zeros(M, (length(h)/M));
% for n = 0:fix(length(h)/M)-1
%     
% %     if (n*M+1 <= length(h))
% %         e0(n+1) = h(n*M+1);
%         el(1,n+1) = h(n*M+1);
% %     end
%     
%     if (n*M+1+1 <= length(h))
% %         e1(n+1) = h((n*M+1+1));
%         el(2,n+1) = h((n*M+1+1));
%     end
%     
%     if (n*M+2 <= length(h))
% %         e2(n+1) = h((n*M+2+1));
%         el(3,n+1) = h((n*M+2+1));
%     end
%     
% end



e0 = h(1:3:end);
e1 = h(2:3:end);
e2 = h(3:3:end);

ee = zeros(1, length(e0)+length(e1)+length(e2));
ee(1:3:end) = e0;
ee(2:3:end) = e1;
ee(3:3:end) = e2;

% ee = resphape(ee,M,length(ee)/3)