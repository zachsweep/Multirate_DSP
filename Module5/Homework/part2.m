%% Problem 2

wp = 0.3*pi;
ws = 0.5*pi;

N = 30;
M = N/2;
alpha = [0.2 0.5];

bb = zeros(length(alpha),M+1);
h = zeros(length(alpha),N+1);
for a=1:length(alpha)
    % Compute Ps Matrix
    syms w;
    cw = cos((0:M)*w).';
    cw = cw*cw.';

    Ps = zeros(M+1,M+1);
    w = linspace(ws,pi,5000);

    for m=0:M
        for n=0:M
            % Evaluate cw(m,n) at each value of w
            c = eval(cw(m+1,n+1));

            % Handling cases where c = 1 or 0
            if c==1
                c = ones(length(w),1);
            elseif c==0
                c = zeros(length(w),1);
            end

            % Compute integral
            Ps(m+1,n+1) = (1/pi)*trapz(w,c);
        end
    end

    % Compute Pp Matrix
    syms w;
    cw = cos((0:M)*w).';
    cw = (1 - cw)*(1 - cw).';

    Pp = zeros(M+1,M+1);
    w = linspace(0,wp,5000);

    for m=0:M
        for n=0:M
            % Evaluate cw(m,n) at each value of w
            c = eval(cw(m+1,n+1));

            % Handling cases where c = 1 or 0
            if c==1
                c = ones(length(w),1);    
            elseif c==0
                c = zeros(length(w),1);
            end

            % Compute Integral
            Pp(m+1,n+1) = (1/pi)*trapz(w,c);
        end
    end

    % Compute P
    P = alpha(a)*Ps + (1-alpha(a))*Pp;

    % Compute Eigen Vectors/Values of P
    [V,D] = eig(P,'vector');

    % Find index of smallest Eigen value in the Eigen value column vector
    ind = find(D==min(D));

    % Find Eigen Vector containing smallest Eigen value using the index
    b = V(:,ind);

    % Re-organize bn to get h(n)
    % h(M) = b(0)
    h(a,M+1) = b(1);
    
    % h(n) = b(n)/2 for n = 1 to M
    h(a,M+2:end) = b(2:M+1).'/2;
    h(a,1:M) = flip(b(2:M+1).'/2);
    
    h(a,:) = h(a,:)/sum(h(a,:));
    
end

% Least Squares Impulse Reponse and Frequency Response
hls = firls(30,[0 .3 0.5 1],[1 1 0 0]);
[Hls,wls] = freqz(hls);

% Magnitude Response for alpha = 0.2
h_02 = h(1,:);
[H_02,w02] = freqz(h_02);

% Magnitude Response for alpha = 0.5
h_05 = h(2,:);
[H_05,w05] = freqz(h_05);


% Plot Impulse Response for alpha=0.2, alpha=0.5, and
% least squares
figure(1)
hold on
stem(h_02);
stem(h_05);
stem(hls);
title('Impulse Response');
xlabel('Samples (n)'); ylabel('Amplitude');
legend('Eigen,alpha=0.2','Eigen,alpha=0.5','Least Squares');

% Zoomed view of some of the impulse response values
figure(2)
hold on
stem(h_02);
stem(h_05);
stem(hls);
axis([15.99 16.01 0.396 0.408]);
title('Zoomed View of h(16)');
xlabel('Samples (n)'); ylabel('Amplitude');
legend('Eigen,alpha=0.2','Eigen,alpha=0.5','Least Squares');


% Plot Magnitude Response for alpha=0.2, alpha=0.5, and
% least squares
figure(3);
hold on
plot(w02/pi,20*log10(abs(H_02)));
plot(w05/pi,20*log10(abs(H_05)));
plot(wls/pi,20*log10(abs(Hls)));
title('Magnitude Response (Eigen Filter VS Least Squares)');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude Response (dB)');
legend('Eigen,alpha=0.2','Eigen,alpha=0.5','Least Squares');
line([0.5 0.5],[-100 10],'color','red','LineStyle','--');
line([0.3 0.3],[-100 10],'color','red','LineStyle','--');

% Zoomed view of Magnitude Response (Passband)
figure(4);
hold on
plot(w02/pi,20*log10(abs(H_02)));
plot(w05/pi,20*log10(abs(H_05)));
plot(wls/pi,20*log10(abs(Hls)));
axis([0 0. -0.1 0.1]);
title('Magnitude Response (Eigen Filter VS Least Squares)');
xlabel('Normalized Frequency (x pi rad/sample)'); ylabel('Magnitude Response (dB)');
legend('Eigen,alpha=0.2','Eigen,alpha=0.5','Least Squares');
line([0.5 0.5],[-100 10],'color','red','LineStyle','--');
line([0.3 0.3],[-100 10],'color','red','LineStyle','--');
