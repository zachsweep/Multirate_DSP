% clear;

M = 4;
Fs = 100e3;
T = 1/Fs;

f = 10e3;

% w = 2*pi*f/Fs*M;


% wm = ws-wp;

y = (f - ((-(M-1):M-1)/(M*T)));
disp(y);


% w0 = (pi/50); % 0.0628
% w1 = (pi/50)+(2*pi/5); % 21/50 1.3195
% w2 = (pi/50)+(4*pi/5); % 41/50 2.5761

% syms w;

% eqn = ((w) - (2*pi*(-(M-1):M-1))/M)

% w = w1;

% format rat;
% eval(eqn)