clear; close('all');

%% Digital filter specifications (lowpass)
ws1 = 0.2*pi;
wp1 = 0.3*pi;
wp2 = 0.7*pi;
ws2 = 0.8*pi;
Gp  = 0.99;
Gs  = 0.01;

% Pre-warp passband/stopband frequencies using the bilinear transformtion
Wp1 = tan(wp1/2);
Ws1 = tan(ws1/2);
Wp2 = tan(wp2/2);
Ws2 = tan(ws2/2);

% Compute passband/stopband ripple parameters
Rp = -20*log10(Gp);
Rs = -20*log10(Gs);
ep = sqrt(1/(Gp*Gp) - 1); es = sqrt(1/(Gs*Gs) - 1);

%% Elliptic response type
% Compute analog filter order and revised passband frequency
[N,Wpnew] = ellipord(Wp1,Ws1,Rp,Rs,'s');
[Bs,As] = ellip(N,Rp,Rs,Wpnew,'s');
% Find s-plane poles/zeros
s_zeros = roots(Bs);
s_poles = roots(As);
Nzeros = length(s_zeros);
Npoles = length(s_poles);
% Plot poles/zeros in s-plane
figure(1)
plot(real(s_zeros),imag(s_zeros),'o',real(s_poles),imag(s_poles),'x')
hold on;
z = axis;
% line([0,0],[z(3),z(4)],'Color',[0 0 0]);
% line([z(1),z(2)],[0,0],'Color',[0 0 0]);
xlabel('Real'); ylabel('Imag'); title('s-plane')
axis('equal')
% Map poles/zeros to z-plane using bilinear transformation
z_zeros = (1 + s_zeros)./(1 - s_zeros);
z_poles = (1 + s_poles)./(1 - s_poles);
% Note that since the numerator and denominator order for H(s) are equal there are
% no zeros at s=inf to map to z=-1
% The following code attempts to handle the case for which the order of the
% numerator and denominator are not equal
if(Npoles ~= Nzeros)
    % For this case we need to add poles or zeros at z=-1
    if(Npoles == Nzeros + 1)
        % Add zero at z=-1
        z_zeros = [z_zeros;-1];
        fprintf('Zero added to H(z) at z=-1 \n')
    end
    % Other cases ???
end
% Plot poles/zeros in z-plane
figure(2)
plot(real(z_zeros),imag(z_zeros),'o',real(z_poles),imag(z_poles),'x')
hold on;
% Plot unit circle
phi = linspace(0,2*pi,1000);
plot(cos(phi),sin(phi),'k-');
hold off
xlabel('Real'); ylabel('Imag'); title('z-plane')
axis('equal')
% Find transfer function from poles/zeros
Bz = poly(z_zeros);
Az = poly(z_poles);
% Adjust DC gain
gain = Gp*sum(Az)/sum(Bz);
Bz = gain*Bz;
% Plot frequency response
w = linspace(0,pi,1000);
H = freqz(Bz,Az,w);
% Calculate group delay
delay = grpdelay(Bz,Az,w);
figure(3) 
subplot(211)
plot(w/pi,abs(H),'-')
xlabel('Normalized Frequency'); ylabel('Magnitude Response');
subplot(212)
plot(w/pi,delay,'-')
xlabel('Normalized Frequency'); ylabel('Group Delay (s)');

%% Alternative approach
% Note that both ellipord and ellip require the passband/stopband
% frequencies normalized to pi
[N,wpnew] = ellipord(wp1/pi,ws1/pi,Rp,Rs);
[Bz,Az] = ellip(N,Rp,Rs,wpnew);
figure(8)
% Plot frequency response
w = linspace(0,pi,1000);
H = freqz(Bz,Az,w);
plot(w/pi,abs(H),'-');