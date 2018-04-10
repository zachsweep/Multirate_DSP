function xq = quantize(x,Nbits)
%
% Function to quantize input x to Nbits (including sign bit)
%

% Use MATLAB uencode/udecode functions
% fullscale = max(abs(x));
% xq = udecode(uencode(x,Nbits,fullscale),Nbits,fullscale);

% Determine full-scale value and quantization step size
% Note that xmax > max(abs(x)) due to the addition of eps
xmax = max(abs(x))+eps;
delta = 2*xmax/(2^Nbits);
% Mid-tread (symmetric range, includes zero and 1 extra quantization level)
% xq = delta*sign(x).*floor(abs(x)/delta + 0.5);
% Mid-tread (asymmetric range, includes zero)
xq = delta*floor(x/delta);
% Mid-rise (symmetric range, no zero value)
% xq = delta*(floor(x/delta) + 0.5);