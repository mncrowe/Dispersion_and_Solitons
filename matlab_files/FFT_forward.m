function f_hat = FFT_forward(f,d,Norm)
% Performs a 1 or 2 dimensional Fourier transform using the 'fft' and 'fftshift' functions
% - f: function (vector), may be 1 or 2 dimensional
% - d: dimension (scalar), use to transform only in a single dimension for a multi-dimension array
% - Norm: normalisation type; 1 - divide by N (default), 0 - unscaled MATLAB tranform

if nargin < 2; d = 0; end
if nargin < 3; Norm = 1; end

if sum(size(f)>1)==1
	f_hat = fftshift(fft(f));
    if Norm == 1; C = 1/length(f); else; C = 1; end
end

if sum(size(f)>1)==2 && d==0
    f_hat = fftshift(fft2(f));
    if Norm == 1; C = 1/numel(f); else; C = 1; end
end

if sum(size(f)>1)>1 && d>0
    f_hat = fftshift(fft(f,[],d),d);
    if Norm == 1; N = size(f); C = 1/N(d); else; C = 1; end
end

if sum(size(f)>1)>2 && d==0
    error('Too many or not enough dimensions.')
end

f_hat = C*f_hat;

end

