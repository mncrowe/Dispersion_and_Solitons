function f = FFT_inverse(f_hat,d,Norm)
% Performs a 1 or 2 dimensional inverse Fourier transform using the 'ifft' and 'fftshift' functions
% - f_hat: function (vector), may be 1 or 2 dimensional
% - d: dimension (scalar), use to transform only in a single dimension for a multi-dimension array
% - Norm: normalisation type; 1 - multiply by N (default), 0 - unscaled MATLAB tranform

if nargin < 2; d = 0; end
if nargin < 3; Norm = 1; end

if sum(size(f_hat)>1)==1
	f = ifft(fftshift(f_hat));
    if Norm == 1; C = length(f_hat); else; C = 1; end
end

if sum(size(f_hat)>1)==2 && d==0
    f = ifft2(fftshift(f_hat));
    if Norm == 1; C = numel(f_hat); else; C = 1; end
end

if sum(size(f_hat)>1)>1 && d>0
    f = ifft(fftshift(f_hat,d),[],d);
    if Norm == 1; N = size(f_hat); C = N(d); else; C = 1; end
end

if sum(size(f_hat)>1)>2 && d==0
    error('Too many or not enough dimensions.')
end

f = C*f;

end

