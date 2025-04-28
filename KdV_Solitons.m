% Solves the KdV equation:

% u_t + 6uu_x + u_xxx = 0

addpath("matlab_files\")
%close all;
clear

% create grids in real space and fourier space
Lx = 10;
Nx = 512;
[x, k] = FFT_grid(Nx, [-Lx Lx]);

% define time domain for plotting animation
T = 5;
Nt = 201;
dt = 0.002;
t = linspace(0, T, Nt);

% define initial condition

f0 = CreateSoliton(x, 4, 0) + CreateSoliton(x, 9, -6);

% define nonlinear term in Fourier space

%dealias = [zeros(ceil(Nx/6), 1); ones(Nx - 2*ceil(Nx/6), 1); zeros(ceil(Nx/6), 1)];
NL = @(f) -1i * 3 * k' .* FFT_forward(FFT_inverse(f).^2); % .* dealias;

% solve using pseudo spectral method

f = real(FFT_inverse(timestep(dt, t, eye(Nx), 1i*diag(k.^3), @(f, t) NL(f), FFT_forward(f0'), 3), 1));

% plot animation
figure;
for i = 1:length(t)
    f_t = f(:, i);
    clf;
    plot(x, f_t, 'k'); xlabel('x'); ylim([0 6]); xlim([-Lx Lx])
    title(['t = ' sig_fig_str(t(i), 3)])
    drawnow;
    pause(0.01)
end

% define soliton

function f = CreateSoliton(x, c, x0)

    f = 1/2 * c * sech(sqrt(c)/2 * (x - x0)).^2;

end