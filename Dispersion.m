% Plots the evolution of a wave-packet for a given dispersion relation
%
% An initial condition of phi0(x) represents our initial wave, this is
% Fourier transformed to get our initial distribution of wavenumbers, k.
% Each wavenumber is then evolved forwards in time using the dispersion
% relation, w = w(k), and the wave at a given time, t, is then calculated
% by an inverse Fourier transform of the evolved wavenumber distribution.

% load FFT_grid, FFT_forward, FFT_inverse, sig_fig_str and sig_fig
addpath("matlab_files\")
close all;
clear

% create grids in real space and fourier space
Lx = 20;
Nx = 512;
[x, k] = FFT_grid(Nx, [-Lx Lx]);

% define time domain for plotting animation
T = 1;
Nt = 101;
t = linspace(0, T, Nt);

% wave parameters, dispersion relation and initial condition
omega = @(k) 20*k;                      % non-dispersive
%omega = @(k) 20*k - k.^3;               % dispersive
phi0 = @(x) exp(-x.^2); %.*exp(1i*x);   % initial oscillations with a gaussian envelope

% calculate phi at a given time using Fourier transform
phi = @(t) real(FFT_inverse(FFT_forward(phi0(x)).*exp(-1i*omega(k)*t)));

% plot animation
figure;
for i = 1:length(t)
    phi_t = phi(t(i));
    clf;
    plot(x, phi_t, 'k'); xlabel('x'); ylim([-1 1])
    title(['t = ' sig_fig_str(t(i), 3)])
    drawnow;
    pause(0.01)
end
    