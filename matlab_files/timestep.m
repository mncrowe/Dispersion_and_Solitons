function [phi,t] = timestep(dt,t,D,L,N,phi0,method,C)
% Solve the system: D*d(phi)/dt = L*phi + N(phi,t) using implicit-explicit (IMEX) Runge-Kutta (RK) schemes
% - dt: timestep, assumed to be constant, (set to min(dt,diff(t)))
% - t: vector of time saves, first and last points define the time interval, t in [t_0 t_end]
% - D: LHS matrix, n x n
% - L: linear RHS matrix, n x n
% - N: nonlinear RHS vector;
%       vector - n x 1 vector
%       function - anonymous function, N = @(phi,t) ...
% - phi0: initial condition for phi
% - method: timestepping method used;
%       1 - 1 stage, IMEX trapezoidal method, 1st order DIRK+ERK scheme (default)
%       2 - 1 stage, RK111, 1st order DIRK+ERK scheme (implicit backward/explicit forward Euler)
%       3 - 2 stage, RK222, 2nd order DIRK+ERK scheme (Ascher et. al. 97)
%       4 - 3 stage, RK332, 2nd order DIRK+ERK scheme (Koto 06)
%       5 - 3 stage, RKSMR, 3-eps order DIRK+ERK scheme (Spalart et. al. 91)
%       6 - 4 stage, RK443, 3nd order DIRK+ERK scheme (Ascher et. al. 97)
%       7 - 1 stage, RK1, 1st order ERK (forward Euler) (not recommended for singular D)
%       8 - 2 stage, RK2, 2nd order ERK scheme (not reccommended for singular D), method depends on C
%       9 - 4 stage, RK4, 4th order ERK scheme (not reccommended for singular D)
%       0 - custom, enter cell array, C, with C{1} = s, C{2} = H, C{3} = A, C{4} = c
% - C: timestepping method parameter, depends on method as;
%       8 - enter alpha for RK2, alpha = 1 (Heun's method), 1/2 (midpoint method, default), 2/3 (Ralston Method)
%       0 - enter cell array of C{1:4} = {s, H, A, c} for custom RK scheme, required for method = 0
%
% ----------------------------------------------------------------------------
% Note: Terminology; DIRK (diagonally implicit Runge-Kutta), ERK (explicit
%       Runge-Kutta), IMEX (implicit-explicit), SMR (Spalart, Moser & Rogers)
%
% Note: H, A, c do not use the full Butcher tableau as sometimes presented.
%       A is the s x s bottom left square block, H is the 2nd to final rows and
%       c is the first s entries only (0 to s-1), so c_0 = 0. This is as any
%       remaining values are zeros or ones and are known. See examples below.
%       Therefore: sum_j A_ij = sum_j H_ij = c_i, for i=1..s.
%       Here c_s = 1 is not included in c. A and H must be lower triangular.
% ----------------------------------------------------------------------------


if nargin < 7; method = 1; end
if method == 8 && nargin < 8; C = 1/2; end

dt = min(min(dt,diff(t)));  % set dt, minimum of input dt and t spacing
NL = length(L);             % size of system of equations
t_grid = t(1):dt:t(end);    % time grid, constant dt
Nt = length(t_grid);        % number of timesteps
save_tol = 1e-2;            % tolerance for saving onto t grid when dt does not divide diff(t)
output = 1;                 % output progress to screen

warning('off','MATLAB:nearlySingularMatrix')    % large matrices usually are ill-conditioned by MATLAB's definition, instead we check for NaNs

% define IMEX Runge-Kutta tableaus, H is s x (s+1), A is s x s, c is 1 x s
if method == 1
    s = 1;
    H = [0.5 0.5];
    A = 1;
    c = 0;
end
if method == 2
    s = 1;
    H = [0 1];
    A = 1;
    c = 0;
end
if method == 3
    s = 2;
    H = [0 1-1/sqrt(2) 0; 0 1/sqrt(2) 1-1/sqrt(2)];
    A = [1-1/sqrt(2) 0; -1/sqrt(2) 1+1/sqrt(2)];
    c = [0 1-1/sqrt(2)];
end
if method == 4
    s = 3;
    H = [0 1 0 0; 0 -1/2 1 0; 0 -1 1 1];
    A = [1 0 0; 1/2 0 0; 0 0 1];
    c = [0 1 1/2];
end
if method == 5
    s = 3;
    H = [29/96 37/160 0 0; 29/96 5/32 5/24 0; 29/96 5/32 3/8 1/6];
    A = [8/15 0 0; 1/4 5/12 0; 1/4 0 3/4];
    c = [0 8/15 2/3];
end
if method == 6
    s = 4;
    H = [0 1/2 0 0 0; 0 1/6 1/2 0 0; 0 -1/2 1/2 1/2 0; 0 3/2 -3/2 1/2 1/2];
    A = [1/2 0 0 0; 11/18 1/18 0 0; 5/6 -5/6 1/2 0; 1/4 7/4 3/4 -7/4];
    c = [0 1/2 2/3 1/2];
end
if method == 7
    s = 1;
    H = [1 0];
    A = 1;
    c = 0;
end
if method == 8
    s = 2;
    H = [C 0 0; 1-1/(2*C) 1/(2*C) 0];
    A = [C 0; 1-1/(2*C) 1/(2*C)];
    c = [0 C];
end
if method == 9
    s = 4;
    H = [1/2 0 0 0 0; 0 1/2 0 0 0; 0 0 1 0 0; 1/6 1/3 1/3 1/6 0];
    A = [1/2 0 0 0; 0 1/2 0 0; 0 0 1 0; 1/6 1/3 1/3 1/6];
    c = [0 1/2 1/2 1];
end
if method == 0
    s = C{1};
    H = C{2};
    A = C{3};
    c = C{4};
end

% decompose LHS implicit terms for fast inversion, convert D and L to sparse
disp('Decomposing matrices...')
for i = 1:s
    disp([num2str(i) ' of ' num2str(s)])
    M{i} = decomposition(D - dt*H(i,i+1)*L);
end
D = sparse(D);
L = sparse(L);

% convert N to a function if constant
if ~isa(N,'function_handle')
    N = @(phi,t) N;
end

% define IC from phi0, set phi = phi0, define K values for RK method
phi_save = phi0.*ones(NL,length(t));
phi = phi0; it_save = 2;
K = zeros(NL,s+1);

% timestep system using RK method
disp('Starting timestepping...')
t1 = tic;
for it = 2:Nt
    t0=t_grid(it-1);
    K(:,1) = phi;
    for i = 1:s                                         % loop through RK stages
        RHS = zeros(NL,1);                              % (re)set RHS to zero
        for j = 1:i                                     % build RHS terms
            RHS = RHS+H(i,j)*L*K(:,j)+A(i,j)*N(K(:,j),t0+c(j)*dt);
        end
        K(:,i+1) = M{i}\(D*K(:,1)+dt*RHS);              % solve implicit system for next K 
    end
    phi = K(:,s+1);                                     % update phi_n -> phi_{n+1}
    if (t_grid(it)+dt*save_tol)>=t(it_save)             % save phi if current time is in t
        phi_save(:,it_save) = phi;
        it_save = it_save+1;
        if output == 1
            disp(['t = ' num2str(t0)])
            counter(toc(t1),(it-1)/(Nt-1))              % print progress to screen
        end
            
    end
    if max(isnan(phi),isinf(phi))==1; error(['NaN error, t = ' num2str(t0)]); end
end
disp(['Timestepping complete, elapsed time = ' num2str(toc(t1))])

phi = phi_save;                                 % define output values

end