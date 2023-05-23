function Qwalk( )
%QWALK Simulation of a Bloch-Zener Quantum Walk

%% Input variables

ACDC = 0;                   % Tag that determines whether calculation is AC(0) or DC(1)

n_T = 3;                    % Number of field oscillations
n_dt = 20;                  % Number of time steps within period of field oscillation
dt = 50;                    % Time step

%E0 = 0.0;                  % Initial energy of the wavepacket
Fmax = 25.0;                % Drop of potential over chain
t = 1.0;                    % Hopping amplitude

a = 2.5*t+1.2*Fmax;         % Defines bounds of the spectrum [-a,a]
M = 2^12;                   % Chain length
Vlat = [0.05];               % Lattice dimerization potential

%% Initialisation
sigsq = 80^2;              % Width of Gaussian
x0 = M/2;                   % Initial position of wavepacket
%k = -acos(1-E0/2);          % Wave packet momentum
x = (1:M)';                 % Chain grid
vtot = zeros(M,n_T*n_dt+1); % Keeps track of wavefunction

% Initial wavepacket phi(x)
v = 1/(sqrt(sqrt(2*pi*sigsq)))*exp(-(x-x0*ones(M,1)).^2/(4*sigsq));
vtot(:,1) = v;

%% Main Calculation

% Create Hamiltonian for positive and negative field
Hp = formMatrix(-Fmax,a,M,Vlat,t);
if ACDC==1      % For DC calculation, Hamiltonian is the same
    Hn = Hp;
else
Hn = formMatrix(Fmax,a,M,Vlat,t);
end

% Start updating wavefunction for n_T field oscillations
for i = 1:n_T
    %i
    % Positive field
    for j = 1:n_dt/2
        v = evolve(v,Hp,a,dt);
        vtot(:,(i-1)*n_dt+j+1) = v;
    end
    
    % Negative field
    for j = n_dt/2+1:n_dt
        v = evolve(v,Hn,a,dt);
        vtot(:,(i-1)*n_dt+j+1) = v;
    end
end

%% Plot the resulting wavefunction trajectory

vtot = conj(vtot).*vtot;
vplot = (vtot(1:2:M,:)+vtot(2:2:M,:))/2;    % Taking average over unit cell
wf = real(vtot);

figure;
contourf(wf,'LineColor','none','LevelList',0:0.001:0.04);
title('Wave-packet trajectory')
xlabel('time step')
ylabel('position')

end

%% FORMMATRIX Builds the Hamiltonian

function [H] = formMatrix(Fin,a,M,Vlat,t)
    % In case of a dimer chain
    if length(Vlat) == 1
        D1 = ((1:M)'*Fin/M+(-1).^(1:M)'*Vlat)/a; %sets the dimerization potential
    
    % In case there are more atoms in the chain
    else
        D1 = ((1:M)*Fin/M+Vlat(mod(0:M-1,length(Vlat))+ones(1,M)))'/a; %sets the on-site potential
    end
    
    D2 = (-t/a)*ones(M,1);  %sets the kinetic energy
    
    H = spdiags([D2 D1 D2],[-1 0 1],M,M);
end

%% EVOLVE Calculates evolution of wavefunction

function [v_out] = evolve(v_in,H,a,dt)

N = 2*(a*dt); %automatically evaluates the number of terms in the Chebyshev expansion

v0 = v_in;                          % |Phi_0> = |Phi(0)>
c = besselj(0,a*dt);                % c_0 = J_0(a*dt)
v_out = c*v0;
v1 = H*v0;                          % |Phi_1> = H*|Phi(0)>
c = (-1i)*besselj(1,a*dt);          % c_1 = (-i)*J_1(a*dt)
v_out = v_out+2*c*v1;               % |Phi(t)> = c_0 + 2*c_1*|Phi_1> + ...

k = 2;
while k < N
    v2 = 2*H*v1-v0;                 % |Phi_{k}> = 2*H*|Phi_{k-1}> - |Phi_{k-2}>
    c = (-1i)^k*besselj(k,a*dt);    % c_{k} = (-i)^k*J_k(a*dt)
    v_out = v_out + 2*c*v2;
    v0 = v1;
    v1 = v2;
    k = k+1;
end
end