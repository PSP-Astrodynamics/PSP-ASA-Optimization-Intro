%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA Optimization Intro
% Planar Rocket Landing Symbolic Dynamics
% Author: Travis Hastreiter
% Created On: 7 November, 2026
% Description: 3DoF planar rocket landing dynamics with changing mass
% Most Recent Change: 14 November, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define physical constants
g = 9.81; % [m / s2]

% Create symbolic variables
syms L I Isp;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
theta = sym("theta", 1);
w = sym("w", 1);
m = sym("m", 1);
x = [r; v; theta; w; m];

T = sym("T", [2, 1]);
u = T; % Control is thrust vector

% Define dynamics
rdot = v;
thetadot = w;
M = cross([-L; 0; 0], [T; 0]);
wdot = M(3) / I;

R = [cos(theta), -sin(theta);
     sin(theta), cos(theta)];
vdot = R * T / m - [0; g];

mdot = -sqrt(T(1)^2 + T(2)^2) / (Isp * g);

xdot = [rdot; vdot; thetadot; wdot; mdot];

A = jacobian(xdot, x);
B = jacobian(xdot, u);

% Create equations of motion function for integration
matlabFunction(xdot,"File","Planar Rocket Landing/Planar Rocket Landing Custom/SymDynamicsPlanar3DoF","Vars", [{t}; {x}; {u}; {L; I; Isp}]);

% Create jacobians for linearization and discretization
matlabFunction(A,"File","Planar Rocket Landing/Planar Rocket Landing Custom/SymAJacobianPlanar3DoF","Vars", [{t}; {x}; {u}; {L; I; Isp}]);
matlabFunction(B,"File","Planar Rocket Landing/Planar Rocket Landing Custom/SymBJacobianPlanar3DoF","Vars", [{t}; {x}; {u}; {L; I; Isp}]);