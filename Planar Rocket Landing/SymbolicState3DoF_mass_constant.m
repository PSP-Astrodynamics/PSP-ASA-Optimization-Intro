%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter + Adarsh Rangayyan
% Created On: 28 April, 2025
% Description: 3DoF rocket landing dynamics with changing mass
% Most Recent Change: 28 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%g = 3.7114e-3; % [km / s2]
g = 9.81e-3; % [km / s2]

syms L I alpha m;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
theta = sym("theta", 1);
w = sym("w", 1);
x = [r;v;theta;w];

thrust = sym("thrust", [2, 1]); % Thrust 
u = [thrust];
thrust_mag = sqrt(u(1) .^ 2 + u(2) .^ 2);



rdot = v;
thetadot = w;
M = cross([-L; 0; 0], [thrust; 0]);
wdot = M(3) / I;

rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
a_T_e = rotationMatrix * thrust/m;
vdot = a_T_e - [0; g];

xdot = [rdot; vdot; thetadot; wdot];

%j_a = jacobian(xdot, x);
%j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","SymDynamics3DoF_mass_constant_noumag","Vars", [{t}; {x}; {u}; {L; I; m}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
%matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
%matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);