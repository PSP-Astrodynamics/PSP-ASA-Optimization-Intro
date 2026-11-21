g = 9.81;

% Dynamics
syms r v [2 1]
syms theta w m t
x = [r; v; theta; w; m];

% Rotation and constants
syms I isp L

% Controls
syms T [2 1]
u = T;

r_dot = v;
Rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
v_dot = [0; -g] + Rot*(T / m);
theta_dot = w;
tor = cross([-L; 0; 0], [T; 0]);
w_dot = tor(3) / I;
m_dot = -sqrt(T(1)^2 + T(2)^2) / (isp*g);

x_dot = [r_dot; v_dot; theta_dot; w_dot; m_dot];

matlabFunction(x_dot,"File","dynamics","Vars", [{t}; {x}; {u}; {L; I; isp}]);

j_x = jacobian(x_dot, x);

j_u = jacobian(x_dot, u);

matlabFunction(j_x,"File","Jacobian_x","Vars", [{t}; {x}; {u}; {L; I; isp}]);
matlabFunction(j_u,"File","Jacobian_u","Vars", [{t}; {x}; {u}; {L; I; isp}]);
