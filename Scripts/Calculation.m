clc;
clear;
close all;


%% 1) Project Details

% 6R Robot manipulator without any 3 censecutive axis

%%% numeric value of dynamic and kinematics parameter will be as follow
%%% they'd been added here for later usage

% kinematics parameter
% Note: all value extract from model in solidworks

l_num = [130, 600, 549, 105, 105, 66]'*1e-3; %all value in meter

a_3_num = 126e-3; % value in meter

d_1_num = 76e-3; % value in meter

% Dynamic parameter

I_tensor = zeros(3, 3, 6);
mass = zeros(6,1);
COMrel2DH = zeros(3,6); % center of mass location relative to DH coordinate

%%% Note: all value for moment of inertia are in g.mm^2. they'll be converted to 
%%% kg.m^2 at the end of this section and all masses value in kg also

% Link 1

I_tensor(:,:,1) = [2.71e+07,	   -1.82,	    17.8
	                  -1.82,	2.59e+07,	2.45e+05
	                   17.8,	2.45e+05,	1.56e+07];
mass(1) = 10645.74;

COMrel2DH(:,1) = [0, 0.71, -0.31]';

%Link 2

I_tensor(:,:,2) = [ 7.18e+07,     -96.9,         0
	                   -96.9,  2.42e+09,         0
	                       0,         0,  2.39e+09];
mass(2) = 43347.13;

COMrel2DH(:,2) = [300, 0, -0.15]';

% Link 3

I_tensor(:,:,3) = [9.18e+08,	 1.55e+07,	     12.2
	               1.55e+07,	 2.15e+07,	    -55.1
	                   12.2,	    -55.1,   9.21e+08];
mass(3) = 23361.11; 

COMrel2DH(:,3) = [313.42, 0, 2.24]';


% Link 4

I_tensor(:,:,4) = [1.03e+07,	    -8.19,	     0.98
	                  -8.19,	 7.82e+06,	 -1.32e+05
	                   0.98,	-1.32e+05,	  9.63e+06];
mass(4) = 6442.29;

COMrel2DH(:,4) = [0, 0.3, -1.36]';

% Link 5

I_tensor(:,:,5) = [1.03e+07,	   -8.19,	    -0.98
	                  -8.19,    7.82e+06,	 1.32e+05
                      -0.98,    1.32e+05,	 9.63e+06];
mass(5) = 6442.29;

COMrel2DH(:,5) = [0.3000,  -1.3600, 0]';


% Link 6

I_tensor(:,:,6) = [1.62e+05,           0,           0
	                      0,    1.51e+05,           0
	                      0,           0,    3.06e+05];
mass(6) = 384.69;

COMrel2DH(:,6) = [0,  0, -5]';

I_tensor = I_tensor*1e-9;
COMrel2DH = COMrel2DH/100;
mass = mass/1000;

%% 1.Forward Kinematics

%% 1.1 FK using DH
q = sym('q', [6, 1], 'real');

d_1 = d_1_num;
l = l_num;
a_3 = a_3_num;

% Denavit-Hartenberg parameter
DH_parameter = [   0,    0,   d_1,          q(1)  
                pi/2,    0,  l(1),   pi/2 + q(2) 
                   0, l(2),  -a_3,          q(3)
                   0, l(3),  l(4),  -pi/2 + q(4)
               -pi/2,    0,  l(5),          q(5)
                pi/2,    0,  l(6),          q(6)];


% initializing an object from SerialRobot Class for automating calculation
robot = SerialRobot(q, DH_parameter);

% transformation of End Efecctor frame to the inertial frame
T_E = robot.T_E;

%% 1.2 FK Verification
%%% verification by using a more geometrical approach(exonential
%%% coordinate)

% axis of rotation of each joint in the initial position

w = [0     0     0     0     0     0
     0    -1    -1    -1     0    -1
     1     0     0     0     1     0];

% an arbitary point on the axis of rotation

p = [0,   0,        0,             0,             0,                  0
     0,   0,        0,             0, a_3-l(1)-l(4),                  0
     0, d_1, d_1+l(2), d_1+l(2)+l(3),             0, d_1+l(2)+l(3)+l(5)];

% velocity of points

v = sym(zeros(3,6));

for  i = 1:6
    v(:,i) = -cross(w(:,i),p(:,i));
end

% twists

zeta = [v;w];

% initial transformation of End Efecctor frame to the inertial frame

T_i =  [1, 0,  0,                  0
        0, 0, -1,  a_3-l(1)-l(4)-l(6)
        0, 1,  0,  d_1+l(2)+l(3)+l(5)
        0, 0,  0,                  1];

% calculation of geometric approach
T_E_geom = sym(eye(4,4));
for i = 1:6
    T_E_geom = simplify( ...
        T_E_geom*SerialRobot.twist2Exp(zeta(:,i),q(i))...
        ,'Steps',40);
end
T_E_geom = T_E_geom*T_i;

FKverification = isequal(expand(T_E_geom),expand(T_E));
%% 1.3 ploting work space


X_EE = matlabFunction(T_E_geom(1,end));
Y_EE = matlabFunction(T_E_geom(2,end));
Z_EE = matlabFunction(T_E_geom(3,end));

q_1 = linspace(0, 2*pi, 15);
q_2 = linspace(-pi/2,pi/2,15);
q_3 = linspace(-2.53,2.53,15); % 145 degree
q_4 = q_1;
q_5 = q_1;
q_6 = q_1;

[q_1,q_2,q_3,q_4,q_5,q_6] = ndgrid(q_1,q_2,q_3,q_4,q_5,q_6);

X = X_EE(q_1,q_2,q_3,q_4,q_5,q_6);
Y = Y_EE(q_1,q_2,q_3,q_4,q_5,q_6);
Z = Z_EE(q_2,q_3,q_4,q_5);

scatter3(X(:),Y(:),Z(:),'.')
grid on
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
axis equal

%% Dynamic

% Dynamic properties calculate using Solidworks

I_tensor = zeros(3, 3, 6);

% all value for moment of inertia are in g.mm^2. they'll be converted to 
% kg.m^2 at the end of script
% Link 1

I_tensor(:,:,1) = [2.71e+07,	   -1.82,	    17.8
	                  -1.82,	2.59e+07,	2.45e+05
	                   17.8,	2.45e+05,	1.56e+07];

%Link 2

I_tensor(:,:,2) = [ 7.18e+07,     -96.9,         0
	                   -96.9,  2.42e+09,         0
	                       0,         0,  2.39e+09];

% Link 3

I_tensor(:,:,3) = [9.18e+08,	 1.55e+07,	     12.2
	               1.55e+07,	 2.15e+07,	    -55.1
	                   12.2,	    -55.1,   9.21e+08];

% Link 4

I_tensor(:,:,4) = [1.03e+07,	    -8.19,	     0.98
	                  -8.19,	 7.82e+06,	 -1.32e+05
	                   0.98,	-1.32e+05,	  9.63e+06];

% Link 5

I_tensor(:,:,5) = [1.03e+07,	   -8.19,	    -0.98
	                  -8.19,    7.82e+06,	 1.32e+05
                      -0.98,    1.32e+05,	 9.63e+06];

% Link 6

I_tensor(:,:,6) = [1.62e+05,           0,           0
	                      0,    1.51e+05,           0
	                      0,           0,    3.06e+05];

I_tensor = I_tensor*1e-9;
