clc;
clear;
close all;

%%

q = sym('q', [6, 1], 'real');

l = sym('l', [6, 1], 'real');

syms d_1 a_3


DH_parameter = [   0,    0,   d_1,          q(1)  
                pi/2,    0,  l(1),   pi/2 + q(2) 
                   0, l(2),  -a_3,          q(3)
                   0, l(3),  l(4),  -pi/2 + q(4)
               -pi/2,    0,  l(5),          q(5)
                pi/2,    0,  l(6),          q(6)];

w = [0     0     0     0     0     0
     0    -1    -1    -1     0    -1
     1     0     0     0     1     0];

p = [0,   0,        0,             0,             0,                  0
     0,   0,        0,             0, a_3-l(1)-l(4),                  0
     0, d_1, d_1+l(2), d_1+l(2)+l(3),             0, d_1+l(2)+l(3)+l(5)];
v = sym(zeros(3,6));

for  i = 1:6
    v(:,i) = -cross(w(:,i),p(:,i));
end
zeta = [v;w];
T_i =  [1, 0,  0,                  0
        0, 0, -1,  a_3-l(1)-l(4)-l(6)
        0, 1,  0,  d_1+l(2)+l(3)+l(5)
        0, 0,  0,                  1];
%%

robot = SerialRobot(q,DH_parameter);
%%

T_E = T_i;
z = sym(eye(4,4));
for i = 1:6
    z = simplify(z*SerialRobot.twist2Exp(zeta(:,i),q(i)),"Steps",20);
end
T_E = z*T_E;
%%
clc
i = 1;
j= 1;



