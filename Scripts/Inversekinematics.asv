clc;
clear;
close all;

%% Exponential coordinate

q = sym('q', [6, 1], 'real');

l = sym('l', [6, 1], 'real');

syms d_1 a_3

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

%% invers kinematics;

T_d = sym(eye(4,4));
T_Ei_inv = SerialRobot.inverseTrasformMatrix(T_i);
T_1 = T_d*T_Ei_inv;
% calculating q_1
p_56 = [0, -l(1)+a_3-l(4), -d_1+l(2)+l(3)+l(5), 1]';
p_1 = T_1*p_56;
O_1 = [0, 0, p_1(3), 1]';
p_2 = [sqrt(p_1(1)^2 + p_1(2)^2 - (l(1)+l(4)-a_3)^2);
        -l(1)-l(4)+a_3;
        p_1(3);
        1];

u_1 = p_2-O_1;
v_1 = p_1-O_1;

s_1 = dot(cross(u_1,v_1 ),[0 0 1 0])/(norm(u_1)*norm(v_1));
c_1 = dot(u_1,v_1)/(norm(u_1)*norm(v_1));
q_1 = atan2(s_1,c_1);

% calculating q_5 q_6
e_z1t_inv = SerialRobot.inverseTrasformMatrix(SerialRobot.twist2Exp(zeta(:,1),q_1));
T_56 = e_z1t_inv*T_1;

M_1 = T_56(1,1);
M_2 = T_56(2,1);
M_3 = T_56(3,1);
N_2 = T_56(2,2);
G_2 = T_56(2,3);

q_5 = atan2(sqrt(1-N_2^2), N_2);

q_6 = atan2(-G_2/sin(q_5), M_2/sin(q_5));

% delta = (cos(q_5)*cos(q_6))^2 + sin(q_6);
% s234 = (M_3*cos(q_5)*cos(q_6) - M_1*sin(q_6))/delta; 
% c234 = (M_1*cos(q_5)*cos(q_6) + M_3*sin(q_6))/delta;

% sum_234 = atan2(s234, c234);

e_z5t_inv = SerialRobot.inverseTrasformMatrix(SerialRobot.twist2Exp(zeta(:,5),q_5));
e_z6t_inv = SerialRobot.inverseTrasformMatrix(SerialRobot.twist2Exp(zeta(:,6),q_6));

T_234 = T_56*e_z6t_inv*e_z5t_inv;

m_1 = T_234(1,1);
m_3 = T_234(3,1);
h_1 = T_234(1,4);
h_3 = T_234(3,4);

sum_234 = atan2(m_3, m_1);

f_1 = (d_1+l(2)+l(3))*sin(sum_234) - h_1;
f_2 = h_3 - d_1 + (d_1+l(2)+l(3))*cos(sum_234);

c3 = (f_1^2 + f_2^2 -l(2)^2-l(3)^2)/(2*l(2)*l(3));

q_3 = atan2(sqrt(1-c3^2),c3);

Delta = (l(2)+l(3)*cos(q_3))^2 + (l(3)*sin(q_3))^2;

s2 = ((l(2)+l(3)*cos(q_3))*f_1 - (l(3)*sin(q_3))*f_2)/Delta;
c2 = ((l(2)+l(3)*cos(q_3))*f_2 + (l(3)*sin(q_3))*f_1)/Delta;

q_2 = atan2(s2,c2);

q_4 = sum_234 - q_2 - q_3;
%%

z = sym(eye(4,4));
for i = 2:4
z = simplify(z*SerialRobot.twist2Exp(zeta(:,i),q(i)),'Steps',70);
end

%%
z

