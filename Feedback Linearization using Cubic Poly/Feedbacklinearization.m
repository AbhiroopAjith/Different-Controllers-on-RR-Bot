clear
clc
clear all


syms a0 a1 a2 a3
syms to tf t
syms qo qodot qf qfdot

a=[a0;a1;a2;a3]
A = [1 to to^2 to^3;0 1 2*to 3*(to^2);1 tf (tf)^2 tf^3;0 1 2*tf 3*(tf^2)]
B= [qo;qodot;qf;qfdot]
a= inv(A)*B

Cub1 = subs(a,[to,tf,qo,qodot,qf,qfdot],[0,10,pi,0,0,0]);
Cub2 = subs(a,[to,tf,qo,qodot,qf,qfdot],[0,10,pi/2,0,0,0]);
Cub1=double(Cub1)
Cub2=double(Cub2)
%Equation is of the form: a0 + a1*t^+a2*t^2+a3*t^3
%The 2 cubic equations are given by 
Cubicpoly1= Cub1(1)+Cub1(2)*t+Cub1(3)*(t^2)+Cub1(4)*(t^3)
Cubicpoly2= Cub2(1)+Cub2(2)*t+Cub2(3)*(t^2)+Cub2(4)*(t^3)

%Manipulator Equation given by
% % M(q) ̈q + C(q,  ̇q)  ̇q + g(q) = τ
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 
syms theta1_dot theta2_dot theta1_ddot theta2_ddot u1 u2 g

eq1= theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
eq2= theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);
eq1=simplify(eq1);
eq2=simplify(eq2);


%Ode Function
T=10;
[t,y] = ode45(@ode_link, [0,T], [(10*pi)/9,(25*pi)/36,0,0]);
figure(1)
subplot(2,2,1)
plot(t,rad2deg(y(:,1)),'linewidth',2);
xlabel('t(sec)');
ylabel('theta1(deg)');
subplot(2,2,2)
plot(t,rad2deg(y(:,2)),'linewidth',2);
xlabel('t(sec)');
ylabel('theta2(deg)');
subplot(2,2,3)
plot(t,rad2deg(y(:,3)),'linewidth',2);
xlabel('t(sec)');
ylabel('theta1_dot(deg/s)');
subplot(2,2,4)
plot(t,rad2deg(y(:,4)),'linewidth',2);
xlabel('t(sec)');
ylabel('theta2_dot(deg/s)');
figure(2)
subplot(2,1,1)

plot(t,(y))



function dz = ode_link(t,z)
m1=1;m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

dz= zeros(4,1);
z=num2cell(z);
[theta1,theta2,theta1_dot,theta2_dot] = deal(z{:});

if abs(theta1)>2*pi
    theta1= mod(theta1,2*pi);

end
if abs(theta2)>2*pi
    theta2= mod(theta2,2*pi);
end


K= [8,4,6,2;6,8,3,1.5];
M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2,m2*r2^2+l1*m2*cos(theta2)*r2+I2;I2+m2*r2^2+m2*r2*l1*cos(theta2),I2+m2*r2^2] ;
C= [0,l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m1*r2*theta2_dot*sin(theta2),0];
G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);-g*m2*r2*sin(theta1+theta2)];
q_d1 = pi- (3*pi*t^2)/100+(pi*t^3)/500;
qdot_d1= -(6*pi*t)/100+(3*t^2*pi)/500;
qddot_d1=-(6*pi)/100+(6*t*pi)/500;
q_d2 = pi/2- (3*pi*t^2)/200+(pi*t^3)/1000;
qdot_d2= -(6*pi*t)/200+(3*t^2*pi)/1000;
qddot_d2=-(6*pi)/200+(6*t*pi)/1000;

% U =double(subs( M*(-K* (([theta1;theta2;theta1_dot;theta2_dot]-[q_d1;qdot_d1;q_d2;qdot_d2])+[qddot_d1;qddot_d2])))+C*[theta1_dot;theta2_dot]+G));
% U= M*((-K*([theta1;theta2;theta1_dot;theta2_dot]-[q_d1;qdot_d1;q_d2;qdot_d2]))+[qddot_d1;qddot_d2])+C*[theta1_dot;theta2_dot]+G));
U = M*(-K* (  [theta1;theta2;theta1_dot;theta2_dot] - [q_d1;qdot_d1;q_d2;qdot_d2] ) +[qddot_d1;qddot_d2])+C*[theta1_dot;theta2_dot]+G;
u1=[U(1,:)];
u2=[U(2,:)];
dz(1)= theta1_dot;
dz(2)=theta2_dot;
dz(3)=(I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dz(4)=-(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);



end
