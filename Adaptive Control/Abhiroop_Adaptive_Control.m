clear
clc
clear all
syms a0 a1 a2 a3
syms to tf t
syms qo qodot qf qfdot

a=[a0;a1;a2;a3];
A = [1 to to^2 to^3;0 1 2*to 3*(to^2);1 tf (tf)^2 tf^3;0 1 2*tf 3*(tf^2)];
B= [qo;qodot;qf;qfdot];
a= inv(A)*B;

Cub1 = subs(a,[to,tf,qo,qodot,qf,qfdot],[0,10,pi,0,0,0]);
Cub2 = subs(a,[to,tf,qo,qodot,qf,qfdot],[0,10,pi/2,0,0,0]);
Cub1=double(Cub1);
Cub2=double(Cub2);
%Equation is of the form: a0 + a1*t^+a2*t^2+a3*t^3
%The 2 cubic equations are given by 
Cubicpoly1= Cub1(1)+Cub1(2)*t+Cub1(3)*(t^2)+Cub1(4)*(t^3);
Cubicpoly2= Cub2(1)+Cub2(2)*t+Cub2(3)*(t^2)+Cub2(4)*(t^3);
q_d1 = pi- (3*pi*t^2)/100+(pi*t^3)/500;
qdot_d1= -(6*pi*t)/100+(3*t^2*pi)/500;
qddot_d1=-(6*pi)/100+(6*t*pi)/500;
q_d2 = pi/2- (3*pi*t^2)/200+(pi*t^3)/1000;
qdot_d2= -(6*pi*t)/200+(3*t^2*pi)/1000;
qddot_d2=-(6*pi)/200+(6*t*pi)/1000;


%Manipulator Equation given by
% % M(q) ̈q + C(q,  ̇q)  ̇q + g(q) = τ
syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 
syms theta1_dot theta2_dot theta1_ddot theta2_ddot u1 u2 
g=9.81
eq1= theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
eq2= theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);
eq1=simplify(eq1);
eq2=simplify(eq2);

%Manipulator Equations
M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2,m2*r2^2+l1*m2*cos(theta2)*r2+I2;I2+m2*r2^2+m2*r2*l1*cos(theta2),I2+m2*r2^2] ;
C= [0,l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m1*r2*theta2_dot*sin(theta2),0];
G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);-g*m2*r2*sin(theta1+theta2)];

%State feedback control:
%Acl=A-BK
AK=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
BK=[0,0;0,0;1,0;0,1];
Eigen=[-3,-3,-4,-4];
K=place(AK,BK,Eigen);
%From K
Kd=[12,0;0,12];
Kp=[7,0;0,7];

Acl=[0,0,1,0;0,0,0,1;Kd,Kp];
Q=eye(4);
P= lyap(Acl',Q);
T=10;
% T=10;
[t,y] = ode45(@ode_link, [0,T], [(10*pi)/9,(25*pi)/36,0,0,1.1798,0.3375,0.2149,1.0875,0.3375]);
theta1=zeros(length(t),1);
theta2=zeros(length(t),1);
theta1_dot=zeros(length(t),1);
theta2_dot=zeros(length(t),1);
theta1_desired =zeros(length(t),1);
theta2_desired =zeros(length(t),1);
theta1_dot_desired =zeros(length(t),1);
theta2_dot_desired =zeros(length(t),1);
Tau1=zeros(length(t),1);
Tau2=zeros(length(t),1);
qd=[q_d1;q_d2];
qd_dot=[qdot_d1;qdot_d2];
qd_ddot=[qddot_d1;qddot_d2];
for i=1:length(t)
    qd_sub = subs(qd,t(i));
    qd_dot_sub = subs(qd_dot,t(i));
    qd_ddot_sub = subs(qd_ddot,t(i));
    theta1_desired(i) = qd_sub(1);
    theta2_desired(i) = qd_sub(2);
    theta1_dot_desired(i) = qd_dot_sub(1);
    theta2_dot_desired(i) = qd_dot_sub(2);  
    
    j=[y(i,1) ; y(i,2)];
    j_d=[theta1_desired(i);theta2_desired(i)];
    jd=[y(i,3) ; y(i,4)];
    jd_d=[theta1_dot_desired(i);theta2_dot_desired(i)];
    x=[j-j_d;jd-jd_d]; %j = theta
    

    v=-K*x+qd_ddot;
    Mmat= [y(i,5)+2*y(i,6)*cos(y(i,2)), y(i,7)+y(i,6)*cos(y(i,2)); y(i,7)+y(i,6)*cos(y(i,2)), y(i,7)];

    Cmat= [-y(i,6)*sin(y(i,2))*y(i,4), -y(i,6)*sin(y(i,2)+y(i,1))*(y(i,4)+y(i,3)); y(i,6)*sin(y(i,2)+y(i,1))*y(i,3),0];
    Gmat= [-y(i,8)*g*sin(y(i,1))-y(i,9)*g*sin(y(i,2)+y(i,1)); -y(i,9)*g*sin(y(i,2)+y(i,1))];

    tau=Mmat*v+Cmat*[y(i,3);y(i,4)]+Gmat;
    tau_sub=subs(tau,t(i));
    theta1(i)=(y(i,1));
    theta2(i)=(y(i,2));
    theta1_dot(i)=(y(i,3));
    theta2_dot(i)=(y(i,4));
    
    
    tau_sub1=subs(tau,[theta1(i),theta2(i),theta1_dot(i),theta2_dot(i)]);
    Tau1(i)=tau_sub1(1);
    Tau2(i)=tau_sub1(2);
    
    
    
    
    
    
    
end





figure(1)
subplot(2,2,1)
plot(t,y(:,1),'linewidth',2);
hold on
plot(t,theta1_desired,'--','linewidth',2);
title('Theta 1 vs Theta 1 desired')
xlabel('Time in (sec)');
ylabel('Angle in (rad)');
legend('Theta1','Theta1 Desired')
hold off
subplot(2,2,2)
plot(t,y(:,2),'linewidth',2);
hold on
plot(t,theta2_desired,'--','linewidth',2);
title('Theta 2 vs Theta 2 desired')
xlabel('Time in (sec)');
ylabel('Angle in (rad)');
legend('Theta2','Theta2 Desired')
hold off
subplot(2,2,3)
plot(t,y(:,3),'linewidth',2);
hold on
plot(t,theta1_dot_desired,'--','linewidth',2);
title('Theta 1 dot vs Theta 1 dot desired')
xlabel('Time in (sec)');
ylabel('Angle in (rad)');
legend('Theta1 dot','Theta1 dot Desired')
hold off
subplot(2,2,4)
plot(t,y(:,4),'linewidth',2);
hold on
plot(t,theta2_dot_desired,'--','linewidth',2);
title('Theta 2 dot vs Theta 2 dot desired')
xlabel('Time in (sec)');
ylabel('Angle in (rad)');
legend('Theta2 dot','Theta2 dot Desired','Location','best')
hold off
figure(2)
subplot(2,2,1)       
plot(t,Tau1,'linewidth',2); 
title('Time vs Torque at joint-1')
xlabel("Time in Seconds")
ylabel("Torque in N-m")


subplot(2,2,2)       
plot(t,Tau2,'linewidth',2); 
title('Time vs Torque at joint-2')
xlabel("Time in Seconds")
ylabel("Torque in N-m")

figure(3)
subplot(2,1,1)
n=[y(:,5),y(:,6),y(:,7),y(:,8),y(:,9)];
plot(t,n,'linewidth',2)
legend('Alpha1','Alpha2','Alpha3','Alpha4','Alpha5')


function dz = ode_link(t,z)
m1=1;m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;
dz= zeros(9,1);
z=num2cell(z);
[theta1,theta2,theta1_dot,theta2_dot,a,b,c,d,f] = deal(z{:});

if abs(theta1)>2*pi
    theta1= mod(theta1,2*pi);

end
if abs(theta2)>2*pi
    theta2= mod(theta2,2*pi);
end
a_hat=cell2mat(z(5));
b_hat=cell2mat(z(6));
c_hat=cell2mat(z(7));
d_hat=cell2mat(z(8));
f_hat=cell2mat(z(9));

M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2,m2*r2^2+l1*m2*cos(theta2)*r2+I2;I2+m2*r2^2+m2*r2*l1*cos(theta2),I2+m2*r2^2] ;
C= [0,l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m1*r2*theta2_dot*sin(theta2),0];
G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);-g*m2*r2*sin(theta1+theta2)];
q_d1 = pi- (3*pi*t^2)/100+(pi*t^3)/500;
qdot_d1= -(6*pi*t)/100+(3*t^2*pi)/500;
qddot_d1=-(6*pi)/100+(6*t*pi)/500;
q_d2 = pi/2- (3*pi*t^2)/200+(pi*t^3)/1000;
qdot_d2= -(6*pi*t)/200+(3*t^2*pi)/1000;
qddot_d2=-(6*pi)/200+(6*t*pi)/1000;
% I1mat=0.063;
% I2mat=0.063;
% m1mat=0.75;
% m2mat=0.75;
% a = I1mat + I2mat + m1mat*r1^2 + m2mat*(l1^2 + r2^2);
% b = m2mat*l1*r2;
% d = I2mat + m2mat*r2^2;
% Mmat= [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
% Cmat= [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
% Gmat= [-m1*g*r1*sin(theta1)-m2*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2*g*r2*sin(theta1+theta2)];
% Mmat= [cell2mat(z(5))+2*cell2mat(z(6))*cos(theta2), cell2mat(z(8))+cell2mat(z(6))*cos(theta2); cell2mat(z(8))+cell2mat(z(6))*cos(theta2), cell2mat(z(8))];
% Cmat= [-cell2mat(z(7))*sin(theta2)*theta2_dot, -cell2mat(z(7))*sin(theta2)*(theta1_dot+theta2_dot); cell2mat(z(7))*sin(theta2)*theta1_dot,0];
% Gmat= [-cell2mat(z(7))*g*cell2mat(z(9))*sin(theta1)-cell2mat(z(9))*g*sin(theta1+theta2); -cell2mat(z(9))*g*sin(theta1+theta2)];
Mmat= [a_hat+2*b_hat*cos(theta2), c_hat+b_hat*cos(theta2); c_hat+b_hat*cos(theta2), c_hat];

Cmat= [-b_hat*sin(theta2)*theta2_dot, -b_hat*sin(theta2+theta1)*(theta2_dot+theta1_dot); b_hat*sin(theta2+theta1)*theta1_dot,0];
Gmat= [-d_hat*g*sin(theta1)-f_hat*g*sin(theta2+theta1); -f_hat*g*sin(theta2+theta1)];


theta=[theta1;theta2];
theta_d=[q_d1;q_d2];
theta_dot=[theta1_dot;theta2_dot];
theta_ddot=[qdot_d1;qdot_d2];
thetaddot= [qddot_d1;qddot_d2];
e=[theta-theta_d];
e_dot=[theta_dot-theta_ddot];

AK=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
BK=[0,0;0,0;1,0;0,1];
Eigen=[-3,-3,-4,-4];
K=place(AK,BK,Eigen);
Acl=[0,0,1,0;0,0,0,1;-K];
% Q=eye(4)*2.1
Q=eye(4);
P= lyap(Acl',Q);
BK=[0,0;0,0;1,0;0,1];
x=[theta-theta_d;theta_dot-theta_ddot];



v=-K*x+thetaddot;
tau=Mmat*v+Cmat*[theta_dot]+Gmat;
ddq=inv(M)*(tau-C*cell2mat(z(1:2))-G);
thh1=ddq(1);
thh2=ddq(2);

% U = M*(-K* (  [theta1;theta2;theta1_dot;theta2_dot] - [q_d1;qdot_d1;q_d2;qdot_d2] ) +[qddot_d1;qddot_d2])+C*[theta1_dot;theta2_dot]+G;
u1=[tau(1,:)];
u2=[tau(2,:)];
disp(u1)
disp(u2)
theta1_ddot=(I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

theta2_ddot=-(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dz(1,1)= theta1_dot;
dz(2,1)=theta2_dot;
dz(3)=theta1_ddot;
dz(4)=theta2_ddot;

% dz(3,1)=th1;
% % dz(4,1)=th2;
% thh1=cell2mat(z(3));
% thh2=cell2mat(z(4));



Y = [thh1, ...
cos(theta2)*(2*thh1 + thh2) - 2*sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
thh2, ...
-sin(theta1)*g, ...
-sin(theta1 + theta2)*g; ...
0, ...
sin(theta2)*theta1_dot^2 + cos(theta2)*thh1, ...
thh1 + thh2, ...
0, ...
-sin(theta1+theta2)*g];
Phi= Mmat \Y;
% Gamma=eye(5)*90;
Gamma=eye(5)*0.5;
dalpha_est= -Gamma\(Phi'*BK'*P*x);

dz(5,1)=dalpha_est(1);
dz(6,1)=dalpha_est(2);
dz(7,1)=dalpha_est(3);
dz(8,1)=dalpha_est(4);
dz(9,1)=dalpha_est(5);

end