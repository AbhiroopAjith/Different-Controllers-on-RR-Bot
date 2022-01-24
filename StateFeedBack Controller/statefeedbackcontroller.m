clear
clc

syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 
syms theta1_dot theta2_dot theta1_ddot theta2_ddot u1 u2 g
m1=1; 
m2= 1;
l1=1;
l2= 1; 
r1=0.45;
r2= 0.45;
I1=0.084;
I2= 0.084;
g= 9.81;
eq1= theta1_ddot*(m2*l1^2 + 2*m2*cos(theta2)*l1*r2 + m1*r1^2 + m2*r2^2 + I1 + I2) - theta2_dot*(l1*m2*r2*theta1_dot*sin(theta2) + l1*m2*r2*sin(theta2)*(theta1_dot + theta2_dot)) - u1 + theta2_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) - g*l1*m1*sin(theta1) - g*m1*r1*sin(theta1);
eq2= theta2_ddot*(m2*r2^2 + I2) - u2 + theta1_ddot*(m2*r2^2 + l1*m2*cos(theta2)*r2 + I2) - g*m2*r2*sin(theta1 + theta2) + l1*m2*r2*theta1_dot*sin(theta2)*(theta1_dot + theta2_dot) - l1*m2*r2*theta1_dot*theta2_dot*sin(theta2);

eq1 = subs(eq1,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]);
eq2 = subs(eq2,[theta1_dot,theta2_dot,theta1_ddot,theta2_ddot,u1,u2],[0,0,0,0,0,0]);

sol=solve([eq1==0 , eq2==0] ,[theta1,theta2]);
sol.theta1
sol.theta2


%state space representation
u=[u1;u2];
x=[theta1,theta2,theta1_dot,theta2_dot];
x1_dot = theta1_dot;
x2_dot= theta2_dot;
x3_dot= (I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

x4_dot= -(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

X_dot = [x1_dot;x2_dot;x3_dot;x4_dot];

A = jacobian(X_dot,x);
B=jacobian(X_dot,u);

A= subs(A,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);
B= subs(B,[theta1,theta2,theta1_dot,theta2_dot],[0,0,0,0]);

A= double(A);
B= double(B);


%Stability
eiga= eig(A) %Unstable

%Stability for the other Equlibrium Points. There are a total of 4
%equlibrium points. We already checked for the stablity of the (0,0,0,0)
%Point
%0,0,0,0
%pi,0,0,0
%0,pi,0,0
%pi,pi,0,0

A1 = jacobian(X_dot,x);
A1= subs(A1,[theta1,theta2,theta1_dot,theta2_dot],[pi,0,0,0]);
A1= double(A1);

eiga1=eig(A1)  %stable at pi,0,0,0

A2 = jacobian(X_dot,x);
A2= subs(A2,[theta1,theta2,theta1_dot,theta2_dot],[0,pi,0,0]);
A2= double(A2);

eiga2=eig(A2) %unstable at 0,pi,0,0


A3 = jacobian(X_dot,x);
A3= subs(A3,[theta1,theta2,theta1_dot,theta2_dot],[pi,pi,0,0]);
A3= double(A3);

eiga3=eig(A3) %unstable pi,pi,0,0



%Check for controllability 
rankC=rank(ctrb(A,B)); %controllable because it is fullrank


syms k11 k12 k13 k14 k21 k22 k23 k24 lambda
%State-Feedback
k=[k11,k12,k13,k14;k21,k22,k23,k24];
Acl = A -(B*k);
AclPoly = simplify(det(Acl-lambda*eye(4)));

%finding k values
lamda = [-1,-2,-1+1i,-1-1i];
A=double(A);
B=double(B);
KValues= place(A,B,lamda)



T=10;
[t,y] = ode45(@ode_link, [0,T], [pi/6,pi/4,0,0]);
Kv= [ 23.5850,5.8875,5.1470,2.6104; 5.8875,4.9875,1.5543,0.9970];
for i = 1:size(y,1)
    u(i)= -( Kv (1,1)*y(i,1) + Kv(1,2)*y(i,2)+ Kv (1,3)*y(i,3) + Kv(1,4)*y(i,4));
    g(i)= -( Kv (2,1)*y(i,1) + Kv(2,2)*y(i,2)+ Kv (2,3)*y(i,3) + Kv(2,4)*y(i,4));
   
end
%Plotting    
 
%overall graph
figure(2)
plot(t,y);

%Seperate graphs
figure(1)
subplot(2,2,1);
plot(t,rad2deg(y(:,1)));
xlabel('t(sec)');
ylabel('theta1(deg)');

subplot(2,2,2);
plot(t,rad2deg(y(:,2)));
xlabel('t(sec)');
ylabel('theta2(deg)');

subplot(2,2,3);
plot(t,rad2deg(y(:,3)));
xlabel('t(sec)');
ylabel('theta1dot(deg/s)');

subplot(2,2,4);
plot(t,rad2deg(y(:,4)));
xlabel('t(sec)');
ylabel('theta2dot(deg/s)');




%Ode Function
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
Kv= [ 23.5850,5.8875,5.1470,2.6104; 5.8875,4.9875,1.5543,0.9970];
U = -Kv* [theta1;theta2;theta1_dot;theta2_dot];


u1=[U(1,:)];
u2=[U(2,:)];
dz(1)= theta1_dot;
dz(2)=theta2_dot;
dz(3)=(I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dz(4)=-(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);



end







