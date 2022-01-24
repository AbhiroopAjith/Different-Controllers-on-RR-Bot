syms m1 m2 theta1 theta2 r1 r2 l1 l2 I1 I2 theta1_dot theta2_dot theta1_ddot theta2_ddot u1 u2 g

Kinetic1 = 0.5*m1*(theta1_dot*r1)^2 + 0.5*I1*theta1_dot^2;

Potential1 = m1*g*r1*cos(theta1);

Kinetic2 = 0.5*m2*(theta1_dot*l1)^2+0.5*m2*(r2*(theta1_dot+theta2_dot))^2 + m2*l1*theta1_dot*r2*(theta1_dot+theta2_dot)*cos(theta2)+0.5*I2*(theta1_dot+theta2_dot)^2;

Potential2= m1*g*l1*cos(theta1)+ m2*g*r2*cos(theta1+theta2);

L = (Kinetic1+Kinetic2)-(Potential1+Potential2);

partialdiff= jacobian(L,[theta1,theta2,theta1_dot,theta2_dot]);

partialdt1 = jacobian(partialdiff(3),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];

partialdt2= jacobian(partialdiff(4),[theta1;theta2;theta1_dot;theta2_dot])*[theta1_dot;theta2_dot;theta1_ddot;theta2_ddot];

eq1= partialdt1-partialdiff(1)-u1;
eq2=partialdt2-partialdiff(2)-u2;

sol=solve([eq1==0 , eq2==0] ,[theta1_ddot,theta2_ddot])
sol.theta1_ddot
sol.theta2_ddot

X = sym('X',[4,1]);
X(1)=theta1
X(2)=theta2
X(3)=theta1_dot
X(4)=theta2_dot

[t,y] = ode45(@ode_link, [0,10], [pi/6,pi/4,0,0]);
plot(t,y);


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
u1=0;u2=0;
dz(1)= theta1_dot;
dz(2)=theta2_dot;
dz(3)=(I2*u1 - I2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

dz(4)=-(I2*u1 - I1*u2 - I2*u2 - l1^2*m2*u2 - m1*r1^2*u2 + m2*r2^2*u1 - m2*r2^2*u2 + l1*m2^2*r2^3*theta1_dot^2*sin(theta2) + l1^3*m2^2*r2*theta1_dot^2*sin(theta2) + l1*m2^2*r2^3*theta2_dot^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + I2*g*l1*m1*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2*r2*u1*cos(theta2) - 2*l1*m2*r2*u2*cos(theta2) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta1_dot^2*sin(theta2) + I2*l1*m2*r2*theta2_dot^2*sin(theta2) + g*l1*m1*m2*r2^2*sin(theta1) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + g*l1^2*m1*m2*r2*cos(theta2)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(- l1^2*m2^2*r2^2*cos(theta2)^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);



end