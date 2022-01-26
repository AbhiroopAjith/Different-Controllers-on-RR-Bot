clear; close; clc;
% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');

JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
joints=rosmessage(JointStates);
tau1.Data = 0;
tau2.Data = 0;


send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;
counter=1;

while(t < 10)
    t = toc;
 m1=1;m2=1; l1=1; l2=1 ;r1=0.45; r2=0.45; g=9.81 ;I2= 0.084; I1= 0.084;

    % read the joint states
    jointData = receive(JointStates);
   x = [jointData.Position(1);jointData.Position(2);jointData.Velocity(1)
        jointData.Velocity(2)];
    theta1=jointData.Position(1);
    theta2=jointData.Position(2);
    theta1_dot=jointData.Velocity(1);
    theta2_dot=jointData.Velocity(2);
    
    if abs(theta1) > 2*pi
        theta1=mod(theta1,2*pi);
    end
    if abs(theta2) > 2*pi
        theta1=mod(theta2,2*pi);
    end
    
%   M = [m2*l1^2+2*m2*cos(theta2)*l2*r2+m1*r1^2+m2*r2^2+I1+I2,m2*r2^2+l1*m2*cos(theta2)*r2+I2;I2+m2*r2^2+m2*r2*l1*cos(theta2),I2+m2*r2^2] ;
% C= [0,l1*m2*r2*theta1_dot*sin(theta2)+l1*m2*r2*sin(theta2)*(theta1_dot+theta2_dot);l1*m2*r2*(theta1_dot+theta2_dot)*sin(theta2)-l1*m1*r2*theta2_dot*sin(theta2),0];
% G = [-g*l1*m2*sin(theta1)-g*m1*r1*sin(theta1)-m2*g*r2*sin(theta1+theta2);-g*m2*r2*sin(theta1+theta2)];
% q_d1 = pi- (3*pi*t^2)/100+(pi*t^3)/500;
% qdot_d1= -(6*pi*t)/100+(3*t^2*pi)/500;
% qddot_d1=-(6*pi)/100+(6*t*pi)/500;
% q_d2 = pi/2- (3*pi*t^2)/2+(pi*t^3)/1000;
% qdot_d2= -(6*pi*t)/200+(3*t^2*pi)/1000;
% qddot_d2=-(6*pi)/200+(6*t*pi)/1000;
qd= [0.0063*t^3 - 0.0942*t^2 + 3.1416;
        0.0031*t^3 - 0.0471*t^2 + 1.5708];
    qd_dot =    [0.0188*t^2 - 0.1885*t;
                0.0094*t^2 - 0.0942*t];
    qd_ddot =   [0.0377*t - 0.1885;
                0.0188*t - 0.0942];
    

I1mat=0.063;
I2mat=0.063;
m1mat=0.75;
m2mat=0.75;
a = I1mat + I2mat + m1mat*r1^2 + m2mat*(l1^2 + r2^2);
b = m2mat*l1*r2;
d = I2mat + m2mat*r2^2;
Mmat= [a+2*b*cos(theta2), d+b*cos(theta2); d+b*cos(theta2), d];
Cmat= [-b*sin(theta2)*theta2_dot, -b*sin(theta2)*(theta1_dot+theta2_dot); b*sin(theta2)*theta1_dot,0];
Gmat= [-m1mat*g*r1'*sin(theta1)-m2mat*g*(l1*sin(theta1)+r2*sin(theta1+theta2)); -m2mat*g*r2*sin(theta1+theta2)];
theta=[theta1;theta2];
% theta_d=[q_d1;q_d2];
theta_d=[qd];
theta_dot=[theta1_dot;theta2_dot];
% theta_ddot=[qdot_d1;qdot_d2];
theta_ddot=[qd_dot];
% thetaddot= [qddot_d1;qddot_d2];
thetaddot=[qd_ddot];
rho=4.8;
Phi=0.015;
Kp=[12,0;0,12];
Kd=[7,0;0,7];
AK=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
BK=[0,0;0,0;1,0;0,1];
Eigen=[-3,-3,-4,-4];
K=place(AK,BK,Eigen);
Acl=[0,0,1,0;0,0,0,1;-K];
Q=eye(4)*1.089;
P= lyap(Acl',Q);
BK=[0,0;0,0;1,0;0,1];
x=[theta-theta_d;theta_dot-theta_ddot];
if Phi>0
    if norm(x'*P*BK)>Phi
        vr= -rho*(x'*P*BK)/(norm(x'*P*BK));
    else
        vr= -rho*(x'*P*BK)/(Phi);
    end
else
    if (norm(x'*P*BK))~=0
         vr=-rho*(x'*P*BK)/(norm(x'*P*BK));
    else
        vr=0;
    end
end

v=-K*x+thetaddot+vr';
tau=Mmat*v+Cmat*[theta_dot]+Gmat;

   tau1.Data = tau(1);
   tau2.Data = tau(2);
   send(j1_effort,tau1);
   send(j2_effort,tau2);
   
   g1(counter)=jointData.Position(1);
   g2(counter)=jointData.Position(2);
   g3(counter)=jointData.Velocity(1);
   g4(counter)=jointData.Velocity(2);
   force1(counter)=tau1.Data;
    
   force2(counter)=tau2.Data;  
   
   time(counter)=t;

   counter=counter+1;
   
   
   
end
tau1.Data=0;
tau2.Data=0;
send(j1_effort,tau1);
send(j2_effort,tau2);
rosshutdown;
theta1_desired =zeros(length(time),1);
theta2_desired =zeros(length(time),1);
theta1_dot_desired =zeros(length(time),1);
theta2_dot_desired =zeros(length(time),1);



for i=1:length(time)
    t= time(i);
    qd = [(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
      (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2];
    qd_dot = [(3*pi*t^2)/500 - (3*pi*t)/50;
          (3*pi*t^2)/1000 - (3*pi*t)/100];
    theta1_desired(i)=qd(1);

    theta2_desired(i) = qd(2);
    theta1_dot_desired(i) = qd_dot(1);
    theta2_dot_desired(i) = qd_dot(2); 
end

figure(1)
subplot(2,2,1);
plot(time,g1);
hold on
plot(time,theta1_desired);
hold off
xlabel('t(sec)');
ylabel('theta1(deg)');

subplot(2,2,2);
plot(time,g2);
hold on
plot(time,theta2_desired);
hold off

xlabel('t(sec)');
ylabel('theta2(deg)');

subplot(2,2,3);
plot(time,g3);
hold on
plot(time,theta1_dot_desired);
hold off

xlabel('t(sec)');
ylabel('theta1dot(deg/s)');

subplot(2,2,4);
plot(time,g4);
hold on
plot(time,theta2_dot_desired);
hold off
xlabel('t(sec)');
ylabel('theta2dot(deg/s)');

figure(2)
subplot(2,2,1);
plot(time,force1);
xlabel('t(sec)');
ylabel('Force1(N)');

figure(2)
subplot(2,2,2);
plot(time,force2);
xlabel('t(sec)');
ylabel('Force2(N)');



