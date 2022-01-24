clear; close; clc;
% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');

JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;


send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(30), deg2rad(45)];
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

U = M*(-K* (  [theta1;theta2;theta1_dot;theta2_dot] - [q_d1;qdot_d1;q_d2;qdot_d2] ) +[qddot_d1;qddot_d2])+C*[theta1_dot;theta2_dot]+G;
    tau1.Data = U(1);
   tau2.Data = U(2);
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
figure(1)
subplot(2,2,1);
plot(time,rad2deg(g1));
xlabel('t(sec)');
ylabel('theta1(deg)');

subplot(2,2,2);
plot(time,rad2deg(g2));
xlabel('t(sec)');
ylabel('theta2(deg)');

subplot(2,2,3);
plot(time,rad2deg(g3));
xlabel('t(sec)');
ylabel('theta1dot(deg/s)');

subplot(2,2,4);
plot(time,rad2deg(g4));
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



%disconnect from roscore
rosshutdown;
