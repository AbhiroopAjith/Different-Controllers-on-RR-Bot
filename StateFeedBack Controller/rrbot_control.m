clear; close; clc;
% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');

JointStates = rossubscriber('/rrbot/joint_states')

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
     
    % read the joint states
    jointData = receive(JointStates);
    x = [jointData.Position(1);jointData.Position(2);jointData.Velocity(1)
        jointData.Velocity(2)];
    k1=[23.5850,5.8875,5.1470,2.6104];
    k2= [5.8875, 4.9875,1.5543,0.9970];
   
   tau1.Data = -k1*x;
   tau2.Data = -k2*x;
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
ylabel('theta2(deg/s)');

subplot(2,2,3);
plot(time,rad2deg(g3));
xlabel('t(sec)');
ylabel('theta1dot(deg)');

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
