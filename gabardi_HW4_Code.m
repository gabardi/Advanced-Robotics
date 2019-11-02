%% ME 739 Homework #4
%% Kaitlyn Gabardi
%04/17/2019

%% Problem 1 
samples = 1000;
% Write a Matlab function which constructs a quintic polynomial

% Input parameters
t0 = 0;        %initial time
tf = 5;        %final time
y0 = -10;      %initial position
yf = 10;       %final position
yd0 = -50;     %initial velocity
ydf = -50;     %final velocity
ydd0 = 0;      %initial acceleration
yddf = 0;      %final acceleration 

%Verification function is correct
[t, y, dy, ddy] = QuinticPolynomial(0,5,-10,10,-50,-50,0,0,samples);

figure(1); 
plot(t,y,'r');
xlabel('Time')
ylabel('Position')
 
figure(2); 
plot(t,dy,'g');
xlabel('Time')
ylabel('Velocity')
 
figure(3);
plot(t,ddy,'b'); 
xlabel('Time')
ylabel('Acceleration')


%% Problem 2

%% Scenario 1: Joint-space trajectory generation
%Time between each successive waypoint (5 sec)

%total number of samples
samples = 1000;

% Define joint-space displacements, v and a which was found using inverse
% kinematics
L = 4; %link length

% joints displacement for the desired steps (1>2>3>4>1). 
% rows (3 total) for each joint

q = [  2.5*L        2.5*L        1.5*L        1.5*L         2.5*L
      3*pi/4         pi/4         pi/4       3*pi/4        3*pi/4
   sqrt(8)*L    sqrt(8)*L    sqrt(8)*L    sqrt(8)*L     sqrt(8)*L];

%joint-space velocity (defined from the problem statement)
dq0 = 0; 
dqf = 0;   

%joint-space acceleration (defined from the problem statement)
ddq = 0; 
ddqf = 0;  

%Defined time vector required to move between successive waypoints
time = linspace(0,20,5);  

% Loop to iterate through the # of joints and trajectories
for i = 1:3 %total joints
    figure();
    for j = 1:4 %total paths
        % generate trajectory displacement, velocity, and acceleration
        
        % q0 = initial displacement and qf = final displacements
        q0 = q(i,j); 
        qf = q(i,j+1); 
        
        % Defines the initial time to final time for each trajectory
        t0 = time(j); 
        tf = time(j+1);  
        
        % QuinticPolynomial function will store above variables in
        [t,y,dy,ddy] = QuinticPolynomial(t0,tf,q0,qf,dq0,dqf,dq0,dqf,samples);
        
        time_plot(j,:,i) = t;  %time
        q_plot(j,:,i) = y;     %position
        dq_plot(j,:,i) = dy;   %velocity
        ddq_plot(j,:,i) = ddy; %acceleraton
        
        subplot(3,1,1); 
        title('Joint Space Displacements as a Fxn of Time');
        grid on; 
        hold on;
        plot(time_plot(j,:,i),q_plot(j,:,i),'b');
        xlabel('Time'); 
        ylabel('Position');
        
        % velocity plot
        subplot(3,1,2); 
        grid on; 
        hold on
        plot(time_plot(j,:,i),dq_plot(j,:,i),'r');
        xlabel('Time'); 
        ylabel('Velocity');
        
        % acceleration plot
        subplot(3,1,3); 
        grid on; 
        hold on;
        plot(time_plot(j,:,i), ddq_plot(j,:,i),'g');
        xlabel('Time'); 
        ylabel('Acceleration');  
        
    end
end


% Using forward kinematics calculate task space displacements
xe = q_plot(:,:,3).*cos(q_plot(:,:,2));
ye = q_plot(:,:,3).*sin(q_plot(:,:,2));
ze = q_plot(:,:,1);

% 2D and 3D plots of task-space
% Bassim Younis credit: Need to define variable names for figures,
% otherwise will generate seperate plots within for loop
f2 = figure;  %figure for 2D plot
f3 = figure;  %figure for 3D plot

% Task-space displacments as a fxn of time
for k = 1:4
     figure(f2);
     title('Task-space Displacements as a Fxn of Time');
    
    % Task-space xe plot
    subplot(3,1,1); 
    grid on; 
    hold on;
    plot(time_plot(k,:,1),xe(k,:),'b');
    xlabel('Time'); 
    ylabel('xe');
    
    % Task-space ye plot
    subplot(3,1,2); 
    grid on; 
    hold on
    plot(time_plot(k,:,1),ye(k,:),'r');
    xlabel('Time'); 
    ylabel('ye');
    
    %Task-space ze plot
    subplot(3,1,3); 
    grid on; 
    hold on;
    plot(time_plot(k,:,1),ze(k,:),'g'); 
    xlabel('Time'); 
    ylabel('ze');
    
    % Use 'plot3' to create task space displacements in 3-D space
    figure(f3);  
    plot3(xe(k,:),ye(k,:),ze(k,:),'b');
    hold on;
end
figure(f3)
xlabel('xe'); 
ylabel('ye'); 
zlabel('ze'); 
grid on;
title('Task Space Displacements for Scenario #1');

%% Scenario 2: Task-space trajectory generation

close all;

%total number of samples
samples = 1000;

% Define joint-space displacements, v and a which was found using inverse
% kinematics
L = 4; %link length

fwd_k = [-2*L     2*L     2*L     -2*L    -2*L
          2*L     2*L     2*L      2*L     2*L
         2.5*L   2.5*L   1.5*L    1.5*L   2.5*L];

%joint-space velocity (defined from the problem statement)
dq0 = 0; 
dqf = 0;   

%joint-space acceleration (defined from the problem statement)
ddq = 0; 
ddqf = 0;  

%Defined time vector required to move between successive waypoints
time = linspace(0,20,5);  

% Loop to iterate through the # of joints and trajectories
for i = 1:3 %total joints
    figure();
    for j = 1:4 %total paths
        % generate trajectory displacement, velocity, and acceleration
        
        % q0 = initial displacement and qf = final displacements
        q0_2 = fwd_k(i,j); 
        qf_2 = fwd_k(i,j+1); 
        
        % Defines the initial time to final time for each trajectory
        t0_2 = time(j); 
        tf_2 = time(j+1);  
        
        % QuinticPolynomial function will store above variables in
        [t,y,dy,ddy] = QuinticPolynomial(t0_2,tf_2,q0_2,qf_2,dq0,dqf,dq0,dqf,samples);
        
        time_plot_2(j,:,i) = t;  %time
        q_plot_2(j,:,i) = y;     %position
        dq_plot_2(j,:,i) = dy;   %velocity
        ddq_plot_2(j,:,i) = ddy; %acceleraton
        
        subplot(3,1,1); 
        title('Joint Space Displacements as a Fxn of Time');
        grid on; 
        hold on;
        plot(time_plot_2(j,:,i),q_plot_2(j,:,i),'b');
        xlabel('Time'); 
        ylabel('Position');
        
        % velocity plot
        subplot(3,1,2); 
        grid on; 
        hold on
        plot(time_plot_2(j,:,i),dq_plot_2(j,:,i),'r');
        xlabel('Time'); 
        ylabel('Velocity');
        
        % acceleration plot
        subplot(3,1,3); 
        grid on; 
        hold on;
        plot(time_plot_2(j,:,i), ddq_plot_2(j,:,i),'g');
        xlabel('Time'); 
        ylabel('Acceleration');  
        
    end
end

% Using forward kinematics calculate task space displacements

xe_2 = q_plot_2(:,:,1);
ye_2 = q_plot_2(:,:,2);
ze_2 = q_plot_2(:,:,3);

% Inverse kinematrics 
q1_2 = ze_2;
q2_2 = atan2(ye_2,xe_2);
q3_2 = sqrt((xe_2).^2 + (ye_2).^2);

% 2D and 3D plots of joint-space
f3 = figure;  %figure for 2D plot of problem 2, part 2
f4 = figure;  %figure for 3D plot of problem 2, part 2

% Task-space displacments as a fxn of time
for k = 1:4
     figure(f3);
     title('Task-space Displacements as a Fxn of Time');
    
    % Task-space xe plot
    subplot(3,1,1); 
    grid on; 
    hold on;
    plot(time_plot_2(k,:,1),q1_2(k,:),'b');
    xlabel('Time'); 
    ylabel('q1 displacement');
    
    % Task-space ye plot
    subplot(3,1,2); 
    grid on; 
    hold on
    plot(time_plot_2(k,:,1),q2_2(k,:),'r');
    xlabel('Time'); 
    ylabel('q2 displacement');
    
    %Task-space ze plot
    subplot(3,1,3); 
    grid on; 
    hold on;
    plot(time_plot_2(k,:,1),q3_2(k,:),'g'); 
    xlabel('Time'); 
    ylabel('q3 displacement');
    
    % Use 'plot3' to create task space displacements in 3-D space
    figure(f4);  
    plot3(xe_2(k,:),ye_2(k,:),ze_2(k,:),'b');
    hold on;
end

figure(f4)
xlabel('xe'); 
ylabel('ye'); 
zlabel('ze'); 
grid on;
title('Task Space Displacements for Scenario #2');

%% Problem 3: 
%-----------------------------------------------------------------------
%  INITIALIZATION
%-----------------------------------------------------------------------
clear all; 
close all; clc

% input parameters
L = 4;                       %manipulator link lengths
alpha = 0.1;                 %gradient descent step size
epsilon = 0.001;             %gradient descent error tolerance
zeta = 100;                  %attractive force scaling

q0 = [L 0 L]';               %initial manipulator joint position
qf = [2*L pi/2 2*L]';        %goal position (joint space)

maxIterations = 1000;        %max gradient descent iterations
q = q0;                      %initialize joint positions
k = 1;                       %initialize counter

%goal position (task space)which was found (see analytical work attatched)
o1g = [  0
       2*L
       2*L];

%-----------------------------------------------------------------------
%  GRADIENT DESCENT
%-----------------------------------------------------------------------
while((q - qf)'*(q - qf) > epsilon)
    %trigometric expresions
    s2 = sin(q(2)); 
    c2 = cos(q(2)); 
   
    %forward kinematics
    o1 = [q(3)*c2
          q(3)*s2
            q(1)];
      
    %attractive forces
    FA1 = -zeta*(o1 - o1g);
    
    %point Jacobians
    Jo1 = [0  -q(3)*s2   c2
           0   q(3)*c2   s2
           1      0       0];
    
    %joint torques
    tau = Jo1'*FA1;
    
    %update q
    q = q + alpha*tau/sqrt(tau'*tau);
    Q(:,k) = q;

    %check convergence
    k = k + 1;
    if k > maxIterations
        break
    end
end

%-----------------------------------------------------------------------
%  RENDERING INITIALIZATION
%-----------------------------------------------------------------------

%----set rendering window view parameters
f_handle = 15;
axis_limits = [-10 10 -10 10 -10 10];
render_view = [1 1 1]; view_up = [0 0 1];
SetRenderingViewParameters(axis_limits,render_view,...
                           view_up,f_handle); 
                       
%----initialize rendering
% link 0 rendering initialization 
L0 = 4; r0 = 1; sides0 = 4; axis0 = 3; norm_L0 = -1.0; 
linkColor0 = [0 .3 .3]; plotFrame0 = 0; 
d0 = CreateLinkRendering(L0,r0,sides0,axis0,norm_L0,linkColor0,... 
    plotFrame0,f_handle); 

% link 1 rendering initialization
L1 = 4; r1 = .75; sides1 = 4; axis1 = 3; norm_L1 = 1.0;
    linkColor1 = [0 0.75 0]; plotFrame1 = 0;
d1 = CreateLinkRendering(L1,r1,sides1,axis1,norm_L1,linkColor1,...
                          plotFrame1,f_handle);

% link 2 rendering initialization
L2 = 4; r2 = .5; sides2 = 12; axis2 = 3; norm_L2 = -1.0;
    linkColor2 = [0.75 0 0]; plotFrame2 = 0;
d2 = CreateLinkRendering(L2,r2,sides2,axis2,norm_L2,linkColor2,...
                          plotFrame2,f_handle);
                      
% link 3 rendering initialization
L3 = 4; r3 = .3; sides3 = 4; axis3 = 1; norm_L3 = 1.0;
    linkColor3 = [0 0 0.75]; plotFrame3 = 0;
d3 = CreateLinkRendering(L3,r3,sides3,axis3,norm_L3,linkColor3,...
                          plotFrame3,f_handle);

%-------------------------------------------------------------------
%  DISPLAY INTERATION RESULTS
%-------------------------------------------------------------------
% Tfixed = [1  0  0  0 
%        0  1  0  0 
%        0  0  1  0 
%        0  0   0  1]; 
% UpdateLink(d1,Tfixed); 

for i = 1:k-1
    % Update frame {1}
    q1 = (Q(1,i)); 
    T_10 = [1   0   0   0
           0   1   0   0
           0   0   1   q1
           0   0   0   1];
       
    % Update frame {2}
    q2 = (Q(2,i));
    T_20 = [-sin(q2)    0    cos(q2)    0
            cos(q2)    0    sin(q2)    0
                 0     1          0    q1
                 0     0          0    1];
        
    % Update frame {3}
    q3 = (Q(3,i));
      T_e0 = [  cos(q2)  -sin(q2)   0   q3*cos(q2)
                sin(q2)   cos(q2)   0   q3*sin(q2)
                   0         0      1       q1
                   0         0      0        1];
    
    UpdateLink(d1,T_10);
    UpdateLink(d2,T_20);
    UpdateLink(d3,T_e0);
    
%     if i == 1;  %pause at start of simulation rendering
%         pause;
%     end
    pause(.02);
      
end


