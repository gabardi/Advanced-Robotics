%% ME/ECE 739 HW #5 
% Kaitlyn Gabardi

%% Problem 1: Joint-space Control
%% Part A: Inverse Dynamics Controller with Full Decoupling of Joints
%% Part B. Inverse Dynamics Controller with Gravity Only Compensation
%% Part C. Decentralized Joint Space Controller with Full Comensation
% Using Average Mass Matrix

%-----------------------------------------------------------------------
%  NUMERICAL INTEGRATION OF DYNAMIC EQUATIONS
%-----------------------------------------------------------------------
% close all; clc; 

% model parameters
modelParameters = InitializeThreeDOF_RRR();

% Create options menu for the three different types of controllers
% User will be able to input desired controller choice
% Will be able to compare and constrast different decoupling compensation

disp('Chose Controller Option:') 
disp('   1 = Controller Design (a): Joint-space control with full decoupling') 
disp('   2 = Controller Design (b): Modified joint-space control with gravity compensation only') 
disp('   3 = Controller Design (c): Decentralized joint-space control') 
modelParameters.controlMethod = input(' '); 
 
%depending on the control type, set different values. The decentralized 
%was found to require a smaller step size for RK-4 to behave well. 
if modelParameters.controlMethod ==1 
    tend = 1;        %simulation run time 
    dT = .005;       %integration step size  
    elseif modelParameters.controlMethod ==2 
    tend = 1;        %simulation run time 
    dT = .003;       %integration step size  
    else 
    tend = 1;        %simulation run time 
    dT = .002;       %integration step size  
end 
 
 
numPts = floor(tend/dT);
q = zeros(3,numPts);         %pre-allocate array memory (to improve sim. speed)
qd = zeros(3,numPts);
t = zeros(1,numPts);

% specify initial variables
q(:,1)  = [0; pi/4; -pi/2];    %initial position
qd(:,1) = [0; 0; 0];           %initial velocity
z = [q(:,1); qd(:,1)];         %initialize the state variables

% specify final variables

qDes(:,1)  = [pi; pi; pi/2];   %final joint space position
qdDes(:,1) = [0; 0; 0];        %final joint velocity
zDes = [qDes; qdDes];          %initialize the state variables
Xend = zeros(3,numPts);        %defined the end effector coordinates
Xvend = zeros(3,numPts);
    
% integrate equations of motion
for i = 1:numPts-1
    % Runge-Kutta 4th order
    k1 = zDot3dof(z,zDes,modelParameters);
    k2 = zDot3dof(z + 0.5*k1*dT,zDes,modelParameters);
    k3 = zDot3dof(z + 0.5*k2*dT,zDes,modelParameters);
    k4 = zDot3dof(z + k3*dT,zDes,modelParameters);
    z = z + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dT;
    
    % store joint position and velocity for post processing
    q(:,i+1) = z(1:3);
    qd(:,i+1) = z(4:6);
    t(1,i+1) = t(1,i) + dT;        
end

%-----------------------------------------------------------------------
%  RENDERING INITIALIZATION
%-----------------------------------------------------------------------
%----set rendering window view parameters
L = modelParameters.L;
L1 = L; 
L2 = L;  
L3 = L;
f_handle = 1;
axis_limits = L*[-2 2 -2 3 -1 3];
render_view = [1 -1 1]; view_up = [0 0 1];
SetRenderingViewParameters(axis_limits,render_view,...
                           view_up,f_handle); 
camproj perspective % turns on perspective

%----initialize rendering

% link 1 rendering initialization
r1 = L1/5; sides1 = 10; axis1 = 2; norm_L1 = 1.0;
    linkColor1 = [0 0.75 0]; plotFrame1 = 0;
d1 =  CreateLinkRendering(L1,r1,sides1,axis1,norm_L1,linkColor1,...
                          plotFrame1,f_handle);

% link 2 rendering initialization
r2 = L2/6; sides2 = 4; axis2 = 1; norm_L2 = 1.0;
    linkColor2 = [0.75 0 0]; plotFrame2 = 0;
d2 =  CreateLinkRendering(L2,r2,sides2,axis2,norm_L2,linkColor2,...
                          plotFrame2,f_handle);

% link 3 rendering initialization
r3 = L2/8; sides3 = 4; axis3 = 1; norm_L3 = 1.0;
    linkColor3 = [0 0 0.75]; plotFrame3 = 0;
d3 =  CreateLinkRendering(L3,r3,sides3,axis3,norm_L3,linkColor3,...
                          plotFrame3,f_handle);
                      
%-------------------------------------------------------------------
%  DISPLAY INTERATION RESULTS
%-------------------------------------------------------------------
for i = 1:numPts
    
    % Update frame {1}
    c = cos(q(1,i)); 
    s = sin(q(1,i)); 
    L = L1;
    T10 = [c  0   s  0
           s  0  -c  0
           0  1   0  L
           0  0   0  1];
       
    % Update frame {2}
    c = cos(q(2,i)); 
    s = sin(q(2,i)); 
    L = L2;
    T21 = [c  -s  0  L*c
           s   c  0  L*s
           0   0  1   0
           0   0  0   1];
       
    % Update frame {3}
    c = cos(q(3,i)); 
    s = sin(q(3,i)); 
    L = L3;
    T32 = [c  -s  0  L*c
           s   c  0  L*s
           0   0  1   0
           0   0  0   1];
    T20 = T10*T21;
    T30 = T20*T32;
    
    % end-effector position in task space 
    Xend(:,i) = T30(1:3,4);
    
    % Plot linear end effector velocity over time by using the basic Jacobian
    % Symbolic code was given and used to compute Jacobian
    
    % Error in original Jacobian, needed to define q1,q2,q3 since Jacobian
    % is currerently symbolic
    q1 = q(1,i); 
    q2 = q(2,i); 
    q3 = q(3,i); 
   
    Jacobian = [L*sin(q1)*sin(q2)*sin(q3) - L*cos(q2)*cos(q3)*sin(q1) - L*cos(q2)*sin(q1),                                                                                                    -cos(q1)*(L*(sin(q2) + 1) - L + L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2)),                                                                                  -cos(q1)*(L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2))
                L*cos(q1)*cos(q2) + L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3),                                                                                                    -sin(q1)*(L*(sin(q2) + 1) - L + L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2)),                                                                                  -sin(q1)*(L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2))
                                                                                        0, cos(q1)*(L*cos(q1)*cos(q2) + L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3)) + sin(q1)*(L*cos(q2)*sin(q1) + L*cos(q2)*cos(q3)*sin(q1) - L*sin(q1)*sin(q2)*sin(q3)), cos(q1)*(L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3)) + sin(q1)*(L*cos(q2)*cos(q3)*sin(q1) - L*sin(q1)*sin(q2)*sin(q3))];
 
    Xvend(:,i) = Jacobian*qd(:,i); 

    % update rendering 
    figure(f_handle); 
    UpdateLink(d1,T10); 
    UpdateLink(d2,T20); 
    UpdateLink(d3,T30); 
    hold on;
    plot3(Xend(1,1:i),Xend(2,1:i),Xend(3,1:i),'r','LineWidth',2); 
 
    %if i == 1;  %pause at start of simulation rendering 
    %    pause; 
    %end 
   
end 

%%%%%%%%%%% Plot joint-space displacements as a function of time %%%%%%%%%%
figure(2); 
plot(t(1:i),q(1,1:i),'b','LineWidth',2); hold on 
plot(t(1:i),q(2,1:i),'r','LineWidth',2); hold on 
plot(t(1:i),q(3,1:i),'k','LineWidth',2); hold off 
title('Joint-Space Displacement as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Position (rad)'); 
grid on 
legend('Joint 1','Joint 2','Joint 3','Location','best'); 
 
%%%%%%%%%%%% Plot joint-space velocities as a function of time %%%%%%%%%%%%
figure(3); 
plot(t(1:i),qd(1,1:i),'b','LineWidth',2); hold on 
plot(t(1:i),qd(2,1:i),'r','LineWidth',2); hold on 
plot(t(1:i),qd(3,1:i),'k','LineWidth',2); hold off 
title('Joint-Space Velocity as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Velocity (r/s)'); 
grid on 
legend('Joint 1','Joint 2','Joint 3','Location','best'); 

 
%%%% Plot operational-point end effector displacement as a fxn of time %%%%
figure(4); 
plot(t,Xend(1,:),'b','LineWidth',2); hold on 
plot(t,Xend(2,:),'r','LineWidth',2); hold on 
plot(t,Xend(3,:),'k','LineWidth',2); hold off 
title('End Effector Displacement as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Xe Ye Ze Position (m)'); 
grid on 
legend('Xe','Ye','Ze','Location','best'); 

%%%%%% Plot operational-point end effector velocity as a fxn of time %%%%%%
figure(5); 
plot(t,Xvend(1,:),'b','LineWidth',2); hold on 
plot(t,Xvend(2,:),'r','LineWidth',2); hold on 
plot(t,Xvend(3,:),'k','LineWidth',2); hold off 
title('End Effector Velocity as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Velocity (m/s)'); 
grid on 
legend('dX_e/dt','dY_e/dt','dZ_e/dt','Location','best'); 

%%%%%%%%%%%% Plot the operational-point dispalcements in 3D %%%%%%%%%%%%%%%
figure(6); 
plot3(Xend(1,:),Xend(2,:),Xend(3,:),'LineWidth',2); hold on; 
plot3(Xend(1,1),Xend(2,1),Xend(3,1),'ro'); 
plot3(Xend(1,end),Xend(2,end),Xend(3,end),'ro'); 
title('Operational-Point 3D Displacements'); 
xlabel('X'); 
ylabel('Y'); 
zlabel('Z'); 
grid on

%% Problem 2: Operational-space Control
%% Part E: Operational-space Controller with Full Decoupling of Joints
%% Part F. Simulation of Response to Step Input 
%% Part G. Operational-space Linear Velocity 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------
%  NUMERICAL INTEGRATION OF DYNAMIC EQUATIONS
%-----------------------------------------------------------------------
close all; clc; 
% model parameters
modelParameters = InitializeThreeDOF_OP_model();
  
% Implementation of operational-space inverse dynamics
% Simulation of the resposne of the system to a unit step command
% Part g will use the velocity limitingh heuristic described in the lecture
% notes

disp('Choose Final Operational Position:') 
disp('   1 = Problem 2(e): [-L -L 0]') 
disp('   2 = Problem 2(f): [-L -L/10 0]') 
disp('   3 = Problem 2(g): [-L -L 0] w/ Velocity Heuristic') 
modelParameters.controlMethod = input(' '); 

tend = 1; %vary run time throughout simulations 
dT = 0.005; 
numPts = floor(tend/dT);
q = zeros(3,numPts);         %pre-allocate array memory (to improve sim. speed)
qd = zeros(3,numPts);
t = zeros(1,numPts);

% specify initial variables
q(:,1)  = [0; pi/4; -pi/2];    %initial position
qd(:,1) = [0; 0; 0];           %initial velocity
z = [q(:,1); qd(:,1)];         %initialize the state variables
L = modelParameters.L; 

% specify final variables
xdDes = [0; 0; 0];             %final joint velocity

Xend = zeros(3,numPts);        %defined the end effector coordinates
Xvend = zeros(3,numPts);       %desired final task space velocity 

if modelParameters.controlMethod ==1 
    xDesired = [0; -L; L; 0; 0; 0;];    %final operational-space position/velocity
    elseif modelParameters.controlMethod ==2 
    xDesired = [-L; -L/10; 0; 0; 0; 0];  %final operational-space position/velocity
    else 
    xDesired = [-L; -L; 0; 0; 0; 0];     %final operational-space position/velocity
    modelParameters.vMax=5;
end 

%-----------------------------------------------------------------------
%  INTEGRATE EQUATIONS OF MOTION 
%-----------------------------------------------------------------------

for i = 1:numPts-1
    modelParameters.k = i;
    % Runge-Kutta 4th order
    k1 = zDot3dof_3_OP(z,xDesired,modelParameters);
    k2 = zDot3dof_3_OP(z + 0.5*k1*dT,xDesired,modelParameters);
    k3 = zDot3dof_3_OP(z + 0.5*k2*dT,xDesired,modelParameters);
    k4 = zDot3dof_3_OP(z + k3*dT,xDesired,modelParameters);
    z = z + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dT;
    
    % store joint position and velocity for post processing
    q(:,i+1) = z(1:3);
    qd(:,i+1) = z(4:6);
    t(1,i+1) = t(1,i) + dT;
    
end

%-----------------------------------------------------------------------
%  RENDERING INITIALIZATION
%-----------------------------------------------------------------------

%----set rendering window view parameters
L = modelParameters.L;
L1 = L;  
L2 = L;  
L3 = L;
f_handle = 1;
axis_limits = L*[-3 2 -3 3 -1 2];
render_view = [1 -1 1]; view_up = [0 0 1];
SetRenderingViewParameters(axis_limits,render_view,...
                           view_up,f_handle); 
camproj perspective % turns on perspective

%----initialize rendering

% link 1 rendering initialization
r1 = L1/5; sides1 = 10; axis1 = 2; norm_L1 = 1.0;
    linkColor1 = [0 0.75 0]; plotFrame1 = 0;
d1 = CreateLinkRendering(L1,r1,sides1,axis1,norm_L1,linkColor1,...
                          plotFrame1,f_handle);

% link 2 rendering initialization
r2 = L2/6; sides2 = 4; axis2 = 1; norm_L2 = 1.0;
    linkColor2 = [0.75 0 0]; plotFrame2 = 0;
d2 = CreateLinkRendering(L2,r2,sides2,axis2,norm_L2,linkColor2,...
                          plotFrame2,f_handle);

% link 3 rendering initialization
r3 = L3/8; sides3 = 4; axis3 = 1; norm_L3 = 1.0;
    linkColor3 = [0 0 0.75]; plotFrame3 = 0;
d3 = CreateLinkRendering(L3,r3,sides3,axis3,norm_L3,linkColor3,...
                          plotFrame3,f_handle);
                      
%-------------------------------------------------------------------
%  DISPLAY INTERATION RESULTS
%-------------------------------------------------------------------
for i = 1:numPts
    % Update frame {1}
    c = cos(q(1,i)); 
    s = sin(q(1,i)); 
    L = L1;
    T10 = [c  0   s  0
           s  0  -c  0
           0  1   0  L
           0  0   0  1];
       
    % Update frame {2}
    c = cos(q(2,i)); 
    s = sin(q(2,i)); 
    L = L2;
    T21 = [c  -s  0  L*c
           s   c  0  L*s
           0   0  1   0
           0   0  0   1];
       
    % Update frame {3}
    c = cos(q(3,i)); 
    s = sin(q(3,i)); 
    L = L3;
    T32 = [c  -s  0  L*c
           s   c  0  L*s
           0   0  1   0
           0   0  0   1];
       
    T20 = T10*T21;
    T30 = T20*T32;
    
    % end-effector position
    Xend(:,i) = T30(1:3,4);
    
     % update rendering        
     figure(f_handle);        
     hold on; 
     plot3(Xend(1,1:i),Xend(2,1:i),Xend(3,1:i),'r','LineWidth',2);        
     UpdateLink(d1,T10);        
     UpdateLink(d2,T20);        
     UpdateLink(d3,T30);        
     
%      if i == 1  %pause at start of simulation rendering            
%          pause;
%      end
     pause(dT);
     
     q1 = q(1,i);        
     q2 = q(2,i);        
     q3 = q(3,i);
 
    % Used the MATLAB symbolic toolbox as suggested in HW
    % Symbolic evaluation of the Jacobian and the time derivative
    % See Three_DOF_symbolic
     Jacobian = [L*sin(q1)*sin(q2)*sin(q3) - L*cos(q2)*cos(q3)*sin(q1) - L*cos(q2)*sin(q1), -cos(q1)*(L*(sin(q2) + 1) - L + L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2)), -cos(q1)*(L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2))
                 L*cos(q1)*cos(q2) + L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3), -sin(q1)*(L*(sin(q2) + 1) - L + L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2)),-sin(q1)*(L*cos(q2)*sin(q3) + L*cos(q3)*sin(q2))
                                                                                         0, cos(q1)*(L*cos(q1)*cos(q2) + L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3)) + sin(q1)*(L*cos(q2)*sin(q1) + L*cos(q2)*cos(q3)*sin(q1) - L*sin(q1)*sin(q2)*sin(q3)), cos(q1)*(L*cos(q1)*cos(q2)*cos(q3) - L*cos(q1)*sin(q2)*sin(q3)) + sin(q1)*(L*cos(q2)*cos(q3)*sin(q1) - L*sin(q1)*sin(q2)*sin(q3))];
     Xvend(:,i) = Jacobian*qd(:,i);
end
modelParameters.controlMethod =1;
 
%%%%%%%%%%% Plot joint-space displacements as a function of time %%%%%%%%%%
figure(2); 
plot(t(1:i),q(1,1:i),'b','LineWidth',2); hold on 
plot(t(1:i),q(2,1:i),'r','LineWidth',2); hold on 
plot(t(1:i),q(3,1:i),'k','LineWidth',2); hold off 
title('Joint Displacement as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Position (rad)'); 
grid on 
legend('Joint 1','Joint 2','Joint 3','Location','best'); 
 
%%%%%%%%%%%% Plot joint-space velocities as a function of time %%%%%%%%%%%%
figure(3); 
plot(t(1:i),qd(1,1:i),'b','LineWidth',2); hold on 
plot(t(1:i),qd(2,1:i),'r','LineWidth',2); hold on 
plot(t(1:i),qd(3,1:i),'k','LineWidth',2); hold off 
title('Joint Velocity as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Velocity (r/s)'); 
grid on 
legend('Joint 1','Joint 2','Joint 3','Location','best'); 

 
%%%% Plot operational-point end effector displacement as a fxn of time %%%%
figure(4); 
plot(t,Xend(1,:),'b','LineWidth',2); hold on 
plot(t,Xend(2,:),'r','LineWidth',2); hold on 
plot(t,Xend(3,:),'k','LineWidth',2); hold off 
title('End Effector Displacement as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Xe Ye Ze Position (m)'); 
grid on 
legend('Xe','Ye','Ze','Location','best'); 

%%%%%% Plot operational-point end effector velocity as a fxn of time %%%%%%
figure(5); 
plot(t,Xvend(1,:),'b','LineWidth',2); hold on 
plot(t,Xvend(2,:),'r','LineWidth',2); hold on 
plot(t,Xvend(3,:),'k','LineWidth',2); hold off 
title('End Effector Velocity as a Function of Time'); 
xlabel('Time (s)'); 
ylabel('Velocity (m/s)'); 
grid on 
legend('dX_e/dt','dY_e/dt','dZ_e/dt','Location','best'); 

%%%%%%%%%%%% Plot the operational-point dispalcements in 3D %%%%%%%%%%%%%%%
figure(6); 
plot3(Xend(1,:),Xend(2,:),Xend(3,:),'LineWidth',2); hold on; 
plot3(Xend(1,1),Xend(2,1),Xend(3,1),'ro'); 
plot3(Xend(1,end),Xend(2,end),Xend(3,end),'ro'); 
title('Operational-Point 3D Displacements'); 
xlabel('X'); 
ylabel('Y'); 
zlabel('Z'); 
grid on