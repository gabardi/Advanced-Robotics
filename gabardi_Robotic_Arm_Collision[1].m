%% ME 739 FINAL PROJECT with Obsticle
% Kaitlyn Gabardi
%5/04/2019

%% Problem 2
%-----------------------------------------------------------------------
%  INITIALIZATION
%-----------------------------------------------------------------------
clear all; close all; clc
% input parameters
L1 = 1;                        %manipulator link lengths
L2 = 1;              
L3 = 1; 

zeta1 = 100;                   %attractive force scaling
zeta2 = 10;                  
zeta3 = 5;

eta1 = 1;                      %repulsive force scaling
eta2 = 1;  
eta3 = 1;

qGoal = [pi/2 pi/2 0]';           %goal position (joint space)

b = [8 7 -5]';                  %obstacle location
q0 = [0 0 0]';                 %initial manipulator joint position
rho0 = 1;                      %radius of influence - repulsive forces

alpha = 0.01;                %gradient descent step size
epsilon = 0.001;             %gradient descent error tolerance
maxIterations = 1000;        %max gradient descent iterations

q = q0;                      %initialize joint positions
k = 1;                       %initialize counter

%trigometric expressions
c1g = cos(qGoal(1));  
s1g = sin(qGoal(1));

c1g2g = cos(qGoal(1)+qGoal(2));
s1g2g = sin(qGoal(1)+qGoal(2));

c1g3g = cos(qGoal(1)+qGoal(2)+qGoal(3));
s1g3g = sin(qGoal(1)+qGoal(2)+qGoal(3));
   
%goal position (task space)
o1g = [L1*c1g
       L1*s1g
         L1];
   
o2g = [L1*c1g + L2*c1g2g
       L1*s1g + L2*s1g2g
              L1];

o3g = [L1*c1g + L2*c1g2g + L3*c1g3g
       L1*s1g + L2*s1g2g + L3*s1g3g
                   L1];

%-----------------------------------------------------------------------
%  GRADIENT DESCENT
%-----------------------------------------------------------------------
while((q - qGoal)'*(q - qGoal) > epsilon)
    %trigometric expresions
    s1 = sin(q(1)); 
    s2 = sin(q(2)); 
    s3 = sin(q(3)); 
    s12 = sin(q(1) + q(2));
    s123 = sin(q(1) + q(2) + q(3));
    
    c1 = cos(q(1)); 
    c2 = cos(q(2)); 
    c3 = cos(q(3)); 
    c12 = cos(q(1) + q(2));
    c123 = cos(q(1) + q(2) + q(3));

    %forward kinematics
    o1 = [L1*c1
          L1*s1
            0];

    o2 = [L2*c1 + L2*c1*s2
          L1*s1 + L2*s1*c2 
              L2 *s2];
     
    o3 = [L1*c1 + L2*c1*c2 - L3*s1*s3 + L3*c1*c2*c3
          L1*s1 + L2*c2*s1 + L3*c1*s3 + L3*c2*c3*s1
                     L2*s2 + L3*c3*s2];
      
    %attractive forces
    FA1 = -zeta1*(o1 - o1g);
    FA2 = -zeta2*(o2 - o2g);
    FA3 = -zeta3*(o3 - o3g);

    %repulsive forces
    %  distance from obstacle
    rho1 = sqrt((o1 - b)'*(o1 - b));
    rho2 = sqrt((o2 - b)'*(o2 - b));
    rho3 = sqrt((o3 - b)'*(o3 - b));
      
    %  gradient of distance
    grad_rho1 = (o1 - b)/rho1;
    grad_rho2 = (o2 - b)/rho2;
    grad_rho3 = (o3 - b)/rho3;
    
    %   repulsive forces
    FR1 = [0; 0; 0]; 
    FR2 = [0; 0; 0]; 
    FR3 = [0; 0; 0];
    
    if rho1 < rho0
        FR1 = eta1*(1/rho1 - 1/rho0)*(1/(rho1^2))*grad_rho1;
    end
    
    if rho2 < rho0
        FR2 = eta2*(1/rho2 - 1/rho0)*(1/(rho2^2))*grad_rho2;
    end
    
    if rho3 < rho0
        FR3 = eta3*(1/rho3 - 1/rho0)*(1/(rho3^2))*grad_rho3;
    end
    
    %sum of forces on each point
    Fo1 = FA1 + FR1;
    Fo2 = FA2 + FR2;
    Fo3 = FA3 + FR3;
    
%point Jacobians
    Jo1 = [-L1*s1   0   0
            L1*c1   0   0
              0     0   0];
   
%     o2 = [L2*c1 + L2*c1*s2
%           L1*s1 + L2*s1*c2 
%               L2*s2];

    Jo2 = [-L1*s1 - L2*s1*s2    L2*c1*c2      0  
            L1*c2 - L2*s1*s2    -L2*s1*s2      0  
                  0                 0         0];
        
%     o3 = [L1*c1 + L2*c1*c2 - L3*s1*s3 + L3*c1*c2*c3
%           L1*s1 + L2*c2*s1 + L3*c1*s3 + L3*c2*c3*s1
%                      L2*s2 + L3*c3*s2];

    Jo3 = [-L1*s1 - L2*s1*c2 - L3*c1*s3 - L3*s1*c2*c3          -L3*c1*s2*c3        -L3*s1*c3 - L3*c1*c2*s3
            L1*c1 + L2*c2*c1 - L3*s1*s3 + L3*c2*c3*c1         -L2*s2*s1 - L3*s2*c3*s1           L3*c1*c3 - L3*c2*s3*s1
                  0             L2*c2           -L3*s3*s2];
    
    %joint torques
    tau = Jo1'*Fo1 + Jo2'*Fo2 + Jo3'*Fo3;
    
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
L = 5;
L1 = L;  
L2 = L;  
L3 = L;
f_handle = 1;
axis_limits = L*[-1 3 -1 4 -1 2.5];
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
                      
% obstacle
r_obs = 1; sides_obs = 15; axis_obs = 3; norm_obs = -1.0;
    linkColor_obs = [0 0 0.75]; plotFrame_obs = 0;
d3_obs = CreateLinkRendering(8,r_obs,sides_obs,axis_obs,norm_obs,linkColor_obs,...
                          plotFrame_obs,f_handle);
To = [1  0  0  b(1)
      0  1  0  b(2)
      0  0  1  b(3)
      0  0  0   1];   
  
%-------------------------------------------------------------------
%  DISPLAY INTERATION RESULTS
%-------------------------------------------------------------------
for i = 1:k-1
    
    % Update frame {1}
    c = cos(Q(1,i)); 
    s = sin(Q(1,i)); 
    L1 = 5;
    
    T10 = [c   0    s   L1*c
           s   0   -c   L1*s
           0   1    0     0
           0   0    0     1];
       
    % Update frame {2}
    c = cos(Q(2,i)); 
    s = sin(Q(2,i)); 
    L2 = 5;
    
    T21 = [c    0   -s   L2*c
           s    0    c   L2*s
           0   -1    0     0
           0    0    0     1];
    
     %transformation from {3} to {2}
    c = cos(Q(3,i)); 
    s = sin(Q(3,i)); 
    L3 = 5;
    
    T32 = [ c    0   -s  L3*c
            s    0    c  L3*s
            0   -1    0   0
            0    0    0   1];
       
    T20 = T10*T21;
    T30 = T20*T32;
    
    UpdateLink(d1,T10);
    UpdateLink(d2,T20);
    UpdateLink(d3,T30);
    UpdateLink(d3_obs,To);
%     if i == 1;  %pause at start of simulation rendering
%         pause;
     F(i) = getframe(gcf);
    pause(.05);
      
end 


% Used to create a video of the simulation to put in powerpoint

video = VideoWriter('3DOF','Uncompressed AVI');
open(video)
writeVideo(video,F);
close(video)

%% Calculations of homoegenous transformation matrices   
    % Update frame {1}
 syms L1 L2 L3 c1 c2 c3 s1 s2 s3   
  
    T10 = [c1   0    s1   L1*c1
           s1   0   -c1   L1*s1
           0    1     0      0
           0    0     0      1];
       
    % Update frame {2}
    
    T21 = [c2    0   -s2   L2*c2
           s2    0    c2   L2*s2
           0    -1    0     0
           0     0    0     1];
    
     %transformation from {3} to {2}
  
    T32 = [ c3    0    -s3  L3*c3
            s3    0     c3  L3*s3
            0    -1      0   0
            0     0      0   1];
       
    T20 = T10*T21;
    T30 = T20*T32;
