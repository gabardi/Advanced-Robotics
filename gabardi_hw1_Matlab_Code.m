%% ME 739 Homework #1
% Kaitlyn Gabardi
%02/16/2019

%% Problem 1

R21 = [1        0         0
       0       1/2   -1*sqrt(3)/2
       0    sqrt(3)/2    1/2];

R31 = [0  0 -1
       0  1  0
       1  0  0];
   
R32 = transpose(R21)*R31

%% Problem 2

%homogenous transformations H10,H21, H20

%H10 solved analytically on page 2

H10 = [0 1  0 0
       0 0 -1 0
      -1 0  0 1
       0 0  0 1];
   
%H21 solved analytically on page 2

H21 = [0 -1  0  1
       0  0 -1  0
       1  0  0 -1
       0  0  0  1];
   
%H20 solved analytically on page 2

H20 = [0 0 -1 0
      -1 0  0 1
       0 1  0 0
       0 0  0 1];

% Verified that H20 found analytically equals H20 found by H10*H21
H20_Proof = H10*H21; 

%% Problem 3

H3_2 = [0  1  0  0
        1  0  0  0
        0  0 -1  2
        0  0  0  1];
    
H2_1 = [1  0  0 -1/2
        0  1  0  1/2
        0  0  1   0
        0  0  0   1];

H1_0 = [1 0 0 0
        0 1 0 1
        0 0 1 1
        0 0 0 1];

H3_0 = H1_0*H2_1*H3_2
H2_0 = H1_0*H2_1

%Part 2: NEed to find the inverse of H32

%Rotation Matrix RT
RT_H32 = [0  1  0 
          1  0  0
          0  0 -1];

neg_RT_H32 = -1*RT_H32;

%Displacment dBA
d32 = [0
       0
       2];
    
d32_neg_RT_H32 = neg_RT_H32*d32;

%Inverse of H32 = H23 *can also do inv(H_32)
H23 = [RT_H32  d32_neg_RT_H32 
        0    0    0     1]


%% Problem 4
%Part 1

%Evaluate the velocity of the point expressed in frame {B}, vB
%VB = HAB*VA

%Rotation Matrix RT
RT_HAB = [-1/2        0    sqrt(3)/2
            0         1        0
       -1*sqrt(3)/2   0      -1/2];

%Define vA
VA = [-2; 4; 2];

%Evaluation of velecity expressed in frame {B}
VB = RT_HAB*VA
    
%Part 2 Evaluating the Magnitube of VB and VA
magnitudeVB = norm(VB)
magnitudeVA = norm(VA)

%% Problem 5
%Part A

% Homogeneous transforms that relate successive frames (DH paramter model)
syms q1 q2 q3 L1 L2 L3 d2  %symbolic variables

%defined variables
c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);

% DH Convention Transformation

T10_A = [c1    0     s1     0
         s1    0    -c1     0
         0     1      0     L1
         0     0      0     1];

T21_A = [c2    -s2    0    L2*c2
         s2     c2    0    L2*s2
         0       0    1      0
         0       0    0      1];

T32_A = [c3    -s3    0    L3*c3
         s3     c3    0    L3*c3
         0      0     1      0
         0      0     0      1];


T_A1 = T10_A*T21_A
T_A2 = T10_A*T21_A*T32_A
%pretty(T_A)

%Part B
% DH Convention Transformation

T10_B = [c1    -s1*cosd(0)    s1*sind(0)   (L1+L2+d2)*(c1)
         s1     c1*cosd(0)   -c1*sind(0)   (L1+L2+d2)*(s1)
         0        sind(0)       cosd(0)          0
         0           0               0           1];

T21_B = [c3    -s3*cosd(0)    s3*sind(0)   (L1+L2+d2)*(c3)
         s3     c3*cosd(0)   -c3*sind(0)   (L1+L2+d2)*(s3)
         0        sind(0)       cosd(0)          0
         0           0             0             1];


T_B = T10_B*T21_B;
pretty(T_B)



%% Problem 6

% Homogeneous transforms as a funciton of joint variables
syms q1 q2 q3 q4 q5 q6  %symbolic variables

%defined variables
c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);
c4 = cos(q4);
c5 = cos(q5);
c6 = cos(q6);

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);
s4 = sin(q4);
s5 = sin(q5);
s6 = sin(q6);

%Homogenous Transforms
H10_6 = [0  s1   c1  0
         1   0    0  3
         0  c1  -s1  0
         0   0    0  1];
   
H21_6 = [c2   0   s2   0
         0    1    0   0
        -s2   0   c2   0
         0    0    0   1];

H32_6 = [c3   s3   0    0
          0    0  -1    0
        -s3   c3   0    5
          0    0   0    1];
   
H43_6 = [1   0   0    0
         0   1   0  (3+q4)
         0   0   1    0
         0   0   0    1];

H54_6 = [c5   s5   0    0
        -s5   c5   0    0
          0    0   1    0
          0    0   0    1];

H65_6 = [-c6   s6   0    0
          0    0    1    3
          s6   c6   0    0
          0    0    0    1];

% H76_6 = [1  0  0 -1
%          0  1  0  0
%          0  0  1  0
%          0  0  0  1];

%Homogenous transformation matrices- frame displacement relative to {0}

H10_6 
H20_6 = H10_6*H21_6
H30_6 = H20_6*H32_6
H40_6 = H30_6*H43_6
H50_6 = H40_6*H54_6
H60_6 = H50_6*H65_6
%H76_6 = H60_6*H76_6

%% Plots of x,y,z position of the end-effector (point E)
%given q1-q6 parameters
% q1 = -pi*sin(t);
% q2 = (pi/4)*(1-cos(t));
% q3 = (pi/4)*sin(t);
% q4 = 1.5*(1-cos(t));
% q5 = (-pi/4)*sin(t);
% q6 = (pi/4)*sin(t);

t = 0:.01:2*pi; %animation of position over given time period
for i = 1:length(t)
    q1 = -pi*sin(t(i));
    q2 = (pi/4)*(1-cos(t(i)));
    q3 = (pi/4)*sin(t(i));
    q4 = 1.5*(1-cos(t(i)));
    q5 = (-pi/4)*sin(t(i));
    q6 = (pi/4)*sin(t(i));
    
    %Homogenous Transforms
    H10_6 = [0   sin(q1)   cos(q1)    0
             1      0         0       3
             0   cos(q1)  -1*sin(q1)  0
             0      0         0       1];
   
    H21_6 = [cos(q2)   0    sin(q2)   0
               0       1       0      0
           -1*sin(q2)  0    cos(q2)   0
               0       0       0      1];

    H32_6 = [cos(q3)   sin(q3)   0    0
               0         0      -1    0
          -1*sin(q3)   cos(q3)   0    5
               0         0       0    1];
   
    H43_6 = [1   0   0    0
             0   1   0  (3+q4)
             0   0   1    0
             0   0   0    1];

    H54_6 = [cos(q5)    sin(q5)   0    0
           -1*sin(q5)   cos(q5)   0    0
               0           0      1    0
               0           0      0    1];

    H65_6 = [-1*cos(q6)  sin(q6)   0    0
                 0         0      -1    3
              sin(q6)    cos(q6)   0    0
                 0         0       0    1];
             
    %H76_6 =  [1  0  0 -1
    %          0  1  0  0
    %          0  0  1  0
    %          0  0  0  1]; Shift in end effector


    H60_6_End_Effector_nonsyms = H10_6*H21_6*H32_6*H43_6*H54_6*H65_6
    
    end_effector_position(:,i) = H60_6_End_Effector_nonsyms(:,4);
    
end

%plot of x-position of end-effector as a function of time
figure(1)
plot(t,end_effector_position(1,:))
xlabel('time')
ylabel('x')
grid on
   
%plot of y-position of end-effector as a function of time
figure(2)
plot(t,end_effector_position(2,:))
xlabel('time')
ylabel('y')
grid on
    
%plot of z-position of end-effector as a function of time
figure(3)
plot(t,end_effector_position(3,:))
xlabel('time')
ylabel('z')
grid on

%% Problem 6 Animation of Manipulator 
figure(4)

%---------------set rendering window view parameters-----------------------
%figure handle
f_handle = 1; 

%axis limits
axis_limits = [-10 10 0 20 -10 10];

%camera position
render_view = [-1 1 -1]; 

%verticle orientation
view_up = [0 1 0];

%initialize rendering view
SetRenderingViewParameters(axis_limits,render_view,view_up,f_handle);

%-----------------------initialize rendering-------------------------------
% base (link 0) rendering initialization
r0 = 2.5; L0 = 1.0; sides0 = 20.0; axis0 = 2.0; norm_L0 = -1.0;
linkColor0 = [.5 .5 .5]; plotFrame0 = 0 ;
d0 = CreateLinkRendering(L0,r0,sides0,axis0,norm_L0,linkColor0,...
    plotFrame0,f_handle);

% link 1 rendering initialization
r1 = 2.0; L1 = 3.0; sides1 = 20.0; axis1 = 1.0; norm_L1 = +1.0;
linkColor1 = [0 1 0]; plotFrame1 = 0;
d1 = CreateLinkRendering(L1,r1,sides1,axis1,norm_L1,linkColor1,...
plotFrame1,f_handle);

% link 2 rendering initialization
r2 = 1.25; L2 = 5.0; sides2 = 4.0; axis2 = 3.0; norm_L2 = -1.0; %change 1 2
linkColor2 = [.5 .5 .5]; plotFrame2 = 0;
d2 = CreateLinkRendering(L2,r2,sides2,axis2,norm_L2,linkColor2,...
plotFrame2,f_handle);

% link 3 rendering initialization
r3 = 1.0; L3 = 3.0; sides3 = 4.0; axis3 = 2.0; norm_L3 = -1.0;
linkColor3 = [0 1 0]; plotFrame3 = 0;
d3 = CreateLinkRendering(L3,r3,sides3,axis3,norm_L3,linkColor3,...
plotFrame3,f_handle);

% link 4 rendering initialization
r4 = .8; L4 = 3.0; sides4 = 4.0; axis4 = 2.0; norm_L4 = +1.0;
linkColor4 = [.5 .5 .5]; plotFrame4 = 0;
d4 = CreateLinkRendering(L4,r4,sides4,axis4,norm_L4,linkColor4,...
plotFrame4,f_handle);

% link 5 rendering initialization
r5 = .75; L5 = 2.0; sides5 = 4.0; axis5 = 2.0; norm_L5 = -1.0;
linkColor5 = [0 1 0]; plotFrame5 = 0;
d5 = CreateLinkRendering(L5,r5,sides5,axis5,norm_L5,linkColor5,...
plotFrame5,f_handle);

% end-effector frame rendering initialization (link 6)
r6 = .5; L6 = 1.0; sides6 = 4.0; axis6 = 3.0; norm_L6 = -1.0;
linkColor6 = [.5 .5 .5]; plotFrame6 = 0;
d6 = CreateLinkRendering(L6,r6,sides6,axis6,norm_L6,linkColor6,...
plotFrame6,f_handle);

%----calculate robot manipulator link displacements as a function of time
tEnd = 2*pi; SamplesPerSec = 20;
t = linspace(0,tEnd,tEnd*SamplesPerSec)';

for i = 1:size(t,1)
    q1 = -pi*sin(t(i));
    q2 = (pi/4)*(1-cos(t(i)));
    q3 = (pi/4)*sin(t(i));
    q4 = 1.5*(1-cos(t(i)));
    q5 = (-pi/4)*sin(t(i));
    q6 = (pi/4)*sin(t(i));
    
%transformation from {1} to {0}
    T10(:,:,i) = [0   sin(q1)   cos(q1)    0
                  1      0         0       3
                  0   cos(q1)  -1*sin(q1)  0
                  0      0         0       1];
              
%transformation from {2} to {1}

    T21(:,:,i) =  [cos(q2)   0    sin(q2)  0
                    0       1       0      0
                -1*sin(q2)  0    cos(q2)   0
                    0       0       0      1];
                
%transformation from {3} to {2}
    T32(:,:,i) =  [cos(q3)    sin(q3)   0    0
                     0          0      -1    0
                 -1*sin(q3)   cos(q3)   0    5
                     0          0       0    1];

%transformation from {4} to {3}
    T43(:,:,i) = [1   0   0    0
                  0   1   0  (3+q4)
                  0   0   1    0
                  0   0   0    1];
              
%transformation from {5} to {4}
    T54(:,:,i) = [cos(q5)    sin(q5)   0    0
                -1*sin(q5)   cos(q5)   0    0
                    0           0      1    0
                    0           0      0    1];


%transformation from the end-effector{6} to {5}
    T65(:,:,i) = [-1*cos(q6)  sin(q6)   0    0
                      0         0      -1    3
                   sin(q6)    cos(q6)   0    0
                      0         0       0    1];
              
    T20(:,:,i) = T10(:,:,i)*T21(:,:,i);
    T30(:,:,i) = T20(:,:,i)*T32(:,:,i);
    T40(:,:,i) = T30(:,:,i)*T43(:,:,i);
    T50(:,:,i) = T40(:,:,i)*T54(:,:,i);
    T60(:,:,i) = T50(:,:,i)*T65(:,:,i);

end

%------loop through the calculated robot motions and update rendering------
for i = 1:size(t,1)
% update the link rendering
UpdateLink(d1,T10(:,:,i));
UpdateLink(d2,T20(:,:,i));
UpdateLink(d3,T30(:,:,i));
UpdateLink(d4,T40(:,:,i));
UpdateLink(d5,T50(:,:,i));
UpdateLink(d6,T60(:,:,i));
% if i == 1
% pause; % to allow resizing of graphics window
% end
pause(0.05)
end
