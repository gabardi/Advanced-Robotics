%% ME 739 Homework #2
% Kaitlyn Gabardi
%03/11/2019

%% Problem 1
% Rotation transformation matrices

% +45 degree rotation about za
R1 = [cosd(45)  -sind(45)  0  0
      sind(45)   cosd(45)  0  0 
         0         0     1  0 
         0         0     0  1];
  
% 30 degree rotation about xb 
R2 = [1    0         0      0
      0  cosd(30) -sind(30)   0
      0  sind(30)  cosd(30)   0
      0    0        0       1];
  
% -45 degree rotation about za 
R3 = [cosd(-45)  -sind(-45)   0   0
      sind(-45)   cosd(-45)   0   0 
         0          0       1   0
         0          0       0   1]; 
  
% +90 degree rotation about yb
R4 = [cosd(90)  0  sind(90)  0
        0      1    0      0
     -sind(90)  0  cosd(90)  0
        0      0    0      1];

% -30 degree roation about xa
R5 = [1     0         0        0
      0  cosd(-30)  -sind(-30)   0
      0  sind(-30)   cosd(-30)   0
      0     0         0        1];
  

T_1 = round(R5*R3*R1*R2*R4)

%% Problem 3

clear all; close all; clc
% Case 1:
d = 0; a = 10; alpha = 0; theta = 0;
T1 = DH2T(d,a,alpha,theta)

% Case 2:
d = 10; a = 0; alpha = 0; theta = 0;
T2 = DH2T(d,a,alpha,theta)

% Case 3:
d = 10; a = 0; alpha = pi; theta = 0;
T3 = DH2T(d,a,alpha,theta)

% Case 4:
d = 0; a = 0; alpha = pi; theta = pi;
T4 = DH2T(d,a,alpha,theta) 

%% Question 4

%Using DH2T function, evaluations of homogeneous transformation matrices 
clear all; close all; clc

%L = 1;

% Case 1:
d = 1; a = 0; alpha = pi/2; theta = 0;
link1_1 = DH2T(d,a,alpha,theta);

d = 0; a = 0; alpha = pi/2; theta = pi/2;
link2_1 = DH2T(d,a,alpha,theta);

d = 1; a = 0; alpha = 0; theta = pi/2;
link3_1 = DH2T(d,a,alpha,theta);

d = 1; a = 0; alpha = pi/2; theta = 0;
link4_1 = DH2T(d,a,alpha,theta);

H40_Configuration_1 = round(link1_1*link2_1*link3_1*link4_1)

% Case 2:
d = 2; a = 0; alpha = pi/2; theta = 0;
link1_2 = DH2T(d,a,alpha,theta);

d = 0; a = 0; alpha = pi/2; theta = pi;
link2_2 = DH2T(d,a,alpha,theta);

d = 2; a = 0; alpha = 0; theta = pi/2;
link3_2 = DH2T(d,a,alpha,theta);

d = 1; a = 0; alpha = pi/2; theta = pi/2;
link4_2 = DH2T(d,a,alpha,theta);

H40_Configuration_2 = round(link1_2*link2_2*link3_2*link4_2)




%% Question 6
%Part 2 Determining Linear Velocity given qdot

% symbolic variables
syms L2 q1 q2

q=[pi; pi/2; L2];

q_dot=[1; 0; 1];

J_linear_velocity=[-(L2+q(3))*sin(q(1))*cos(q(2))  -(L2+q(3))*cos(q(1))*sin(q(2))   cos(q(1))*cos(q(2))
                    (L2+q(3))*cos(q(1))*cos(q(2))   -(L2+q(3))*sin(q(1))*sin(q(2))  sin(q(1))*cos(q(2))
                                  0                       (L2+q(3))*cos(q(2))             sin(q(2))];
%evaluation of linear velocity
linear_velocity = J_linear_velocity*q_dot

%evaluation of singularities 
det(J_linear_velocity)

%% Problem 7
%Part 1: Forward Kinematics on attached paper

%Part 2: 
L = 1;

syms c1 c2 c3 s1 s2 s3

T10 = [c1  0   s1  L*c1
       s1  0  -c1  L*s1
       0   1    0   0
       0   0    0   1];

T21 = [c2  0  -s2  L*c2
       s2  0   c2  L*s2
       0  -1    0    0
       0   0    0    1];
   
T32 = [c3  0  -s3  L*c3
       s3  0   c2  L*s3 
       0  -1    0    0
       0   0    0    1];

T20 = T10*T21;
T30 = T10*T21*T32;

%z-unit vectors expressed in {0}
z00 = transpose([0 0 1]);
z10 = transpose([s1 -c1 0]);
z20 = transpose([-c1*s1 -s1*s2 c2]);
z30 = transpose([(-c2*s1)-(c1*c2*s3) (c1*c2)-(c2*s1*s3) -s2*s3]); 

%frame origins expressed in {0}
o10 = transpose([L*c1 L*s1 0]);
o20 = transpose([c1+c1*c2 s1+c2*s1 s2]);
o30 = transpose([c1+c1*c2-s1*s3+c1*c2*c3 s1+c2*s1+c1*s3+c2*c3*s1 s2+c3*s2]);

%Calculation of JV1
JV1 = cross(z00,o30);

%Calculation of JV2
origin_JV2 = o30 - o10;
JV2 = cross(z10,origin_JV2);

%Calculation of JV3
origin_JV3 = o30 - o20;
JV3 = cross(z20,origin_JV3);

%Basic jacobian realting joint space velocities and task space velocities
J0 = [JV1 JV2 JV3 
      z00 z10 z20]
    
%% Part 3: Evaluation of required joint torques
%Evaluation of all three configurations 
L = 1;
q1 = 0;
q2 = 0;
q3 = 30;

%defined variables
c1 = cosd(q1);
c2 = cosd(q2);
c3 = cosd(q3);

s1 = sind(q1);
s2 = sind(q2);
s3 = sind(q3);

T10 = [c1  0   s1  L*c1
       s1  0  -c1  L*s1
       0   1    0   0
       0   0    0   1];

T21 = [c2  0  -s2  L*c2
       s2  0   c2  L*s2
       0  -1    0    0
       0   0    0    1];
   
T32 = [c3  0  -s3  L*c3
       s3  0   c2  L*s3 
       0  -1    0    0
       0   0    0    1];

T20 = T10*T21;
T30 = T10*T21*T32;

%z-unit vectors expressed in {0}
z00 = transpose([0 0 1]);
z10 = transpose([s1 -c1 0]);
z20 = transpose([-c1*s1 -s1*s2 c2]);
z30 = transpose([(-c2*s1)-(c1*c2*s3) (c1*c2)-(c2*s1*s3) -s2*s3]); 

%frame origins expressed in {0}
o10 = transpose([L*c1 L*s1 0]);
o20 = transpose([c1+c1*c2 s1+c2*s1 s2]);
o30 = transpose([c1+c1*c2-s1*s3+c1*c2*c3 s1+c2*s1+c1*s3+c2*c3*s1 s2+c3*s2]);

%Calculation of JV1
JV1 = cross(z00,o30);

%Calculation of JV2
origin_JV2 = o30 - o10;
JV2 = cross(z10,origin_JV2);

%Calculation of JV3
origin_JV3 = o30 - o20;
JV3 = cross(z20,origin_JV3);

%Basic jacobian realting joint space velocities and task space velocities
J0_T = transpose([JV1 JV2 JV3 
                  z00 z10 z20])
                 
%Using transpose of jacobian to find joint torques
syms fx Tz %symbolic variables
F = transpose([fx 0 0 0 0 Tz]);
Joint_Torque = J0_T*F

%% Part 4: Evaluation of linear velocity Jacobian, Jv1
%Evaluation of all three configurations 
L = 1;
syms q1 q2 q3

%defined variables
c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);

T10 = [c1  0   s1  L*c1
       s1  0  -c1  L*s1
       0   1    0   0
       0   0    0   1];

T21 = [c2  0  -s2  L*c2
       s2  0   c2  L*s2
       0  -1    0    0
       0   0    0    1];
   
T32 = [c3  0  -s3  L*c3
       s3  0   c2  L*s3 
       0  -1    0    0
       0   0    0    1];

T20 = T10*T21;
T30 = T10*T21*T32;

%z-unit vectors expressed in {0}
z00 = transpose([0 0 1]);
z10 = transpose([s1 -c1 0]);
z20 = transpose([-c1*s1 -s1*s2 c2]);
z30 = transpose([(-c2*s1)-(c1*c2*s3) (c1*c2)-(c2*s1*s3) -s2*s3]); 

%frame origins expressed in {0}
o10 = transpose([L*c1 L*s1 0]);
o20 = transpose([c1+c1*c2 s1+c2*s1 s2]);
o30 = transpose([c1+c1*c2-s1*s3+c1*c2*c3 s1+c2*s1+c1*s3+c2*c3*s1 s2+c3*s2]);

%Calculation of JV1
JV1 = cross(z00,o30);

%Calculation of JV2
origin_JV2 = o30 - o10;
JV2 = cross(z10,origin_JV2);

%Calculation of JV3
origin_JV3 = o30 - o20;
JV3 = cross(z20,origin_JV3);

%Basic jacobian realting joint space velocities and task space velocities

% R01 = [c1   s1   0
%         0    0   1
%        s1  -c1   0];
   
%JV1_linearv_frame1 = simplify(R01*JV0)

JV1 = [-s3           (-c3*s2 - s2)     -c2*s3
        0             (c3*c2 + c2)     -s2*s3
  (-c3*c2 - c2 -1)          0            -c3];

%determining singularity 
determinant = det(JV1);
det_simplified = simplify(determinant)

%Part 7: Finding the analytical Jacobian 

L = 1;
q1_a = 0;
q2_a = pi/2;
q3_a = pi/2;


%defined variables
c1 = cos(q1_a);
c2 = cos(q2_a);
c3 = cos(q3_a);

s1 = sin(q1_a);
s2 = sin(q2_a);
s3 = sin(q3_a);

E = [2   3    0
     1   1   -1
     0   0    1];


J_0 = [-s1 - (c2*s1) - (c1*s3) - (c2*c3*s1),-c1*(s2 + c3*s2),(-c3*s1*(s2^2)) - c2*((c1*s3) + (c2*c3*s1))
       c1 + (c1*c2) - (s1*s3) + (c1*c2*c3), -s1*(s2 + (c3*s2)),(c1*c3*s1*s2) - c2*((s1*s3) - (c1*c2*c3))
        0, c1*((c1*c2) - (s1*s3) + (c1*c2*c3)) + s1*((c2*s1) + (c1*s3) + (c2*c3*s1)), -c1*s1*((c1*s3) + (c2*c3*s1)) - s1*s2*((s1*s3) - (c1*c2*c3))];

Ja = E*J_0

%% Problem 8
%Evaluating the inverse kinematics of 6 DOF manipulator

%Kinematic Decoupling: Solving for q1,q2,q3

syms q1 q2 q3

%defined variables
c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);

Rz_q1 = [1   0   0
          0   1   0
          0   0   1];

Rz_q2 = [c2   -s2   0
         s2    c2   0
          0     0   1];
    
Rx_q3 = [1   0   0
         0   1   0
         0   0   1];
     
R30 = Rz_q1*Rz_q2*Rx_q3

%Inverse Euler Angles: Solving for q4,q5,q6

syms q4 q5 q6

%defined variables
c4 = cos(q4);
c5 = cos(q5);
c6 = cos(q6);

s4 = sin(q4);
s5 = sin(q5);
s6 = sin(q6);

Rz_q4 = [c4   -s4   0
         s4    c4   0
          0     0   1];

Ry_q5 = [c5   0   s5
          0   1    0
        -s5   0   c5];
    
Rx_q6 = [1   0    0
         0  c6  -s6
         0  s6   c6];
     
R63 = Rz_q4*Ry_q5*Rx_q6






