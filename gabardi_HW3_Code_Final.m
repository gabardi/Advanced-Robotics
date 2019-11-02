%% ME 739 Homework #3
% Kaitlyn Gabardi
%03/31/2019

%% Problem 1 

%% Part 1: Find the homogenous transformation matrices symbolically

% Homogenous transformation matrices
syms L  
q1 = 90;
q2 = 90;

%defined variables
c1 = cosd(q1);
c2 = cosd(q2);

s1 = sind(q1);
s2 = sind(q2);


T10 = [c1  0   s1   L*c1
       s1  0  -c1   L*s1
       0   1    0    0
       0   0    0    1];
   
T21 = [c2  -s2   0  L*c2
       s2   c2   0  L*s2
        0    0   1    0
        0    0   0    1];
    
Tc10 = [c1   0   s1   .5*L*c1
        s1   0  -c1   .5*L*s1
        0    1    0      0
        0    0    0      1];
   
Tc21 = [c2   -s2    0   .5*L*c2
        s2    c2    0   .5*L*s2
         0     0    1      0
         0     0    0      1];
    
T20 = T10*T21
Tc20 = T10*Tc21

%% Part 2: Verify matrices evaluated above are correct

T10_ver = [0  0  1  0
           1  0  0  L
           0  1  0  0
           0  0  0  1];
       
T21_ver = [0  -1  0  0
           1   0  0  L
           0   0  1  0
           0   0  0  1];
       
T20_ver = T10_ver*T21_ver

if T20 - T20_ver == 0
   T20 - T20_ver
   fprintf('verified matched matrices\n')
end 

%% Part 3: Evaluation of the mass matri f the manipulator 

% Homogenous transformation matrices
syms L q1 q2 

%defined variables
c1 = cos(q1);
c2 = cos(q2);

s1 = sin(q1);
s2 = sin(q2);


T10_sym = [c1  0   s1   L*c1
           s1  0  -c1   L*s1
           0   1    0    0
           0   0    0    1];
   
T21_sym = [c2  -s2   0  L*c2
           s2   c2   0  L*s2
            0    0   1    0
            0    0   0    1];
    
Tc10_sym = [c1   0   s1   .5*L*c1
            s1   0  -c1   .5*L*s1
             0    1    0      0
             0    0    0      1];
   
Tc21_sym = [c2   -s2    0   .5*L*c2
            s2    c2    0   .5*L*s2
             0     0    1      0
             0     0    0      1];
    
T20_sym = T10_sym*T21_sym
Tc20_sym = T10_sym*Tc21_sym

%Computation of the linear Jacobian components
z00 = [0 
       0 
       1];
   
z10 = [s1
      -c1
        0];
    
o10 = [L*c1
       L*s1
        0];
    
oc10 = [.5*L*c1
        .5*L*s1
           0];
    
oc20 = [L*c1*(1 + .5*c2)
        L*s1*(1 + .5*c2)
            L*.5*s2];
       
%Center of mass for each link
%Jvc1

Jvc1_cross = cross(z00,oc10);

Jvc1 = [-.5*L*s1   0
         .5*L*c1   0
            0      0];
        
%Jvc2
Jvc2_cross_1 = cross(z00,oc20);
Jvc2_cross_2 = cross(z10,(oc20-o10));
Jvc2_cross_2_simp = simplify(Jvc2_cross_2);

Jvc2 = [Jvc2_cross_1 Jvc2_cross_2_simp]

%Angular velocity
Jw1 = [0  0
       0  0
       1  0];
   
Jw2 = [0   s1
       0  -c1
       1    0];
   
%Verification of Matrices 
%center of link matrix of link 1 and link 2
center_link = [-.5*L
                  0
                  0
                  1]; %center of mass of link 1&2 to origin of frame

%Center of Link 1
center_link_1 = T10_sym*center_link

%Center of Link 2
center_link_2 = T20_sym*center_link

%Used this information to vertify Jacobians matched (see analytical work)

%% Part 4: Evaluate the mass matrix D(q) in terms of mass and geometric properties
syms q1 q2
%defined variables
c1 = cos(q1);
c2 = cos(q2);

s1 = sin(q1);
s2 = sin(q2);


R10 = [c1  0   s1 
       s1  0  -c1  
       0   1    0];

R20 = [c1*c2   -c1*s2    s1
       c2*s1   -s1*s2   -c1
         s2       c2      0];
     
syms L m Ia

z00 = [0 
       0 
       1];
   
z10 = [s1
      -c1
        0];
    
o10 = [L*c1
       L*s1
        0];
    
oc10 = [.5*L*c1
        .5*L*s1
           0];
    
oc20 = [L*c1*(1 + .5*c2)
        L*s1*(1 + .5*c2)
            L*.5*s2];
       
%Center of mass for each link
%Jvc1

Jvc1_cross = cross(z00,oc10);

Jvc1 = [-.5*L*s1   0
         .5*L*c1   0
            0      0];
        
%Jvc2
Jvc2_cross_1 = cross(z00,oc20);
Jvc2_cross_2 = cross(z10,(oc20-o10));
Jvc2_cross_2_simp = simplify(Jvc2_cross_2);

Jvc2 = [Jvc2_cross_1 Jvc2_cross_2_simp];

%Angular velocity
Jw1 = [0  0
       0  0
       1  0];
   
Jw2 = [0   s1
       0  -c1
       1    0];
   
Ia = [.5*Ia   0    0
        0    Ia    0
        0     0   Ia];

d1_q = simplify((m*transpose(Jvc1)*Jvc1 + transpose(Jw1)*R10*Ia*transpose(R10)*Jw1));

d2_q = simplify((m*transpose(Jvc2)*Jvc2 + transpose(Jw2)*R20*Ia*transpose(R20)*Jw2));

D_q = simplify(d1_q + d2_q)


%% Part 5: Evaluate the centrifugal and Coriolis inertial terms
syms q1 q2 Ia L m

%defined variables
c1 = cos(q1);
c2 = cos(q2);

s1 = sin(q1);
s2 = sin(q2);

% Mass Matrix

D = [(3*Ia)/2 + (5*L^2*m)/4 + (Ia*cos(q2)^2)/2 + L^2*m*cos(q2) + (L^2*m*cos(q2)^2)/4,              0
                                                                                     0, (m*L^2)/4 + Ia];

%Evaluate patrial derivatives of mass matrix
d111 = diff(D(1,1),q1);
d112 = diff(D(1,1),q2);
d121 = diff(D(1,2),q1);
d122 = diff(D(1,2),q2);
d211 = diff(D(2,1),q1);
d212 = diff(D(2,1),q2);
d221 = diff(D(2,2),q1);
d222 = diff(D(2,2),q2);

% evaluate the Christoffel symbols
b112 = (1/2)*(d112 + d121 - d121);
b212 = (1/2)*(d212 + d221 - d122);
b111 = (1/2)*(d111 + d111 - d111);
b122 = (1/2)*(d122 + d122 - d221);
b211 = (1/2)*(d211 + d211 - d112);
b222 = (1/2)*(d222 + d222 - d222);

% form the Coriolis and centrifugal matrices
B = [2*b112; 2*b212];

C = [b111 b122
     b211 b222];

%Find B(q)[qdot qdot]: Coriolis Terms
pretty(B)

%Find C(q)[qdot^2]: Centrifugal Terms
pretty(C)

                   
% NOTE: This matches my analytical work, although in code, cos(q2)*sin(q2)
%       represents 1/2sin(q2), therefore then matches results

%% Part 6: Evalutate the gravity vector, G(q). 
syms L m Ia q1 q2 g

%In frame {0} the gravity vector is given as:
g = [g; 0; 0];

%defined variables
c1 = cos(q1);
c2 = cos(q2);

s1 = sin(q1);
s2 = sin(q2);


R10 = [c1  0   s1 
       s1  0  -c1  
       0   1    0];

R20 = [c1*c2   -c1*s2    s1
       c2*s1   -s1*s2   -c1
         s2       c2      0];

z00 = [0 
       0 
       1];
   
z10 = [s1
      -c1
        0];
    
o10 = [L*c1
       L*s1
        0];
    
oc10 = [.5*L*c1
        .5*L*s1
           0];
    
oc20 = [L*c1*(1 + .5*c2)
        L*s1*(1 + .5*c2)
            L*.5*s2];
       
%Center of mass for each link
%Jvc1

Jvc1_cross = cross(z00,oc10);

Jvc1 = [-.5*L*s1   0
         .5*L*c1   0
            0      0];
        
%Jvc2
Jvc2_cross_1 = cross(z00,oc20);
Jvc2_cross_2 = cross(z10,(oc20-o10));
Jvc2_cross_2_simp = simplify(Jvc2_cross_2);

Jvc2 = [Jvc2_cross_1 Jvc2_cross_2_simp];

G1 = [transpose(Jvc1)*m*g];
G2 = [transpose(Jvc2)*m*g];

G_total = -[G1 + G2];
simplify(G_total)

%% Part 7: Form the complete equations of motion 
 %Wrote out complete form in analytical work
 
%% Problem 2

%% Part 1: Write a MATLAB numerical simulation of the manipulator
%Implement a 4th order Runge-Kutta numerical integration of EOM

%Geometric parameters
%L = 2;    %meters
%m = 10;   %kg
%Ia = 5;   %km/m^s
%g = 9.81; %m/s^2

%Simulation parameters
%delta_t = 0.01; %integration time step in seconds
%t_final = 20;   %simulated time duration

% model parameters 
modelParameters.g = 9.81; % gravitational constant [m/s^2] 
modelParameters.m = 10;   % link mass [kg] 
modelParameters.L = 2;    % link length [m] 
modelParameters.Ia = 5;   % inertia [kg/m^2] 
 
%%%%%%%%%%%%%%%%%% initialize integration variables %%%%%%%%%%%%%%%%%%%%%%%
dT = .01;                 %integration step size 
tend = 20;                %simulation run time 
numPts = floor(tend/dT); 
q = zeros(2,numPts);      %pre-allocate array 
dq = zeros(2,numPts); 
t = zeros(1,numPts); 
q(:,1) = [pi; 0];         %initial position when q1=180 and q2=0
qd(:,1) = [0.1; 0.1];     %initial joint velocity of 0.1 r/s

%initialize the state variables 
z = [q(:,1); qd(:,1)];  %initialize the state variables 
 
%%%%%%%%%%%%%%%%%%% integrate equations of motion %%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numPts-1 
    % Runge-Kutta 4th order 
    k1 = zDot2dof(z,modelParameters); 
    k2 = zDot2dof(z + 0.5*k1*dT,modelParameters); 
    k3 = zDot2dof(z + 0.5*k2*dT,modelParameters); 
    k4 = zDot2dof(z + k3*dT,modelParameters); 
    z  = z + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*dT; 
 
    % store joint position and velocity for post processing 
    q(:,i+1)  = z(1:2); 
    qd(:,i+1) = z(3:4); 
    t(1,i+1)  = t(1,i) + dT; 
end 
 
%%%%%%%%%%%%%%%%%% setting rendering window view parameters %%%%%%%%%%%%%%%

L = modelParameters.L; 
L1 = L; 
L2 = L; 
f_handle = 1; 
axis_limits = L*[-2 2 -2 2 -1 1]; 
render_view = [-1 -1 1]; 
view_up = [-1 0 0];
SetRenderingViewParameters(axis_limits,render_view,view_up,f_handle); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%% initialize rendering %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% link 1 rendering initialization 
r1 = L1/8; sides1 = 10; axis1 = 1; norm_L1 = 1.0; 
linkColor1 = [0 0.75 0]; plotFrame1 = 0; 
d1 =  CreateLinkRendering(L1,r1,sides1,axis1,norm_L1,linkColor1,... 
    plotFrame1,f_handle); 
 
% link 2 rendering initialization 
r2 = L2/8; sides2 = 10; axis2 = 1; norm_L2 = 1.0; 
linkColor2 = [0.75 0 0]; plotFrame2 = 0; 
d2 =  CreateLinkRendering(L2,r2,sides2,axis2,norm_L2,linkColor2,... 
    plotFrame2,f_handle); 
 
for i = 1:numPts 
    % Update frame {1} 
    c = cos(q(1,i)); 
    s = sin(q(1,i)); 
    L = L1; 
    
    T10 = [c  0   s  L*c 
           s  0  -c  L*s 
           0  1   0  0 
           0  0   0  1]; 
 
    % Update frame {2} 
    c = cos(q(2,i)); 
    s = sin(q(2,i)); 
    L = L2; 
    
    T21 = [c  -s  0  L*c 
           s   c  0  L*s 
           0   0  1   0 
           0   0  0   1]; 
       
    T20 = T10*T21; 
 
    UpdateLink(d1,T10); 
    UpdateLink(d2,T20); 
 
    if i == 1  %pause at start of simulation rendering 
        pause; 
    end 
    pause(dT); 
end 
 
%% Part 2: System Energy Plots of PE, KE and Total Energy

for i = 1:numPts 
    % kinetic energy 
    Q = q(:,i); 
    Qd = qd(:,i); 
    D = Dmatrix_2DOF(Q,modelParameters); 
    T(i) = (1/2)*Qd'*D*Qd; 
 
    % potential energy 
    % Link 1: 
    c = cos(q(1,i)); 
    s = sin(q(1,i)); 
    L = L1; 
    
    T10 = [c  0    s   L*c 
           s  0   -c   L*s 
           0  1    0    0 
           0  0    0    1]; 
        
    L = L1/2; 
    T10c = [c   0    s   L*c 
            s   0   -c   L*s 
            0   1    0    0 
            0   0    0    1]; 
 
    % Link 2: 
    c = cos(q(2,i)); 
    s = sin(q(2,i)); 
    L = L2; 
    T21  = [c  -s   0  L*c 
            s   c   0  L*s 
            0   0   1   0 
            0   0   0   1]; 
       
    L = L2/2; 
    T21c = [c  -s  0  L*c 
            s   c  0  L*s 
            0   0  1   0 
            0   0  0   1]; 
 
    T20 = T10*T21; 
    T20c = T10*T21c; 
 
    % assign center of mass position vectors 
    rc1 = T10c(1:3,4); 
    rc2 = T20c(1:3,4); 
    g = modelParameters.g*[1; 0; 0]; 
    m = modelParameters.m; 
 
    % calculate the gravitational potential energy 
    V(i) = -(m*g'*rc1 + m*g'*rc2); 
 
    % calculate the total system energy 
    E(i) = V(i) + T(i); 
end 
 
%Create plot that shows KE, PE, and total energy as a function of time 
figure; 
plot(t,E,'b',t,V,'g',t,T,'r'); 
xlabel('Time (s)');
ylabel('Energy (J)'); 
legend('Total Energy', 'Potential Energy', 'Kinetic Energy','location','SouthEast') 
title('System Energy as a Function of Time'); 
axis([0 20 -1000 1000])


 
 