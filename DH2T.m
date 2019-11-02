function [T] = DH2T(d,a,alpha,theta)
%DH2T function takes in a matrix of relationship between links and DH
%parameters 
 
% Build up Transformation Matrix
% Initialization
T = eye(4);

 T = [cos(theta)     -sin(theta)*cos(alpha)     sin(theta)*sin(alpha)    a*cos(theta)
      sin(theta)      cos(theta)*cos(alpha)    -cos(theta)*sin(alpha)    a*sin(theta)
          0                sin(alpha)                cos(alpha)              d
          0                    0                         0                   1];

end %end of function