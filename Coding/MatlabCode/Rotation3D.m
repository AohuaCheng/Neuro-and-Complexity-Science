function [NewC]=Rotation3D(C,T,Origin)
% T defines rotation angles
% C is the original coordinate
% Origin is the origin coordinate of rotation

%% Define the rotation matrix
T1=T(1);
T2=T(2);
T3=T(3); %% Here T1,T2, and T3 are rotation angles
Rx = [1,0,0; 0,cos(T1),-sin(T1); 0,sin(T1),cos(T1)]; %% Rotation in X direction
Ry = [cos(T2),0,sin(T2); 0,1,0; -sin(T2),0,cos(T2)]; %% Rotation in Y direction
Rz = [cos(T3),-sin(T3),0; sin(T3),cos(T3),0; 0,0,1]; %% Rotation in Z direction
R=Rx*Ry*Rz; %% Rotation matrix

% A=T(1);
% B=T(2);
% C=T(1);
% R=[cos(A)*cos(B), cos(A)*sin(B)*sin(C)-sin(A)*cos(C), cos(A)*sin(B)*cos(C)+sin(A)*sin(C); sin(A)*cos(B), sin(A)*sin(B)*sin(C)+cos(A)*cos(C), sin(A)*sin(B)*cos(C)-cos(A)*sin(C); -sin(B), cos(B)*sin(C), cos(B)*cos(C)];

%% Transform the coordinate
NewC=C*R+Origin;