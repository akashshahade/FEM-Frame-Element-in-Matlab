
%-------------------------------------------------------------------%
% Software by -  Akash S. Shahade
%
% Title - Program for FE Analysis of Frame.
%-------------------------------------------------------------------%

clc;clear;
fprintf(' \nFEM ANALYSIS OF FRAME. \n\n');

%------------------------------------------------------%
% STEP 01 - PRE-PROCESSING
%------------------------------------------------------%

area = 4*10^(-3);          % Area (m^2)
moi = 2.133*10^(-6);       % Moment of Inertia (m^4)
E = 70*10^9;               % Young's Modulus (kN/m^2)

ae = area*E;
ei=E*moi;

L(1) = 2;            % Length of element 01 (m)
L(2) = 1;            % Length of element 02 (m)

ele_nod=[1 2;2 3];                   % Elements are connected with these nodes
nod_coor=[0 0;0 2;1 2];              % Node coordinates
num_ele=size(ele_nod,1);             % Number of elements
ele_dof=[1 2 3 4 5 6;4 5 6 7 8 9];   % D.O.F  associated with Nodes [1 2;2 3]
num_nod=3;                           % Number of Nodes
dof = 3;                             % D.O.F per node

displacement = zeros(dof*num_nod,1);    % Zero Matrix for Displacement
force = zeros(dof*num_nod,1);           % Zero Matrix for Force
stiffness = zeros(dof*num_nod);         % Zero Matrix for Stiffness

t1 = [0 1 0 0 0 0;...                   % Transformation Matrix for Element 01
    -1 0 0 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 0 1 0;...
    0 0 0 -1 0 0;...
    0 0 0 0 0 1];

t2 = [1 0 0 0 0 0;...                   % Transformation Matrix for Element 02
    0 1 0 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1];

%------------------------------------------------------%
% Stiffness matrix calculation & ASSEMBLY
%------------------------------------------------------%

for e=1:num_ele                               % For 1 to Number of elements
    
    ke = [ (ae/L(e)) 0 0 (-ae/L(e)) 0 0;...
        0 ((12*ei)/(L(e))^3) ((6*ei)/(L(e)^2)) 0 ((12*ei)/(L(e))^3) ((6*ei)/(L(e)^2));...
        0 ((6*ei)/(L(e)^2)) (4*ei/L(e)) 0 ((-6*ei)/(L(e)^2)) (2*ei/L(e));...
        (-ae/L(e)) 0 0 (ae/L(e)) 0 0;...
        0 ((-12*ei)/(L(e))^3) ((-6*ei)/(L(e)^2)) 0 ((12*ei)/(L(e))^3) ((-6*ei)/(L(e)^2));...
        0 ((6*ei)/(L(e)^2)) (2*ei/L(e)) 0 ((-6*ei)/(L(e)^2)) (4*ei/L(e))];

    if e == 1
        k = t1.'*ke*t1;                    % Multiplication with Transformation Matrix
    else
        k = t2.'*ke*t2;
    end
    
    
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);

    for i=1:6
        for j=1:6
                                              % Assembly of Global Matrix
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);

        end
    end
end

force(4) = 17320;
force(5) = 10000;
force(9) = -15000;

fprintf('Global Stiffness Matrix: \n');
disp(stiffness);
fprintf('\n Global Load Vector: \n');
disp(force);
fprintf('\n------------------------------------------------\n');


%------------------------------------------------------%
% Boundary Conditions
%------------------------------------------------------%

fixed_dof = [1 2 3 8];               % Constrained D.O.F.
k=stiffness;
k(fixed_dof,:)=[];                   % Eliminating Rows
k(:,fixed_dof)=[];                   % Eliminating Columns

f=force;
f(fixed_dof,:)=[];                   % Eliminating Rows

%------------------------------------------------------%
% STEP 02 - SOLVE
%------------------------------------------------------%

q = k\f ;

%------------------------------------------------------%
% STEP 03 - POST-PROCESSING
%------------------------------------------------------%

displacement=[0;0;0;q;0];                 % Displacement Vector

reaction = stiffness*displacement - force;  % Calculate Reaction Forces

Node = [1;1;1;2;2;2;3;3;3];
qxy = {'q1x';'q1y';'Theta1';'q2x';'q2y';'Theta2';'q3x';'q3y';'Theta3'};
Q = displacement;
T=table(Node,qxy,Q);
disp(T);

fprintf('\n------------------------------------------------\n');

F = {'R1x';'R1y';'M1';'R2x';'R2y';'M2';'R3x';'R3y';'M3'};
Reaction = reaction;
J = table(F,Reaction);
disp(J);
fprintf('END OF PROGRAM.\n');
