function [X1,res1,flag,fval] = L1norm(A,b,lb,ub,f,intcon)
%% Solves least squares for matrices A,b using three different norms and
%% output residuals.
%% X1: L1 norm best estimate.
%% res1: L1 norm residual (AX1 - b)
%A
[m,n] = size(A);

%% Computation for L1 norm.
%% Lets prepare the matrices for solving linprog.

%%constraint
M1=[ A -eye(m); 
    -A -eye(m);
    zeros(m,n) -eye(m)];
%%RHS consts
r1=[b; -b; zeros(m,1)];
%% objective.
f1 = [f; zeros(m,1)];


localstarttime = cputime;
%% Solving for the L1 minima
options = optimoptions('intlinprog','CutGeneration','none','BranchingRule','maxpscost','Display','off');

[X11,fval,flag]= intlinprog(f1,intcon,[],[],M1,r1,lb,ub,options);
localendtime = cputime;
localrunning = localendtime - localstarttime

%% Lets drop the extra variables
X1=X11(1:n);
%% Compute residuals
res1 = A*X1- b
end