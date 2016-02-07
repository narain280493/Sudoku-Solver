function [FinalOptimalSolution,Finalobjval,FinalEXITFLAG,bb]=Branch(f,A,b,Aeq,beq,lb,ub,InitialOptimalSolution,Initialobjval,Integers,bound) 
[Optimalx,objval,EXITFLAG]=linprog(f,A,b,Aeq,beq,lb,ub,[]); 

Infeasible=-2
Unbounded = -3
Nanval = - 4
PDInfeasible = -5
SmallMag = -7

Optimalflag=1

if ismember(EXITFLAG,Infeasible)==1 | ismember(EXITFLAG,Unbounded)==1 | ismember(EXITFLAG,Nanval)==1 | ismember(EXITFLAG,PDInfeasible)==1 | ismember(EXITFLAG,SmallMag)==1 |  objval > bound  
    FinalOptimalSolution=InitialOptimalSolution; Finalobjval=Initialobjval; FinalEXITFLAG=EXITFLAG; bb=bound;
    return;
end


residual = (Optimalx(Integers)-round(Optimalx(Integers)))
absoluteresidual = abs(residual)
resindex=[]

j=1
for i=1:numel(absoluteresidual)
    if absoluteresidual(i) > 1e-10
        resindex(j)= i
        j=j+1
        
    end;
end; 

[x y]=size(resindex)

if x==0 && y==0
    FinalEXITFLAG=1;        
    if objval < bound   
        Optimalx(Integers)=round(Optimalx(Integers));
        FinalOptimalSolution=Optimalx;        
        Finalobjval=objval;
        bb=objval;
    else
        FinalOptimalSolution=InitialOptimalSolution; 
        Finalobjval=Initialobjval;
        bb=bound;
    end
    return
end
[u v]=size(A)
if u==0 && v==0
    [r s]=size(Aeq);
else
    r=u;
    s=v;
end

branchvariable=Integers(resindex(1));
branchvaluef=floor(InitialOptimalSolution(Integers(resindex(1))));
branchvaluec=-ceil(InitialOptimalSolution(Integers(resindex(1))));

% [f g]= size(A)
% k=1
% A1=A
% while k<=s
%     A1(f+1,k)=0
%      k=k+1
%  end 
A1=[A ; zeros(1,s)];
[j k]=size(A1);
k=branchvariable
A1(j,k)=1;
b1=[b;branchvaluef];

[u v]=size(A)
if u==0 && v==0
    [r s]=size(Aeq);
else
    r=u;
    s=v;
end
% [f g]= size(A)
%  k=1
%  A2=A
% while k<=s
%     A2(f+1,k)=0
%      k=k+1
%  end 
A2=[A ;zeros(1,s)];
[j k]=size(A2);
k=branchvariable
A2(j,k)=-1;
b2=[b; branchvaluec];

[Optimalx1,objval1,EXITFLAG1,bound1]=Branch(f,A1,b1,Aeq,beq,lb,ub,Optimalx,objval,Integers,bound);
FinalEXITFLAG=EXITFLAG1;
if EXITFLAG1 >0 & bound1<bound 
   FinalOptimalSolution=Optimalx1;
   Finalobjval=objval1;
   bound=bound1;
   bb=bound1;
else
    FinalOptimalSolution=Optimalx;
    Finalobjval=objval;
    bb=bound;
end
    
[Optimalx2,objval2,EXITFLAG2,bound2]=Branch(f,A2,b2,Aeq,beq,lb,ub,Optimalx,objval,Integers,bound);

if EXITFLAG2 >0 & bound2<bound 
    FinalEXITFLAG=EXITFLAG2;
    FinalOptimalSolution=Optimalx2;
    Finalobjval=objval2;
    bb=bound2;
end