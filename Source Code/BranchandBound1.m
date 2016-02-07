function [FinalOptimalsolution,Finalobjval,FinalExitFlag]=BranchandBound(c,A,b,Aeq,beq,lb,ub,Integers)
bound=inf;
[InitialOptimalSolution,Intobjval]=linprog(c,A,b,Aeq,beq,lb,ub,[]); 
[FinalOptimalsolution,Finalobjval,FinalExitFlag,bb]=Branch(c,A,b,Aeq,beq,lb,ub,InitialOptimalSolution,Intobjval,Integers,bound);