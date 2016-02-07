function [] = cuttingPlaneForSudoku(f,Aeq,beq,lb,ub,NNN)
    A = [Aeq;-Aeq;-eye(NNN);eye(NNN)];
    b = [beq;-beq;-reshape(lb,NNN,1);reshape(ub,NNN,1)];
    [m,n] = size(A);
    I = eye(m);
    A = [A I];
    f = [f;zeros(m,1)];
    nonbas = (1:n)';
    bas = (n+1:n+m)';
    isDual = 0;
    objValue = 0;
    pivotT(A,b,f,bas,nonbas,objValue,isDual);
end

