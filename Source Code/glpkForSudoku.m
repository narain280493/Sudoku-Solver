function [xopt, fmin, status, extra] = glpkForSudoku(f,Aeq,beq,lb,ub,NNN)
    size(Aeq)
    param.msglev = 1;
    cType = repmat('S',1,size(Aeq,1));
    vType = repmat('I',1,NNN);
    lb = reshape(lb,NNN,1);
    ub = reshape(ub,NNN,1);
    [xopt, fmin, status, extra] = glpk(f,Aeq,beq,lb,ub,cType,vType,-1,param);
end
