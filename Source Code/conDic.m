function [A,b,c,bas,nonbas] = conDic(bas,nonbas,A,b,c,isDual)
    Ab = A(:,bas);
    Ai = A(:,nonbas);
    cb = c(bas,:)';
    ci = c(nonbas,:)';
    pII = (Ab'\cb')';
    c = (ci - pII*Ai)';
    A = (Ab\Ai);
    b = Ab\b;
    if(isDual)
        [A,b,c,bas,nonbas] = dual(A,b,c,bas,nonbas);
    end
end

