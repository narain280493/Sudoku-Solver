function [bas1,nonbas1,objval,isFinal,isUnbounded,val] = pivotT(A,b,c,bas,nonbas,val,isDual)

[bas1,nonbas1,objval,isFinal,isUnbounded] = pivot(A,b,c,bas,nonbas);
while (isFinal ~= 1 && isUnbounded ~= 1)
    [bas1,nonbas1,objval,isFinal,isUnbounded] = pivot(A,b,c,bas1,nonbas1);
end

if(isDual)
    objval = objval * -1;
end
val = val + objval;

[A,b,c,bas,nonbas] = conDic(sort(bas1),sort(nonbas1),A,b,c,isDual);
[A,b,c,bas,nonbas,isDone] = nextDic(bas,nonbas,A,b,c);

if(~isDone)
    if(isempty(find(c>0, 1)))
        [A, b, c, bas, nonbas] = dual(A,b,c,bas,nonbas);
        [m,n] = size(A);
        I = eye(m);
        temp = zeros(m, (m+n));
        temp(:,nonbas) = A;
        temp(:,bas) = I;
        A = temp;
        temp = zeros((m+n),1);
        temp(nonbas) = c;
        c = temp;
        isDual = 1;
    end
    pivotT(A,b,c,bas,nonbas,val,isDual);
else
    fprintf('\nFeasible Solution:');
    pos = find(bas <= length(nonbas));
    for i = 1:length(pos)
        fprintf('\nx(%i) = %+5.2f' ,pos(i),b(pos(i)));
    end
    fprintf('\n\nFeasible Value:\n%f \n' , val);
end
end