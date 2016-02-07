function [A,b,c,bas,nonbas,isDone] = nextDic(bas,nonbas,A,b,c)

    pos = find(bas <= length(nonbas));
    check = b(pos);
 
    check = single(check);
    nonintvalues = find((check - floor(check))>1E-08, 1);
    if(isempty(nonintvalues))
        isDone = 1;
    else
        isDone = 0;
        pos = randsample(pos,1);
        temp = (floor(A(pos,:)) - A(pos,:));
        if(isempty(find(abs(temp) > 0,1)))
            nextDic(bas,nonbas,A,b,c);
        elseif(sum(ismember(A, temp, 'rows')))
            nextDic(bas,nonbas,A,b,c);    
        else
            A = [A; temp];
            b = [b; (floor(b(pos,1)) - b(pos,1))];
            bas = [bas; (length(bas) + length(nonbas) + 1)];
        end
    end
end

