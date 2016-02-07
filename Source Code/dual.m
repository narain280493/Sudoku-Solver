function [dualA, dualb, dualc, dualbas, dualnonbas] = dual(A,b,c,bas,nonbas)
    dualA = -(A)';
    dualb = -(c);
    dualc = -(b);
    [dualnonbas, dualbas] = deal(bas, nonbas); 
end