function [sudokuPuzzle,statusFlag] = ruleModel(clueMatrix,sz,n,choice)

globalstarttime = cputime;
%sz = 2;% for a 9*9 sudoku 
sz2 = sz ^ 2;
draw(clueMatrix,sz) %representing the partially filled sudoku puzzle

[a,b] = size(clueMatrix);
if b~=3 % The matrix dimensions must be atleast n*3. Atleast 17 clues are needed at minimum, to solve the sudoku uniquely.
    error('Incorrect clue matrix dimensions. The matrix dimensions must be N*3')
end

%% Constraint 1
% The numbers in the clue matrix must lie between 1 and 9. The numbers must also be an integer!
if sum([any(rem(clueMatrix,1)~=0),any(clueMatrix < 1),any(clueMatrix > 9)])
    error('Error: Numbers should be between 1 and 9. The numbers must be integers! ')
end

%% The rules of Sudoku:

% number of variables in the LP form. 
% Eg. 729 variables for a 9*9 sudoku puzzle. Since we are converting a 9*9
% puzzle into a 9*9*9 matrix.
N = n^3;

% Total number of constraints. There are n^2 constraints for one rule.
% Since there are 4 main rules in sudoku, we get 4*n^2 constraints totally.
M = 4*n^2;

% Initial the Aeq matrix
Aeq = zeros(M,N);

% A non-constant objective function helps speed up the solver.
%f = (1:N)';
%size(f)
f = ones(N,1);

% Setting the lower and upper bounds as 0 and 1 as we are considering a
% binary integer linear programming model where each variable lies between
% 0 and 1
lb = zeros(n,n,n); 
ub = lb+1; 

%% Constraint 2: Rule 1

% in each row i of the 2-D matrix, there is exactly one value out of each of
% the digits from 1 to 9.

iterator = 1;
j=1;
while(j <= n)
    k =1;
    while(k <= n)
        Atemp = zeros(n,n,n);
        for i =1:n
            Atemp(i,j,k) = 1;
        end
        Aeq(iterator,:) = Atemp(:)';
        iterator = iterator + 1;
        k = k+1;
    end
    j = j+1;
end

%% Constraint 3: Rule 2

% in each column j of the 2-D matrix, there is exactly one value out of each of
% the digits from 1 to 9.

i=1;
while(i <= n)
    k =1;
    while(k <= n)
        Atemp = lb;
        for j=1:n
            Atemp(i,j,k) = 1;
        end
        Aeq(iterator,:) = Atemp(:)';
        iterator= iterator + 1;
        k = k+1;
    end
    i = i+1;
end

%% Constraint 4: Rule 3

% in each depth k of the 2-D matrix, there is exactly one value out of each of
% the digits from 1 to 9.

i=1;
while(i <= n)
    j =1;
    while(j <= n)
        Atemp = lb;
        for k=1:n
            Atemp(i,j,k) = 1;
        end
        Aeq(iterator,:) = Atemp(:)';
        iterator= iterator + 1;
        j = j+1;
    end
    i = i+1;
end
%% Constraint 5: Rule 4

% in each square (0,3,9) or (0,2,4) of the matrix, there is exactly one value out of each of
% the digits from 1 to 9.

U=0;
while U < sz2
    V=0;
    while V < sz2
        for k=1:n
            Atemp = lb;
            Atemp(U+(1:sz),V+(1:sz),k) = 1;
            Aeq(iterator,:) = Atemp(:)';
            iterator = iterator + 1;
        end
        V = V + sz;
    end
    U = U + sz;
end

%% Constraint 6: 

%Getting the size of the cluematrix and adding that to the size of the Aeq
%and beq matrix. 
clueSize=size(clueMatrix,1);

% Keeping the lower bound of the indexes of the clues as 1. Upper bound is
% also 1. So we are indirectly fixing the clues in the main puzzle.

for i = 1:size(clueMatrix,1) % size(B,1) returns number of rows
    lb(clueMatrix(i,1),clueMatrix(i,2),clueMatrix(i,3)) = 1;
end

% All N integers must be integers.
intcon = 1:N;

if(choice == 1)
   disp('Using intlinprog and L1-norm')
   
% The clue matrix can be included as constraints in the Aeq and beq matrix.
% This is because those particular indexes MUST be 1. This can also be done
% by simply putting the lower bound of those indexes as 1.

   for i=1:clueSize
    Aeq(4*n^2+i,clueMatrix(i,1)+(clueMatrix(i,2)-1)*n+(clueMatrix(i,3)-1)*n^2)=1;
   end
   beq = ones(M+clueSize,1); 
   [x1,res1,flag,fval]=L1norm(Aeq,beq,lb,ub,f,intcon);
   fval
end
%  Setting the options in intlinprog.

%% All parameters are available here. Call Branch and Bound from here.
%% Calling intlinprog solver. The rules are represented in the Aeq
% and beq matrices, and the clues are ones in the lb array.

if(choice == 2)
    disp('Branch and bound using linprog')
   % options = optimoptions('intlinprog','CutGeneration','none','BranchingRule','maxpscost','Display','off');
    starttime = cputime;
    %[x1,fval,flag] = intlinprog(f,intcon,[],[],Aeq,beq,lb,ub,options);
    beq = ones(M,1); 
    [x1,fval,flag] = BranchandBound1(f,[],[],Aeq,beq,lb,ub,intcon);
    stoptime = cputime;
    runningtime = stoptime - starttime
    fval %The optimal objective function value. Not needed in the case of a sudoku.
end
 
if(choice == 3)
    disp('Using cutting plane method')
    cuttingPlaneForSudoku(f,Aeq,beq,lb,ub,729)
end
%{
if(choice == 4)
   starttime = cputime;
  [x1, fval,flag, extra] = glpkForSudoku(f,Aeq,beq,lb,ub,729)
  stoptime = cputime;
  runningtime = stoptime - starttime
  fval
end
%}
%% The solution is a vector of 729*1. It is converted into a displayable form.
%{
% If we get flag =1, it means we have an optimal solution.
if flag ==1 
    %Currently x is a 729*1 vector. We convert it back to a 9*9*9 vector
    %for our convenience.
    x = reshape(x,n,n,n);  
    depth = ones(size(x));
    
    % multiplier for each depth k
    for k = 2:n
        depth(:,:,k) = k; 
    end
    %y
    % multiply each entry by its depth
    SudokuPuzzle = x.*depth;
    %SudokuPuzzle
    answer = zeros(n,n);
   
    % Compressing the 3D matrix to a 2D matrix.
    for i=1:n
        for j=1:n
            for k=1:n
                if(SudokuPuzzle(i,j,k)~=0)
                    answer(i,j)=SudokuPuzzle(i,j,k);
                end
            end
        end
    end
    answer        
    %SudokuPuzzle = sum(SudokuPuzzle,2); 
    %SudokuPuzzle
    
% Solved puzzle.
    draw(answer,sz)
    globalstoptime = cputime;
    globalrunningtime = globalstoptime - globalstarttime
    
else
    SudokuPuzzle = [];
    error('The sudoku puzzle could not be solved using the given clues')
end
%}
%size(x1)
if flag ==1 || stat == 5
    %Currently x is a 729*1 vector. We convert it back to a 9*9*9 vector
    %for our convenience.
    x1 = reshape(x1,n,n,n);  
    depth = ones(size(x1));
    %size(x1)
    % multiplier for each depth k
    for k = 2:n
        depth(:,:,k) = k; 
    end
    %y
    % multiply each entry by its depth
    SudokuPuzzle = x1.*depth;
    %SudokuPuzzle
    answer = zeros(n,n);
   
    % Compressing the 3D matrix to a 2D matrix.
    for i=1:n
        for j=1:n
            for k=1:n
                if(SudokuPuzzle(i,j,k)~=0)
                    answer(i,j)=SudokuPuzzle(i,j,k);
                end
            end
        end
    end
   
% Solved puzzle.
    draw(answer,sz)
    globalstoptime = cputime;
    globalrunningtime = globalstoptime - globalstarttime
    
else
    SudokuPuzzle = [];
    error('The sudoku puzzle could not be solved using the given clues')
end