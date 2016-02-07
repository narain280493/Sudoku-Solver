function draw(B,n)

nn = n^2;
if (n == 2)
    k=2 ;
    ycord = 4.5;
else
    k=4;
    ycord = 9.5;
end
% k differs based on the type of solver. Got the value of k by trial and
% error.

figure;
hold on;
axis off;
axis equal

%Outside frame
rectangle('Position',[0 0 nn nn],'LineWidth',3) 

% External dark lines
rectangle('Position',[n,0,n,nn],'LineWidth',2) 
rectangle('Position',[0,n,nn,n],'LineWidth',2) 

% Internal Horizontal Lines
rectangle('Position',[0,1,nn,1],'LineWidth',1)
rectangle('Position',[0,k,nn,1],'LineWidth',1)
rectangle('Position',[0,nn-(n-1),nn,1],'LineWidth',1)

%Internal Vertical Lines
rectangle('Position',[1,0,1,nn],'LineWidth',1)
rectangle('Position',[k,0,1,nn],'LineWidth',1)
rectangle('Position',[nn-(n-1),0,1,nn],'LineWidth',1)

if size(B,2) == nn % nn columns
    %disp('inside')
    [SM,SN] = meshgrid(1:nn); % make i,j entries %%
    B = [SN(:),SM(:),B(:)];% i,j,k rows
end

for ii = 1:size(B,1)
    text(B(ii,2)-0.5,ycord-B(ii,1),num2str(B(ii,3)))
end

hold off

end
