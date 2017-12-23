function Hhat=ConstrainHeights(Nx,Nobs,Hobs,iobs)

%iobs == index of x measured by each observation (j=1...Nobs)

f=[zeros(Nx,1)
   ones(Nobs,1)];

J=zeros(Nobs,Nx);
for i=1:Nobs,
    J(i,iobs(i))=1;
end

F=zeros(Nx-1,Nx);
F(1:Nx-1,1:Nx-1)=F(1:Nx-1,1:Nx-1)+eye(Nx-1);
F(1:Nx-1,2:Nx)=F(1:Nx-1,2:Nx)-eye(Nx-1);

B=[ J -eye(Nobs);
   -J -eye(Nobs);
   F  zeros(Nx-1,Nobs);];

c=[Hobs;
  -Hobs;
  zeros(Nx-1,1);];

% zguess=[zeros(Nx,1)
%         Hguess'];
    
zhat=linprog(f,B,c);

Hhat=zhat(1:Nx);

return