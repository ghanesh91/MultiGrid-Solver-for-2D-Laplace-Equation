function [eps]=residual(uout,RHS,a,b,c,Nx,Ny)
eps(1:Ny,1:Nx)=0;
%Compute residual \epsilon=Lu-R
for i=2:Nx-1
    for j=2:Ny-1
        eps(j,i)=a*(uout(j,i-1)+uout(j,i+1))+c*(uout(j-1,i)+uout(j+1,i))+b*uout(j,i)-RHS(j,i);
    end
end

end