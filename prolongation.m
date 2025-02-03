function [uoutnew]=prolongation(uout)
[Ny,Nx]=size(uout);
Nxnew=(Nx-1)*2+1;
Nynew=(Ny-1)*2+1;

%Prolongation step
uoutnew(1:Nynew,1:Nxnew)=0;
uoutnew(1:Nynew,1:Nxnew)=interp2(uout);

    
end

