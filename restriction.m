function [epsnew]=restriction(eps,Nx,Ny)
Nxnew=(Nx-1)/2+1;
Nynew=(Ny-1)/2+1;

epsnew(1:Nynew,1:Nxnew)=0;
%Perform restriction by copying every alternate points in the grid
epsnew(1:Nynew,1:Nxnew)=eps(1:2:Ny,1:2:Nx);
end