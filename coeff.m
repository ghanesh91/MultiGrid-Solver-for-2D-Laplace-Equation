function [invA,invB,A,B,a,b,c]=coeff(n,Nx,Ny,Lx,Ly)
 %Number of points at n level
 ny=((Ny-1)/2^(n-1))+1;
 nx=((Nx-1)/2^(n-1))+1;   
 dx=Lx/(nx-1);
 dy=Ly/(ny-1);
 
 %Creating matrices for Line-GS and computing its inverse.
 a=(1/dx^2);b=-((2/dx^2)+(2/dy^2));c=(1/dy^2);
 noA = nx-2;
 nOnes = ones(noA, 1) ;
 A = diag(b * nOnes, 0) + diag(c*nOnes(1:noA-1), -1) + diag(c*nOnes(1:noA-1), 1); 
 B = diag(b * nOnes, 0) + diag(a*nOnes(1:noA-1), -1) + diag(a*nOnes(1:noA-1), 1);
 invA=inv(A);
 invB=inv(B);
end
