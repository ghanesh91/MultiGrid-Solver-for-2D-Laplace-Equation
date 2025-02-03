function [uout]=fine_to_coarse(uin,RHS,maxiter,invA,invB,A,B,a,b,c,nx,ny,itype,omg,Nlvlmax,Nlvlmin,uout,eps,epsnew)

%Fine to coarse 
      for n=Nlvlmin:Nlvlmax
         fprintf('na=%d \n',n)
         %Solve for Lu=R iteratively "maxiter" times
         uout{n}=iterative_solve(uin{n},RHS{n},maxiter,invA{n},invB{n},A{n},B{n},a(n),b(n),c(n),nx(n),ny(n),itype,omg);
         %Compute the residual from the smoothened solution
         eps{n}=residual(uout{n},RHS{n},a(n),b(n),c(n),nx(n),ny(n));
         %Define the new RHS for the next level
         if(n~=Nlvlmax)
            epsnew{n+1}=restriction(eps{n},nx(n),ny(n));
            RHS{n+1}=-epsnew{n+1};
            uin{n+1}(1:ny(n+1),1:nx(n+1))=0;
         end   
      end
 
 
end