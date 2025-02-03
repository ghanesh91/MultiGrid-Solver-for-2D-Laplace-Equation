function [uout]=coarse_to_fine(Nlvlmax,uout,uoutnew,Nlvlmin)
      %Coarse to fine
      for n=Nlvlmax:-1:Nlvlmin
         fprintf('nb=%d \n',n-1)
         uoutnew{n-1}=prolongation(uout{n});
         uout{n-1}=uout{n-1}+uoutnew{n-1};    
      end
      

end