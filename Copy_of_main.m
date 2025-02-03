clear
%% Multigrid parameters
Nlvlmax=4;      %Number of multigrid levels
maxiter=6;     %Maximum iterations in iterative solver
itype=3;        %Type of iterative solve
Ncycle=100;      %Number of MG cycles
omg=1;          %SOR parameter
tot_iter=Ncycle*maxiter*Nlvlmax;
%itype=1 Point Gauss-Seidel red black 
%itype=2 Point Gauss-Seidel normal

%% Grid definition for each level
Lx=2*pi;Ly=2*pi;
Nx=257;Ny=257; 
nx(1)=Nx;
ny(1)=Ny;
%Number of grid points at other levels
for i=2:Nlvlmax
    ny(i)=((ny(1)-1)/2^(i-1))+1;
    nx(i)=((nx(1)-1)/2^(i-1))+1;
end

%% Initialise arrays on each level
error(1:Ncycle)=0;
for i=1:Nlvlmax
    uin{i}(1:ny(i),1:nx(i))=0;
    uout{i}(1:ny(i),1:nx(i))=0;
    uoutnew{i}(1:ny(i),1:nx(i))=0;
    eps{i}(1:ny(i),1:nx(i))=0;
    epsnew{i}(1:ny(i),1:nx(i))=0;
    RHS{i}(1:ny(i),1:nx(i))=0;
end

%% Initial condition (random noise+ boundary condition)
for n=1:Nlvlmax
  dx=Lx/(nx(n)-1);dy=Ly/(ny(n)-1);
  x=0:dx:Lx;
  y=0:dy:Ly;
  if (n==1)
   %uin{n}= -1 + (1+1)*rand(ny(1),nx(1));
   u=load('init_cond.mat');
   uin{n}(1:ny(n),1:nx(n))=u.u(1:ny(n),1:nx(n));
   uin{n}(1:ny(n),1)=sin(4*y);
   uin{n}(1:ny(n),nx(n))=0;
   uin{n}(1,1:nx(n))=sin(4*x);
   uin{n}(ny(n),1:nx(n))=0;
  end   
end

%% Calculate co-efficients
a(1:Nlvlmax)=0;b(1:Nlvlmax)=0;c(1:Nlvlmax)=0;
for n=1:Nlvlmax
    A{n}=0;B{n}=0;invA{n}=0;invB{n}=0;
    [invA{n},invB{n},A{n},B{n},a(n),b(n),c(n)]=coeff(n,Nx,Ny,Lx,Ly);
end
%% Solving using Non-MG
if (itype==4 || itype==5)
  n=1;tol=1e-5;counternMG=0;
  [uout{n},counternMG,tnonMG,errornMG]=NonMG(uin{n},RHS{n},invA{n},invB{n},a(n),b(n),c(n),Nx,Ny,itype,omg,tol);
  semilogy(errornMG)
end

%% Solving using MG
if (itype~=4 || itype~=5)
   %% V-cycle
   ts=cputime;

   for Ncyl=1:Ncycle
       Ncyl
      %Fine to coarse 
      for n=1:Nlvlmax
         n;
         uout{n}=iterative_solve(uin{n},RHS{n},maxiter,invA{n},invB{n},A{n},B{n},a(n),b(n),c(n),nx(n),ny(n),itype,omg);
         if (n==1)
            uf{n}=uout{n}; 
         end
         eps{n}=residual(uout{n},RHS{n},a(n),b(n),c(n),nx(n),ny(n));
         if(n~=Nlvlmax)
            epsnew{n+1}=restriction(eps{n},nx(n),ny(n));
            RHS{n+1}=-epsnew{n+1};
            uin{n+1}(1:ny(n+1),1:nx(n+1))=0;
         end   
      end
 
      %Coarse to fine
      for n=Nlvlmax:-1:2
         n;
         uoutnew{n-1}=prolongation(uout{n});
         uout{n-1}=uout{n-1}+uoutnew{n-1};
    
      end
      if (Ncycle>1)
         uin{1}=uout{1};
      end
      %convergence analysis
      error(Ncyl)=norm(residual(uout{1},RHS{1},a(1),b(1),c(1),nx(1),ny(1)));
%       if (error(Ncyl)<1e-5)
%           Ncyl
%           break
%       end
   end
end

te=cputime-ts
%% Analysis
figure(3)
Lx=2*pi;Ly=2*pi;
dx=Lx/(nx(1)-1);dy=Ly/(ny(1)-1);
x=0:dx:2*pi;
y=0:dy:2*pi;
[X,Y]=meshgrid(x,y);
surf(x,y,uout{1},'linestyle','none')
set(gca, 'CameraPosition', [2*pi 2*pi 0.25]);
xlabel('$x$','interpreter','latex','fontsize',16)
ylabel('$y$','interpreter','latex','fontsize',16)
title('$u$','interpreter','latex','fontsize',16)
set(gcf,'Color','w')
set(gca,'fontsize',16,'fontname','times')
colorbar
colormap(flipud(gray))

figure(1)
semilogy(error);hold on
xlabel('Number of cycles','interpreter','latex','fontsize',16)
ylabel('$\epsilon (|\!|\nabla^2u-R|\!|)$','interpreter','latex','fontsize',16)
title('$u$','interpreter','latex','fontsize',16)
set(gcf,'Color','w')
set(gca,'fontsize',16,'fontname','times')

