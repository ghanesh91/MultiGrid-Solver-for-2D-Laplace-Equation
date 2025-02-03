%%**************MULTIGIRD PROJECT******************
%%*****SUBMITTED BY GHANESH NARASIMHAN*************
%%*******NUMERICAL METHODS (FALL 2018)*************

clear
%% Multigrid parameters
Nlvlmax=8;            %Number of multigrid levels
Ncycle=500;           %Number of MG cycles
%tot_iter=2500;       %total iterations=Ncycle*maxiter*Nlvlmax;
maxiter=5;%tot_iter/(Nlvlmax*Ncycle);%Maximum iterations in iterative solver

omg=1.00;            %SOR parameter
Vcycle=1;            %Flag for V-Cycle
Wcycle=0;            %Flag for W-Cycle

        
%Type of smoother
itype=3;
%
%itype=1 Point Gauss-Seidel red black + SOR + MG
%itype=2 Point Gauss-Seidel normal    + SOR + MG
%itype=3 Line  Gauss-Seidel (ADI)     + SOR + MG
%itype=4 Line  Gauss-Seidel (ADI)     + SOR + Non-MG
%itype=5 Point Gauss-Seidel           + SOR + Non-MG

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
%% Solving Lu=R using Non-MG
if (itype>=4)
  n=1;tol=1e-5;counternMG=0;
  [uout{n},counternMG,tnonMG,errornMG]=NonMG(uin{n},RHS{n},invA{n},invB{n}...
                                             ,a(n),b(n),c(n),Nx,Ny,itype,omg,tol);
  semilogy(errornMG)
end

%% Solving Lu=R using MG (V-cycle)
if (Vcycle==1 && itype<4)
   tsV=cputime;
   for Ncyl=1:Ncycle
       Ncyl
      [uout]=fine_to_coarse(uin,RHS,maxiter,invA,invB,A,B,a,b,c,nx,ny,itype,omg,Nlvlmax,1,uout,eps,epsnew);
      [uout]=coarse_to_fine(Nlvlmax,uout,uoutnew,2);
      uin{1}=uout{1};
      %residue for testing convergence 
      error(Ncyl)=norm(residual(uout{1},RHS{1},a(1),b(1),c(1),nx(1),ny(1)));
%       if (error(Ncyl)<1e-5)
%           Ncyl
%           teV=cputime-tsV;
%           break
%       end
   end
end

%% Solving Lu=R using MG (W-cycle)
if (Wcycle==1 && itype<4 && Nlvlmax>3)
   tsW=cputime;
   for Ncyl=1:Ncycle
       Ncyl
      [uout]=fine_to_coarse(uin,RHS,maxiter,invA,invB,A,B,a,b,c,nx,ny,itype,...
                            omg,Nlvlmax,1,uout,eps,epsnew);
      [uout]=coarse_to_fine(Nlvlmax,uout,uoutnew,Nlvlmax);
      uin{Nlvlmax-1}=uout{Nlvlmax-1};
      
      [uout]=fine_to_coarse(uin,RHS,maxiter,invA,invB,A,B,a,b,c,nx,ny,itype,...
                            omg,Nlvlmax,Nlvlmax-1,uout,eps,epsnew);
      [uout]=coarse_to_fine(Nlvlmax,uout,uoutnew,Nlvlmax-1);
      uin{Nlvlmax-2}=uout{Nlvlmax-2};
      
      
      
      [uout]=fine_to_coarse(uin,RHS,maxiter,invA,invB,A,B,a,b,c,nx,ny,itype,omg,...
                            Nlvlmax,Nlvlmax-2,uout,eps,epsnew);
      [uout]=coarse_to_fine(Nlvlmax,uout,uoutnew,2);
      uin{1}=uout{1};    
      %residue for testing convergence 
      error(Ncyl)=norm(residual(uout{1},RHS{1},a(1),b(1),c(1),nx(1),ny(1)));
   end
   teW=cputime-tsW;
end


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
ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
title('$u$ (Line GS)','interpreter','latex','fontsize',16)
set(gcf,'Color','w')
set(gca,'fontsize',16,'fontname','times')
colorbar
colormap(flipud(gray))

figure(2)
semilogy(1:5:500,error(1:5:end));hold on
xlabel('Number of cycles','interpreter','latex','fontsize',16)
ylabel('$\epsilon =|\!|\nabla^2u-R|\!|$','interpreter','latex','fontsize',16)
title('Line \ GS','interpreter','latex','fontsize',16)
set(gcf,'Color','w')
set(gca,'fontsize',16,'fontname','times')

