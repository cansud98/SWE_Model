%RK4/LeapFrog+CIS schemes for SWE-1D
clear all, clc, close all
% parameters
xstart =          -5; %Length of left side
xend   =           5; %Length of right side
L      = xend-xstart; %Length of domain (m)
T      =         300; % duration of simulation (s)
dx     =        0.05; %grid spacing (m)
dt     =        0.01; %time step (s)
f      =      0.0001; %coriolis for 43N for exp
g      =           1; %gravitational force (m/s2)
H      =        0.01; %m
perturb=         0.2;  

Nt = floor(T/(dt)+1); %number of time steps
Nx = floor(L/(dx)+1); % number of grid points

scheme = input(['Which time-stepping scheme do you prefer? \n[1] Leap Frog Scheme \n[2] Runge-Kutta 4th order \nYour selection:']);
       if scheme == 1
           fprintf('Building the model with Leap Frog time-stepping scheme...');
       elseif scheme == 2
           fprintf('Building the model with Runge-Kutta 4th order time-stepping scheme...');
       else
           scheme = input(['Please select one of them only. \n[1] Leap Frog Scheme \n[2] Runge-Kutta 2nd order \n[3] Runge-Kutta 4th order \nYour selection:']);
       end

%Creating empty matrix
u=zeros(Nx,Nt);
v=zeros(Nx,Nt);
h=zeros(Nx,Nt);

%Initializing variables
for j = 1 : Nx
    x=xstart:dx:xend;
    t=0:dt:T;
    u(j, 1) = 0;
    v(j, 1) = 0;
    h(j, 1) = 1 + perturb*exp(-x(j)^2);
end
t(1)=0;
for n=2:Nt
   t(n)=(n-1)*dt;
   %n
   for j=1:Nx

        if j == 1
                xb = Nx;
                xc = j;
                xf = j+1;
        elseif j == Nx
                xb = j-1;
                xc = j;
                xf = 1;
        else 
                xb = j-1;
                xc = j;
                xf = j+1;
        end

       %Calculations
       hx(j,n-1) = (h(xf,n-1)-h(xb,n-1))/(2*dx);
       if j == 1 | Nx
           hx(j,n) = 0;
       end
       ux(j,n-1) = (u(xf,n-1)-u(xb,n-1))/(2*dx);

       %Summing up everything
       Fh(j,n-1)=-H*ux(j,n-1);
       Fu(j,n-1)=-g*hx(j,n-1)+f*v(j,n-1);
       Fv(j,n-1)=-f*u(j,n-1);

       if scheme == 1

       if n==2  | n/20==floor(n/20)
            %n
            u(j,n)=u(j,n-1)+dt*(Fu(j,n-1));
            v(j,n)=v(j,n-1)+dt*(Fv(j,n-1));
            h(j,n)=h(j,n-1)+dt*(Fh(j,n-1));
       else
            %n
            u(j,n)=u(j,n-2)+2*dt*(Fu(j,n-1));
            v(j,n)=v(j,n-2)+2*dt*(Fv(j,n-1));
            h(j,n)=h(j,n-2)+2*dt*(Fh(j,n-1));
       end

       elseif scheme == 2

       %RK-4th order step1,2,3,and4 for "u"
       ku1(j,n-1)=Fu(j,n-1);
       ku2(j,n-1)=Fu(j,n-1)+ku1(j,n-1)*dt/2;
       ku3(j,n-1)=Fu(j,n-1)+ku2(j,n-1)*dt/2;
       ku4(j,n-1)=Fu(j,n-1)+ku3(j,n-1)*dt;

       u(j,n)=u(j,n-1)+(ku1(j,n-1)+2*ku2(j,n-1)+2*ku3(j,n-1)+ku4(j,n-1))*dt/6;

       %RK-4th order step1,2,3,and4 for "h"
       kh1(j,n-1)=Fh(j,n-1);
       kh2(j,n-1)=Fh(j,n-1)+kh1(j,n-1)*dt/2;
       kh3(j,n-1)=Fh(j,n-1)+kh2(j,n-1)*dt/2;
       kh4(j,n-1)=Fh(j,n-1)+kh3(j,n-1)*dt;

       h(j,n)=h(j,n-1)+(kh1(j,n-1)+2*kh2(j,n-1)+2*kh3(j,n-1)+kh4(j,n-1))*dt/6;

       %RK-4th order step1,2,3,and4 for "v"
       kv1(j,n-1)=Fv(j,n-1);
       kv2(j,n-1)=Fv(j,n-1)+kv1(j,n-1)*dt/2;
       kv3(j,n-1)=Fv(j,n-1)+kv2(j,n-1)*dt/2;
       kv4(j,n-1)=Fv(j,n-1)+kv3(j,n-1)*dt;

       v(j,n)=v(j,n-1)+(kv1(j,n-1)+2*kv2(j,n-1)+2*kv3(j,n-1)+kv4(j,n-1))*dt/6;
       end
       end
end

%%
%% If you don't want to see it as a video delete following lines.
%%

%Movie
ifreq=100;
icount=1;
figure
for ii=1:ifreq:n
    hold off
    plot(x,h(:,ii),'b');ylim([0.8 1.4]);

                if scheme == 1
                    title('Shallow Water Equations w/LF');
                elseif scheme == 2
                    title('Shallow Water Equations w/RK4');
                end

    xlabel('Distance(m)')
    ylabel('Height(m)')
    Fmov(icount)=getframe;
    icount=icount+1;
end
