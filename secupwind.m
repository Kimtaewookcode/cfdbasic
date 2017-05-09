%Second order upwind scheme,1-D convection-diffusion steady state, no source but only boundary
%condition
%Author : Taewook Kim
%LINKEDIN : www.linkedin.com/in/kimtw
%GITHUB : github.com/Kimtaewookcode
%Email : kimtaewook87@gmail.com
%Ref:"Lars Davidson: Numerical Methods for Turbulent Flow chapter5
%%%%%%%%
clear
clc
clf
L=1;%[m] length of geometry
%rho=1;%[kg/m3]
ksi=0.1;%[kg/ms]
u=2.5;%[m/s]%2.5;
ngrid=20;%numberof grid%20,5
dx=L/ngrid;%delta x

F=2.5;%coefficient
D=ksi/dx;%coefficient

phi0=1;%boundary condition at x=0
phil=0;%boundary condition at x=L

%%%%%%%%
array=zeros(ngrid);
bound=zeros(ngrid,1);
%%%%%%%%%
 
for i=1:ngrid
    if i==1
        aw=0;
        ae=D;
        ap=3/2*F+3*D;
        su=(3/2*F+2*D)*phi0;
        array(i,i)=ap;
        array(i,i+1)=-ae;
        bound(i,1)=su;
    elseif i==2
        aw=D+2*F;
        ae=D;
        ap=3/2*F+2*D;
        su=-F/2*phi0;%no source
        array(i,i-1)=-aw;
        array(i,i)=ap;
        array(i,i+1)=-ae;
        bound(i,1)=su;
        
    elseif i<ngrid
        aw=D+2*F;
        ae=D;
        ap=3/2*F+2*D;
        aww=-F/2;
        su=0;%no source
        array(i,i-2)=-aww;
        array(i,i-1)=-aw;
        array(i,i)=ap;
        array(i,i+1)=-ae;
        bound(i,1)=su;
        else
        aw=D+2*F;
        ae=0;
        ap=3/2*F+3*D;
        aww=-F/2;
        su=2*D*phil;
        array(i,i-2)=-aww;
        array(i,i-1)=-aw;
        array(i,i)=ap;
        bound(i,1)=su;
    end
end
 phi=array\bound;
 
 xarr=zeros(ngrid,1);
 for j=1:ngrid
     xarr(j,1)=dx/2+dx*(j-1);
 end
 
 %%%%%%%analytical
 dx1=dx/10;
 x1=[0:dx1:L];
 anaphi=1+(1-exp(25*x1))/(7.20*10^10);
 
 plot(xarr,phi,'x',x1,anaphi,'-')
 