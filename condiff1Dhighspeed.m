%Central differencing scheme,1-D convection-diffusion steady state, no source but only boundary
%condition
%Author : Taewook Kim
%LINKEDIN : www.linkedin.com/in/kimtw
%GITHUB : github.com/Kimtaewookcode
%Email : kimtaewook87@gmail.com
%Ref:"An introduction to Computational Fluid Dynamics:HKversteeg and
%WMalalasekera example 5.1
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
        ae=D-F/2;
        sp=-(2*D+F);
        ap=aw+ae-sp;
        su=-sp*phi0;
        array(i,i)=ap;
        array(i,i+1)=-ae;
        bound(i,1)=su;
    elseif i<ngrid
        aw=D+F/2;
        ae=D-F/2;
        sp=0;
        ap=aw+ae-sp;
        su=0;%no source
        array(i,i-1)=-aw;
        array(i,i)=ap;
        array(i,i+1)=-ae;
        bound(i,1)=su;

    else
        aw=D+F/2;
        ae=0;
        sp=-(2*D-F);
        ap=aw+ae-sp;
        su=-sp*phil;
        
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
 ylim([0 1.5])
 