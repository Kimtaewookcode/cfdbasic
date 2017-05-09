%SIMPLE algorithm,2-D nozzel steady state, 
%Author : Taewook Kim
%LINKEDIN : www.linkedin.com/in/kimtw
%GITHUB : github.com/Kimtaewookcode
%Email : kimtaewook87@gmail.com
%Ref:"An introduction to Computational Fluid Dynamics:HKversteeg and
%WMalalasekera example 6.2
%%%%%%%%
clear all
clc 
clf

rho=1;%density
L=2;%length of geometry
nnode=5;%numberof node
nnodev=nnode-1;%number of cell
dx=L/nnodev;%delta x, uniform
AA=0.5;%Area at A
AE=0.1;%Area at E
iteration=100;%numberof maximum iteration
tor=10e-05;%maximum sum of residual 
und=0.8;%underrelaxation factor
xdarrp=zeros(1,nnode);%matrix of distance for pressure points
for i=1:nnode
xdarrp(1,i)=-dx+dx*i;%matrix of x distances for pressure points
end
Axp=zeros(1,nnode);%Areas of xdarrp
Axp=0.5-(AA-AE)/L*xdarrp;

xdarrv=zeros(1,nnodev);%matrix of x distances for velocity points
for i=1:nnodev
xdarrv(1,i)=xdarrp(1,i)+dx/2;
end
Axv=0.5-(AA-AE)/L*xdarrv;%matrix of areas of velocity points

p0=10;%pressure at the inlet[pa]
pE=0;%pressure at the outlet
mdot=1;%[kg/s]
u_ini=mdot/rho./Axv;
%%%%psuedo initial velocity
u_ini_p=mdot/rho./Axp;
p_ini=p0-p0/L*xdarrp;%%%%p guess-linear guess
d=zeros(1,nnodev);%parameter that is used for pressure corrector
%%%%%
array=zeros(nnodev);
bound=zeros(nnodev,1);
residuals=zeros(1,nnodev);
%%%%%
for j=1:iteration

for i=1:nnodev
    if i==1
        Fw=rho*u_ini_p(1,i)*Axp(1,i);
        Fe=rho*(u_ini(1,i)+u_ini(1,i+1))/2*Axp(1,i+1);
        aW=0;
        aE=0;
        aP=Fe+Fw*0.5*(Axv(1,i)/Axp(1,i))^2;
        Su=(p0-p_ini(1,i+1))*Axv(1,i)+Fw*Axv(1,i)/Axp(1,i)*u_ini(1,i);%u_ini is the velocity at previous iteration
        d(1,i)=Axv(1,i)/aP;

        array(i,i)=aP;
        array(i,i+1)=-aE;
        bound(i,1)=Su;
        residuals(1,i)=abs(aP*u_ini(1,i)-Su);

    elseif i<nnodev
        Fw=rho*(u_ini(1,i-1)+u_ini(1,i))/2*Axp(1,i);
        Fe=rho*(u_ini(1,i)+u_ini(1,i+1))/2*Axp(1,i+1);
        aW=Fw;
        aE=0;
        aP=aW+aE+(Fe-Fw);
        Su=(p_ini(1,i)-p_ini(1,i+1))*Axv(1,i);
        d(1,i)=Axv(1,i)/aP;
        
        array(i,i-1)=-aW;
        array(i,i)=aP;
        array(i,i+1)=-aE;
        bound(i,1)=Su;
        
        residuals(1,i)=abs(aP*u_ini(1,i)-aW*u_ini(1,i-1)-Su);

    else
        Fw=rho*(u_ini(1,i-1)+u_ini(1,i))/2*Axp(1,i);
        Fe=mdot;
        aW=Fw;
        aE=0;
        aP=aW+aE+(Fe-Fw);
        Su=(p_ini(1,i)-p_ini(1,i+1))*Axv(1,i);
        d(1,i)=Axv(1,i)/aP;

        array(i,i-1)=-aW;
        array(i,i)=aP;
        bound(i,1)=Su;
        residuals(1,i)=abs(aP*u_ini(1,i)-aW*u_ini(1,i-1)-Su);
    end
    
end
    u=array\bound;%%%%new velocity%%%%%
    
    %%%%pressure corrector
    arrayp=zeros(nnode-2,nnode);
    boundp=zeros(nnode-2,1);
    for i=2:nnode-1
            aW=rho*d(1,i-1)*Axv(1,i-1);
            aE=rho*d(1,i)*Axv(1,i);
            Fw=rho*u(i-1)*Axv(1,i-1);
            Fe=rho*u(i)*Axv(1,i);
            aP=aW+aE;
            b=Fw-Fe;
            arrayp(i-1,i-1)=-aW;
            arrayp(i-1,i)=aP;
            arrayp(i-1,i+1)=-aE;
            boundp(i-1,1)=b;
            
    end
arrayp(:,[1,nnode])=0;%pressure correction at the first node and last node is 0
pcorr=arrayp\boundp;%%correction pressure
p=p_ini+pcorr';%%corrected pressure
for i=1:nnodev
    u(i,1)=u(i,1)+d(1,i)*(pcorr(i,1)-pcorr(i+1,1));%corrected velocity
end
p(1,1)=p0-0.5*rho*(u(1,1)*Axv(1,1)/Axp(1,1))^2;%corrcted pressure at the first node
%underrelaxation
u_ini=u'*und+u_ini*(1-und); 
p_ini=p*und+p_ini*(1-und);
%%%%%%%
mdot=rho*u_ini(1,1)*Axv(1,1);
u_ini_p=mdot/rho./Axp;%u_ini_p updated
resi=sum(residuals(:));%sume of residuals
if resi<tor
    disp(j)
    break%simulation stops when resi is smaller than tor, and shows iterations
end
end
disp(mdot)
xber= [0 : 0.01: 2];
Aber=0.5-(0.5-0.1)/L*xber;
mreal=((2*(p0-0)*(rho*AE)^2)/rho)^0.5;
preal=p0-(0.5*rho*mreal^2)./((rho*Aber).^2);
vreal=mreal./Aber/rho;
%plot(xdarrv,mdot,'x-',xdarrv,mreal,'.-')
plot(xdarrp,p_ini,'x-',xber,preal,'.-')
xlim([0 2])
figure;plot(xdarrv,u_ini,'x-',xber,vreal,'.-')
xlim([0 2])