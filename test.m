clear all;

B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
umin=[1;1;1;1]*(-20)*pi/180;
umax=[1;1;1;1]*20*pi/180;
[m,k] = size(B);
ep=0.1;
yd=[0.1;0.00000;0];
ye=[0.000;0.000;0];
% u_inv = B*Constrain(pinv(B)*(yd+ye),umin,umax)
% u1=SB_LPCA(yd,[0;0;0],B,ep*ones(4,1),Constrain(u_inv,umin,umax),umin,umax,5e2)
u1 = SB_LPCA(yd,ye,B,ep*ones(4,1),Constrain(pinv(B)*yd,umin,umax),umin,umax,5e2)
% u2 = B*SB_LPCA(yd,ye,B,ep*ones(4,1),zeros(4,1),umin,umax,5e2)%Constrain(pinv(B)*(yd+ye),umin,umax)
% u3 = B*DP_LPCA(yd+ye,B,umin,umax,5e2)
% u4= B* DPscaled_LPCA((yd+ye),B,umin,umax,5e2)



% u1=DB_LPCA(yd,B,ones(3,1),zeros(4,1),ones(4,1),ones(3,1)*0.5,umin,umax,5e2)
% u2=DBinf_LPCA(yd,B,ones(3,1),zeros(4,1),ones(4,1),ones(3,1)*0.5,umin,umax,5e2)
% u3 = MO_LPCA(yd,B,zeros(4,1),0.1, ones(3,1)*20,umin,umax,5e2)
%===================================∑˘÷µ≤‚ ‘==================================================
N=50;
x=zeros(k,(N+1)^2);
u=zeros(k,1);
[X,Y,Z] = sphere(N);
% t=0:0.01:1;
% X=0.3*sin(pi*t);
% Y=0.3*cos(pi*t);
% Z=0*sin(pi*t);
x1=zeros(k,(N+1)^2);
u1=zeros(k,1);
x2=zeros(k,(N+1)^2);
u2=zeros(k,1);

for i=1:(N+1)^2%%length(X)
v=0.5*[X(i);Y(i);Z(i)];% –Èƒ‚÷∏¡Ó
%=====================================
% u=pinv(B)*v;
% x(:,i)=Constrain(u,umin,umax);

% u1=DP_LPCA(v,B,umin,umax,100);
% x1(:,i)=Constrain(u,umin,umax);

% u1=DPscaled_LPCA(v,B,umin,umax,100);
% x1(:,i)=Constrain(u1,umin,umax);

u2=SB_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),umin,umax,5e2);
x2(:,i)=Constrain(u2,umin,umax);
end
U=B*x;
U1=B*x1;
U2=B*x2;

figure(1),
% plot3(U(1,:),U(2,:),U(3,:),'b*');
% hold on;
% plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
% hold on;
% plot3(ye(1,1),ye(2,1),ye(3,1),'b*');
% hold on;
% plot3(yd(1,1),yd(2,1),yd(3,1),'g*');
% hold on;
plot3(U2(1,:),U2(2,:),U2(3,:),'k*');


