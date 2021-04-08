clear all;
close all;
B=[-0.5   0       0.5   0;
    0  -0.5    0       0.5;
    0.25   0.25   0.25   0.25];
uMin=[1;1;1;1]*(-20)*pi/180;
uMax=[1;1;1;1]*20*pi/180;
% B =[
% 
%     0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
%     1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
%    -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
% uMin =[
% 
%    -0.9599;
%    -0.9599;
%    -0.5236;
%    -0.5236;
%    -0.5236;
%    -0.5236;
%    -0.5236];
% uMax =[
% 
%     0.4363;
%     0.4363;
%     0.5236;
%     0.5236;
%     0.5236;
%     0.5236;
%     0.5236];
% yd =[
% 
%    -0.0645;
%    -0.1007;
%    -0.0000];
%    
%    
% 
% uMin =[
% 
%    -0.0318;
%    -0.0349;
%    -0.0262;
%    -0.0262;
%    -0.0262;
%    -0.0262;
%    -0.0175];
%    
%    uMax =[
% 
%     0.0349;
%     0.0349;
%          0;
%          0;
%     0.0262;
%     0.0262;
%     0.0175];
% uMin =[
%     0.4673;
%     0.4840;
%     0.6981;
%     0.6981];
% uMax =[-0.2309;
%    -0.2141;
%          0;
%          0];
[m,k] = size(B);
ep=0.1;
yd=[0;0.2;0]; 
% ye=[0;0;1.2];
u_inv = Constrain(pinv(B)*(yd),uMin,uMax)
% u1=SBprio_LPCA(yd,ye,B,ep*ones(4,1),Constrain(u_inv,uMin,uMax),uMin,uMax,5e2)
% u1 = SBprio_LPCA(yd,ye,B,ep*ones(4,1),Constrain(pinv(B)*yd,uMin,uMax),uMin,uMax,5e2)
% u1 = SBprio_LPCA(yd,ye,B,ep*ones(4,1),zeros(4,1),uMin,uMax,5e2)
% u2 = B*SBprio_LPCA(yd,[0;0;0],B,ep*ones(7,1),zeros(7,1),uMin,uMax,5e2)%Constrain(pinv(B)*(yd+ye),uMin,uMax)
u3 = DP_LPCA(yd,B,uMin,uMax,5e2)
u4=  DPscaled_LPCA((yd),B,uMin,uMax,5e2)
u5=SB_LPCA(yd,B,ep*ones(k,1),zeros(k,1),uMin,uMax,5e2)


% B*DBcaLP1f_sol(yd,B,ones(3,1),ones(3,1)*20,zeros(7,1),uMin,uMax,m,k,5e2)
% u1=B*DB_LPCA(yd,B,ones(3,1),zeros(7,1),ep*ones(7,1),ones(3,1)*20,uMin,uMax,5e2)
% u2=DBinf_LPCA(yd,B,ones(3,1),zeros(4,1),ones(4,1),ones(3,1)*0.5,uMin,uMax,5e2)
% u3 = MO_LPCA(yd,B,zeros(4,1),0.1, ones(3,1)*20,uMin,uMax,5e2)
%===================================��ֵ����==================================================
% N=100;
% x=zeros(k,(N+1)^2);
% u=zeros(k,1);
% [X,Y,Z] = sphere(N);
% load('M_des');
% % t=0:0.01:1;
% % X=0.3*sin(pi*t);
% % Y=0.3*cos(pi*t);
% % Z=0*sin(pi*t);
% x1=zeros(k,(N+1)^2);
% u1=zeros(k,1);
% % x2=zeros(k,(N+1)^2);
% % u2=zeros(k,1);
% 
% for i=1:(N+1)^2%%length(X)%length(M_des(1:1000,1))%
% v=10*[X(i);Y(i);Z(i)];% ����ָ��M_des(i,:)'%
% %=====================================
% % u=pinv(B)*v;
% % x(:,i)=Constrain(u,uMin,uMax);
% 
% u1=DP_LPCA(v,B,uMin,uMax,100);
% x1(:,i)=Constrain(u,uMin,uMax);
% 
% % u1=DPscaled_LPCA(v,B,uMin,uMax,100);
% % x1(:,i)=Constrain(u1,uMin,uMax);
% 
% % u2=SBprio_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),uMin,uMax,5e2);
% % x2(:,i)=Constrain(u2,uMin,uMax);
% end
% U=B*x;
% U1=B*x1;
% % U2=B*x2;
% % V=(yd+ye);
% figure(1),
% % plot3(U(1,:),U(2,:),U(3,:),'b*');
% % hold on;
% plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
% % hold on;
% % plot3(ye(1,1),ye(2,1),ye(3,1),'b*');
% % hold on;
% % plot3(V(1,1),V(2,1),V(3,1),'g*');
% % hold on;
% % plot3(u2(1,:),u2(2,:),u2(3,:),'k*');
% % 
% % hold on;
% % plot3(u4(1,:),u4(2,:),u4(3,:),'g>');