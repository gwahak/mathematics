function [A] = Jacobian(uHa)
%Jacobian matrix를 구해주는 모듈이다.
gamma=1.4;
u=uHa(1,1);
H=uHa(2,1);
a=uHa(3,1);
A=zeros(3,3);
A(1,2)=1;
A(2,:)=[0.5*(gamma-3)*u^2 (3-gamma)*u gamma-1];
A(3,:)=[u*(0.5*(gamma-1)*u^2-H) H-(gamma-1)*u^2 gamma*u]; 
end