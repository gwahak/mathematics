function [R, Lambda, invR, eigs] = decompose(uHa)
%한점에서의 [u;H;a]의 상태를 받아서 그 상태에서의 Jacobian matrix를 고윳값 분해하면 나오는 
%right eigenvector matrix R, 고유값을 담은 대각행렬, 
%left eigenvector matrix인 inv(R), 고유값을 담은 벡터 eigs를 계산해준다. 
%  Jacobian matrix들의 고유값만 빠르게 구하고 싶으면 genLambda.m을 이용한다.
gamma=1.4;
u=uHa(1,1);
H=uHa(2,1);
a=uHa(3,1);
R=zeros(3,3);
R(:,1)=[1;u-a;H-u*a];
R(:,2)=[1;u;0.5*u^2];
R(:,3)=[1;u+a;H+u*a];

Lambda=[u-a 0 0; 0 u 0; 0 0 u+a];
invR=zeros(3,3);
invR(1,:)=[((gamma-1)/4 * (u/a)^2 +(u/a)/2) (-(gamma-1)/a *u/(a^2)-0.5/a) (gamma-1)/(2*a^2)];
invR(2,:)=[(1-(gamma-1)/2 * (u/a)^2) ((gamma-1)*u/(a^2)) ((1-gamma)/(a^2))];
invR(3,:)=[((gamma-1)/4 * (u/a)^2 -(u/a)/2) (-(gamma-1)/a *u/(a^2)+0.5/a) (gamma-1)/(2*a^2)];
eigs=[u-a;u;u+a];
end
