function [R, Lambda, invR, eigs] = decompose(uHa)
%���������� [u;H;a]�� ���¸� �޾Ƽ� �� ���¿����� Jacobian matrix�� ������ �����ϸ� ������ 
%right eigenvector matrix R, �������� ���� �밢���, 
%left eigenvector matrix�� inv(R), �������� ���� ���� eigs�� ������ش�. 
%  Jacobian matrix���� �������� ������ ���ϰ� ������ genLambda.m�� �̿��Ѵ�.
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
