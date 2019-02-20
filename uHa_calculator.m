function [uHa] = uHa_calculator(U)
%�� �������� �ʱⰪ U=[rho; rho*u; E(=���� ���Ǵ� ��ü ������)]�� �Է����� �޾Ƽ�
%�� ������ [u(=�ӵ�,H(=specific enthalpy),a(=sound speed))�� ���Ͽ��ش�.
%��� �������� u, H, a�� �ѹ��� ��ȯ�����ִ� ���� uHa_calculator2.m ��⿡ �ִ�.
gamma=1.4;
uHa(1,1)=U(2,1)/U(1,1);
uHa(2,1)=(U(3,1)+(gamma-1)*(U(3,1)-0.5*U(2,1)^2/U(1,1)))/U(1,1);
uHa(3,1)=sqrt((gamma-1)*(uHa(2,1)-uHa(1,1)^2/2));
end