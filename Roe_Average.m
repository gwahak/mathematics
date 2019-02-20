function [uHa] = Roe_Average( U)
%�� �������� �ʱⰪ U=[rho; rho*u; E(=���� ���Ǵ� ��ü ������)]�� �Է����� �޾Ƽ�
%�� ����� ������ ������ [u(=�ӵ�,H(=specific enthalpy),a(=sound speed))�� Roe-Average ������� ���Ͽ��ش�.
%��� ���� �������� u, H, a�� ����, ���������� u, H, a�� �ѹ��� ���ؽ����ִ� ���� uHa_calculator_half.m ��⿡ �ִ�.
gamma=1.4;
rhoLsq=sqrt(U(1,1));
rhoRsq=sqrt(U(1,2));
HL=(U(3,1)+(gamma-1)*(U(3,1)-0.5*U(2,1)^2/U(1,1)))/U(1,1);
HR=(U(3,2)+(gamma-1)*(U(3,2)-0.5*U(2,2)^2/U(1,2)))/U(1,2);
uHa(1,1)=(U(2,1)/rhoLsq+U(2,2)/rhoRsq)/(rhoLsq+rhoRsq);
uHa(2,1)=(HL*rhoLsq+HR*rhoRsq)/(rhoLsq+rhoRsq);
uHa(3,1)=sqrt((gamma-1)*(uHa(2,1)-uHa(1,1)^2/2));
end