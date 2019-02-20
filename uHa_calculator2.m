function [uHa] = uHa_calculators(U)
%��� �������� �ʱⰪ U=[rho; rho*u; E(=���� ���Ǵ� ��ü ������)]�� �Է����� �޾Ƽ�
%��� ������ [u(=�ӵ�,H(=specific enthalpy),a(=sound speed))�� ���Ͽ��ش�.
%�� �������� u, H, a�� ��ȯ�����ִ� ���� uHa_calculator.m ��⿡ �ִ�.
gamma=1.4;
uHa(1,:)=U(2,:)./U(1,:);
uHa(2,:)=(U(3,:)+(gamma-1)*(U(3,:)-0.5*U(2,:).^2./U(1,:)))./U(1,:);
uHa(3,:)=sqrt((gamma-1)*(uHa(2,:)-uHa(1,:).^2/2));
end