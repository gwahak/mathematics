function [uHa_half] = uHa_calculator_half(U)
%��� �������� �ʱⰪ U=[rho; rho*u; E(=���� ���Ǵ� ��ü ������)]�� �Է����� �޾Ƽ�
%������ ������ [u(=�ӵ�,H(=specific enthalpy),a(=sound speed))�� Roe-Average ������� ���Ͽ��ش�.
gamma=1.4;
[~,N] = size(U);
uHa_half=zeros(3,N+1);
uHa_half(:,1)=uHa_calculator(U(:,1));
uHa_half(:,N+1)=uHa_calculator(U(:,N));
for i=1:N-1
uHa_half(:,i+1)=Roe_Average(U(:,i:i+1));
end
end