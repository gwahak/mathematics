function [uHa_half] = uHa_calculator_half(U)
%모든 점에서의 초기값 U=[rho; rho*u; E(=단위 부피당 전체 에너지)]를 입력으로 받아서
%반정수 점에서 [u(=속도,H(=specific enthalpy),a(=sound speed))를 Roe-Average 방법으로 구하여준다.
gamma=1.4;
[~,N] = size(U);
uHa_half=zeros(3,N+1);
uHa_half(:,1)=uHa_calculator(U(:,1));
uHa_half(:,N+1)=uHa_calculator(U(:,N));
for i=1:N-1
uHa_half(:,i+1)=Roe_Average(U(:,i:i+1));
end
end