function [uHa] = Roe_Average( U)
%두 점에서의 초기값 U=[rho; rho*u; E(=단위 부피당 전체 에너지)]를 입력으로 받아서
%정 가운데의 반정수 점에서 [u(=속도,H(=specific enthalpy),a(=sound speed))를 Roe-Average 방법으로 구하여준다.
%모든 정수 점에서의 u, H, a에 대해, 반정수점의 u, H, a를 한번에 구해시켜주는 것은 uHa_calculator_half.m 모듈에 있다.
gamma=1.4;
rhoLsq=sqrt(U(1,1));
rhoRsq=sqrt(U(1,2));
HL=(U(3,1)+(gamma-1)*(U(3,1)-0.5*U(2,1)^2/U(1,1)))/U(1,1);
HR=(U(3,2)+(gamma-1)*(U(3,2)-0.5*U(2,2)^2/U(1,2)))/U(1,2);
uHa(1,1)=(U(2,1)/rhoLsq+U(2,2)/rhoRsq)/(rhoLsq+rhoRsq);
uHa(2,1)=(HL*rhoLsq+HR*rhoRsq)/(rhoLsq+rhoRsq);
uHa(3,1)=sqrt((gamma-1)*(uHa(2,1)-uHa(1,1)^2/2));
end