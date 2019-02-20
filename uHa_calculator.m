function [uHa] = uHa_calculator(U)
%한 점에서의 초기값 U=[rho; rho*u; E(=단위 부피당 전체 에너지)]를 입력으로 받아서
%그 점에서 [u(=속도,H(=specific enthalpy),a(=sound speed))를 구하여준다.
%모든 점에서의 u, H, a를 한번에 변환시켜주는 것은 uHa_calculator2.m 모듈에 있다.
gamma=1.4;
uHa(1,1)=U(2,1)/U(1,1);
uHa(2,1)=(U(3,1)+(gamma-1)*(U(3,1)-0.5*U(2,1)^2/U(1,1)))/U(1,1);
uHa(3,1)=sqrt((gamma-1)*(uHa(2,1)-uHa(1,1)^2/2));
end