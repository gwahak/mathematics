function [uHa] = uHa_calculators(U)
%모든 점에서의 초기값 U=[rho; rho*u; E(=단위 부피당 전체 에너지)]를 입력으로 받아서
%모든 점에서 [u(=속도,H(=specific enthalpy),a(=sound speed))를 구하여준다.
%한 점에서의 u, H, a만 변환시켜주는 것은 uHa_calculator.m 모듈에 있다.
gamma=1.4;
uHa(1,:)=U(2,:)./U(1,:);
uHa(2,:)=(U(3,:)+(gamma-1)*(U(3,:)-0.5*U(2,:).^2./U(1,:)))./U(1,:);
uHa(3,:)=sqrt((gamma-1)*(uHa(2,:)-uHa(1,:).^2/2));
end