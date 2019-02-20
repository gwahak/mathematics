function [F] = build_Flux( U )
%U로부터 Flux를 구해주는 모듈이다.
gamma=1.4;
rho=U(1,:);
rhou=U(2,:);
totE=U(3,:);
F(1,:)=rhou;
F(2,:)=(gamma-1).*(totE-0.5*rhou.^2./rho)+rhou.^2./rho;%(gamma-1)*totE+rhou.^2./rho;
F(3,:)=(totE+(gamma-1).*(totE-0.5*rhou.^2./rho)).*(rhou./rho);
end
