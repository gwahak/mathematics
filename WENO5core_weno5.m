function [dF] = WENO5core_weno5(U,N,h,wenooption)
%WENO5CORE_weno5 모듈은 flux의 derivative를 WENO-5(WENO-JS) 방법으로 계산해준다. 
%WENO5CORE 모듈에서 불필요한 다른 연산이 시행되지 않도록 일부 부분을 제거 및 수정 한 모듈이다.
F=build_Flux(U);
uHa_half=zeros(3,N+1);
uHa_half(:,1)=uHa_calculator(U(:,1));
uHa_half(:,N+1)=uHa_calculator(U(:,N));
for i=1:N-1
uHa_half(:,i+1)=Roe_Average(U(:,i:i+1));
end
for i=1:N+1
Jacobian_half(:,:,i)=Jacobian(uHa_half(:,i));
[R_half(:,:,i),~, invR_half(:,:,i), eigs_half(:,i)]=decompose(uHa_half(:,i));
end
for i=1:N
uHa_int(:,i)=uHa_calculator(U(:,i));
end
eigs_int=genLambda(uHa_int);

f_weno5=zeros(3,N+1);

max_eig=max(abs(eigs_int(:,:)),[],2);
f_weno5(:,1)=WENO5module_weno5(R_half(:,:,1),F(:,[1 1 1 1:3]),U(:,[1 1 1 1:3]),max_eig,invR_half(:,:,1));
f_weno5(:,2)=WENO5module_weno5(R_half(:,:,2),F(:,[1 1 1:4]),U(:,[1 1 1:4]),max_eig,invR_half(:,:,2));
f_weno5(:,3)=WENO5module_weno5(R_half(:,:,3),F(:,[1 1:5]),U(:,[1 1:5]),max_eig,invR_half(:,:,3));

for i=1:N-5
f_weno5(:,i+3)=WENO5module_weno5(R_half(:,:,i+3),F(:,i:i+5),U(:,i:i+5),max_eig,invR_half(:,:,i+3));
end

f_weno5(:,N-1)=WENO5module_weno5(R_half(:,:,N-1),F(:,[N-4:N N]),U(:,[N-4:N N]),max_eig,invR_half(:,:,N-1));
f_weno5(:,N)=WENO5module_weno5(R_half(:,:,N),F(:,[N-3:N N N]),U(:,[N-3:N N N]),max_eig,invR_half(:,:,N));
f_weno5(:,N+1)=WENO5module_weno5(R_half(:,:,N+1),F(:,[N-2:N N N N]),U(:,[N-2:N N N N]),max_eig,invR_half(:,:,N+1));

f_half=f_weno5;

for i=1:N+1
F_half(:,i)=R_half(:,:,i)*f_half(:,i);
end
dF=(F_half(:,2:N+1)-F_half(:,1:N))/h;
end