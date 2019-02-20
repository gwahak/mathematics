function [dF] = WENO5core(U,N,h,wenooption)
%WENO5CORE 모듈은 flux의 derivative를 hybrid WENO 방법으로 계산해준다. 
%전체 point 개수, x축에 대한 node 간격 h, 유한차분법 scheme을 선택할 지에 대한 option을 받아서 
%weno-5, 6th order central scheme, hybrid 방법에 대해 정수 index를 가진 점에서의 flux를 각각 계산한뒤
%option에서 정한 방법으로 구한 flux에 대해 derivative를 계산한다.
F=build_Flux(U);
uHa_half=zeros(3,N+1);%반정수에서의 u, H, a를 구하기 위해 배열을 생성한다.
uHa_half(:,1)=uHa_calculator(U(:,1));
uHa_half(:,N+1)=uHa_calculator(U(:,N));%경계 조건에 맞게 index N+0.5 에서의 u,H,a를 계산한다. 
for i=1:N-1%Roe_Average 방법에 따라 index i+0.5 에서의 u,H,a를 계산한다.
uHa_half(:,i+1)=Roe_Average(U(:,i:i+1));
end
%반정수 index 에서의 Jacobian, eigenvalues, eigenvector, left eigenvector를 계산한다.
for i=1:N+1
Jacobian_half(:,:,i)=Jacobian(uHa_half(:,i));
[R_half(:,:,i),~, invR_half(:,:,i), eigs_half(:,i)]=decompose(uHa_half(:,i));
end
%Rusanov splliting을 수행하는데 사용될 lambdamax를 계산하는 단계이다.
for i=1:N
uHa_int(:,i)=uHa_calculator(U(:,i));
end
eigs_int=genLambda(uHa_int);
f_hybrid=zeros(3,N+1);
f_weno5=zeros(3,N+1);
f_central=zeros(3,N+1);
max_eig=max(abs(eigs_int(:,:)),[],2);
%WENO5modulehybrid 모듈에 값을 먹여준다. Boundary condition에 맞게 좌우측 각 3개의 값은 반복문에 넣지 않고 따로 계산한다. 
[f_hybrid(:,1),f_weno5(:,1),f_central(:,1)]=WENO5modulehybrid(R_half(:,:,1),F(:,[1 1 1 1:3]),U(:,[1 1 1 1:3]),max_eig,invR_half(:,:,1));
[f_hybrid(:,2),f_weno5(:,2),f_central(:,2)]=WENO5modulehybrid(R_half(:,:,2),F(:,[1 1 1:4]),U(:,[1 1 1:4]),max_eig,invR_half(:,:,2));
[f_hybrid(:,3),f_weno5(:,3),f_central(:,3)]=WENO5modulehybrid(R_half(:,:,3),F(:,[1 1:5]),U(:,[1 1:5]),max_eig,invR_half(:,:,3));

for i=1:N-5
[f_hybrid(:,i+3),f_weno5(:,i+3),f_central(:,i+3)]=WENO5modulehybrid(R_half(:,:,i+3),F(:,i:i+5),U(:,i:i+5),max_eig,invR_half(:,:,i+3));
end

[f_hybrid(:,N-1),f_weno5(:,N-1),f_central(:,N-1)]=WENO5modulehybrid(R_half(:,:,N-1),F(:,[N-4:N N]),U(:,[N-4:N N]),max_eig,invR_half(:,:,N-1));
[f_hybrid(:,N),f_weno5(:,N),f_central(:,N)]=WENO5modulehybrid(R_half(:,:,N),F(:,[N-3:N N N]),U(:,[N-3:N N N]),max_eig,invR_half(:,:,N));
[f_hybrid(:,N+1),f_weno5(:,N+1),f_central(:,N+1)]=WENO5modulehybrid(R_half(:,:,N+1),F(:,[N-2:N N N N]),U(:,[N-2:N N N N]),max_eig,invR_half(:,:,N+1));

switch wenooption
    case 'hybrid' % hybrid scheme 
        f_half=f_hybrid;
    case 'weno5' % weno5 scheme
        f_half=f_weno5;
    case 'central' % central scheme
        f_half=f_central;
    otherwise
        f_half=f_hybrid;
end
%projection된 flux를 복원한다.
for i=1:N+1
F_half(:,i)=R_half(:,:,i)*f_half(:,i);
end
dF=(F_half(:,2:N+1)-F_half(:,1:N))/h;
end