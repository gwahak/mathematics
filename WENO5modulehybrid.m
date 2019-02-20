function [f_half,f_weno5,f_central] = WENO5modulehybrid(R_half,Flux,Que,lamMax,invR_half);

qk=invR_half*Que;%characteristic space로 flux와 Q를 projection한다.
fk=invR_half*Flux;
f_neg=0.5*(fk-lamMax.*qk);
f_pos=0.5*(fk+lamMax.*qk);%Rusanov splitting
f_weno5=WENO5module2(f_neg,-1)+WENO5module2(f_pos,+1);
f_central=WENO5modulecentral(fk);%central scheme
epsilon2=1e-6;%ACM-switch 시작
r_hat_phalf=abs(abs(qk(:,4))-abs(qk(:,3)))./abs(abs(qk(:,4))+abs(qk(:,3))+epsilon2);
r_hat_mhalf=abs(abs(qk(:,3))-abs(qk(:,2)))./abs(abs(qk(:,3))+abs(qk(:,2))+epsilon2);
r=max([r_hat_phalf r_hat_mhalf],[],2);%ACM switch를 통해 weight를 구했다.
sigma=1-r;
f_half=sigma.*f_central+r.*f_weno5;%hybrid된 결과를 계산한다.
end