function f_weno5 = WENO5module_weno5(R_half,Flux,Que,lamMax,invR_half);
%계산 시간을 비교하기 위해 WENO5modulehybrid.m에서 WENO-5(WENO-JS)만 수행하도록 군더더기를 삭제해 코드를
%수정한 모듈이다. 나머지는 WENO5modulehybrid.m과 같다.
qk=invR_half*Que;
fk=invR_half*Flux;
f_neg=0.5*(fk-lamMax.*qk);
f_pos=0.5*(fk+lamMax.*qk);
f_weno5=WENO5module2(f_neg,-1)+WENO5module2(f_pos,+1);
end
