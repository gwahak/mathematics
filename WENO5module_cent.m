function [f_central] = WENO5module_cent(R_half,Flux,Que,lamMax,invR_half);
%계산 시간을 비교하기 위해 WENO5modulehybrid.m에서 central scheme만 수행하도록 군더더기를 삭제해 코드를
%수정한 모듈이다. 나머지는 WENO5modulehybrid.m과 같다.
qk=invR_half*Que;
fk=invR_half*Flux;
f_central=WENO5modulecentral(fk);
end

