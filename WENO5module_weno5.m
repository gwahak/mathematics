function f_weno5 = WENO5module_weno5(R_half,Flux,Que,lamMax,invR_half);
%��� �ð��� ���ϱ� ���� WENO5modulehybrid.m���� WENO-5(WENO-JS)�� �����ϵ��� �������⸦ ������ �ڵ带
%������ ����̴�. �������� WENO5modulehybrid.m�� ����.
qk=invR_half*Que;
fk=invR_half*Flux;
f_neg=0.5*(fk-lamMax.*qk);
f_pos=0.5*(fk+lamMax.*qk);
f_weno5=WENO5module2(f_neg,-1)+WENO5module2(f_pos,+1);
end
