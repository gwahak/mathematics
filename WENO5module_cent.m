function [f_central] = WENO5module_cent(R_half,Flux,Que,lamMax,invR_half);
%��� �ð��� ���ϱ� ���� WENO5modulehybrid.m���� central scheme�� �����ϵ��� �������⸦ ������ �ڵ带
%������ ����̴�. �������� WENO5modulehybrid.m�� ����.
qk=invR_half*Que;
fk=invR_half*Flux;
f_central=WENO5modulecentral(fk);
end

