function [Lambda] = genLambda(uHa)
%[u;H;a]�� ���¸� �޾Ƽ� �� ���¿����� Jacobian matrix�� �������� ������ ���Ѵ�.
% Rusanov spllitting�� �ſ� �����ϴ�.
u=uHa(1,:);
H=uHa(2,:);
a=uHa(3,:);
Lambda=[u-a;u;u+a];
end
