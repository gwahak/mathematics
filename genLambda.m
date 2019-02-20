function [Lambda] = genLambda(uHa)
%[u;H;a]의 상태를 받아서 그 상태에서의 Jacobian matrix의 고유값을 빠르게 구한다.
% Rusanov spllitting에 매우 유용하다.
u=uHa(1,:);
H=uHa(2,:);
a=uHa(3,:);
Lambda=[u-a;u;u+a];
end
