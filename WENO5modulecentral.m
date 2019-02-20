function [f] = WENO5modulecentral(g)
f=(g(:,1)-8*g(:,2)+37*g(:,3)+37*g(:,4)-8*g(:,5)+g(:,6))/60;
end