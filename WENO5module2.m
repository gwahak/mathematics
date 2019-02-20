function res = WENO5module2(g,dir)
% we can obtain three ENO fluxes with the 3rd order accuracy as:
if dir == -1
 g= fliplr(g);
end
gm2 = g(:,1);
gm1  = g(:,2);
g0 = g(:,3);
gp1  = g(:,4);
gp2 = g(:,5);

g03 = (2*gm2 - 7*gm1 + 11*g0)/6;
g13 = ( -gm1 + 5*g0  + 2*gp1)/6;
g23 = (2*g0  + 5*gp1 - gp2 )/6;
% Efficient Implementation of Weighted ENO Schemes
c03=0.1; c13= 0.6; c23= 0.3; epweno=1e-7;

% Smooth Indicators (Beta factors)
beta0 = 0.25*(gm2-4*gm1+3*g0).^2+ 13/12*(gm2-2*gm1+g0).^2 ;  
beta1 = 0.25*(gm1-gp1).^2+ 13/12*(gm1 -2*g0 +gp1).^2 ;
beta2 = 0.25*(3*g0-4*gp1+gp2).^2+ 13/12*(g0-2*gp1+gp2).^2 ;

W00 = c03 ./ (epweno + beta0).^2; 
W01 = c13 ./ (epweno + beta1).^2; 
W02 = c23 ./ (epweno + beta2).^2;
Wsum=W00+W01+W02;
W0=W00./Wsum; W1=W01./Wsum;; W2=W02./Wsum;

res= g03.*W0 + g13.*W1 + g23.*W2;
end

