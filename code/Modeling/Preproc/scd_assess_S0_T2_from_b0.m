function [S0, T2, D] = scd_assess_S0_T2_from_b0(scheme, data, bmin, bmax)

% Estimate S0 and T2
% when bvalue = 0 we have : S0_experimental = S0*exp(-T2/TE)   (1)
%
% least square problem y = A*X :
%       - observations y : S0_experimental
%       - matrix for the equation (1) after linearisation ( "ln( (1) )" )
%       - vector of parameters X : for T2 and S0


% extract b0 values
indexb0 = scd_scheme2bvecsbvals(scheme)>=bmin & scd_scheme2bvecsbvals(scheme)<=bmax;
schemeb0 = scheme(indexb0, :);
datab0 = data(indexb0);

% extract TE and observations y
y = log(datab0(:));
TE = schemeb0(:, 7);
N = length(TE);

if length(unique(TE))>1
    
    A = [ones(N,1), -TE -scd_scheme2bvecsbvals(schemeb0)];
    
    % least square problem resolution
    %xMS = (A'*A)^-1*A'*y;
    xMS = lsqlin(A,double(y),[],[],[],[],[0 1/200 1/3000],[inf inf inf],[],optimoptions('lsqlin','Display','off'));
    
    
    S0 = exp(xMS(1));
    T2 = 1/xMS(2);
    D = 1/xMS(3);
else
    error('You need b=0 acquired at different echo times')
end

% 
% figure; hold on
% plot(datab0,'r*')
% plot(S0.*exp(-TE./T2),'b-*')
% legend('datab0 experimental', 'least square')

end