function Y = zcb_integrand(Z, kappa, theta, sigma, r0, T)
% ZCB_INTEGRAND - integrand function for zero coupon bonds
%   Ground truth for the price of zero coupon bonds.
%   The input points are columns of Z that need to
%   have the standard normal distribution.

% Toni Karvonen, 2017
  
  [d N] = size(Z);
  dt = T/(d+1);
  Z  = sqrt(dt)*Z;
  r = r0*ones(1,N);
  rsum = r0*ones(1,N);
    
  for k = 1:d
    r = r + kappa*(theta - r)*dt + sigma*Z(k,:);
    rsum = rsum + r;
  end
  Y = exp(-dt*rsum);
  
end
