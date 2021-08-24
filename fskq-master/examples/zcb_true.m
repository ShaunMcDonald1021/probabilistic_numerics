function P = zcb_true(kappa, theta, sigma, r0, T, d)
% ZCB_TRUE - zero coupon bond ground truth
%   Ground truth for the price of zero coupon bonds.

% Toni Karvonen, 2017

  dt = T/d;
  
  g = 0;
  kappadt = 1-kappa*dt;
  
  for k = 1:d-1
    b = (1-kappadt^k)/(kappa*dt);
    g = g + b*kappa*theta*dt - (b*sigma*dt)^2/2;
  end
  b = (1-kappadt^d)/(1-kappadt);
  
  P = exp(-(g+b*r0)*T/d);
  
end
