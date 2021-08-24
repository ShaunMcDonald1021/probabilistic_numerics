function [k kmean Ikmean kwol kdl] = kq_kernel(kernel, l, d, distribution)
% KQ_KERNEL - select kernel for kernel quadrature
%   Given a kernel, its length-scale l, the 
%   dimension d and the integration distribution, 
%   [k kmean Ikmean] = kq_kernel(kernel, l, d, distribution)
%   generates the kernel function, kernel mean function
%   and the initial KQ worst-case error.
%
% INPUT
%   - kernel        the kernel, possible values:
%                     'gauss' Gaussian (the a.k.a. SE, RBF) kernel
%   - l             kernel length-scale, a positive real
%   - d             dimension
%   - distribution  integration distribution
%
% OUTPUT
%   - k             kernel function
%   - kmean         kernel mean function
%   - Ikmean        initial kernel quadrature WCE
%   - kwol          kernel withouth length-scale for length-scale fitting
%   - kdl           kernel derivative w.r.t. length-scale for length-scale fitting

% Toni Karvonen, 2017
  
  % GAUSSIAN KERNEL
  % possible distributions:
  %   'normal'    - standard normal distribution on R^d
  %   'uniform'   - normalizing uniform distribution on [-1,1]^d  
  if strcmp(kernel, 'gauss')
    k    = @(r) exp(-r.^2/(2*l^2));
    kwol = @(r) exp(-r.^2/2);
    kdl  = @(r,l) r.^2 .* exp(-r.^2/(2*l.^2))./l.^3;
    if strcmp(distribution, 'normal')
      kmean   = @(x) (l^2 / (1+l^2))^(d/2) * exp( -norm(x)^2 /(2*(1+l^2)) );
      Ikmean  = (l^2/(2+l^2))^(d/2);
    elseif strcmp(distribution, 'uniform')
      kmean   = @(x) (pi*l^2/8)^(d/2) * ...
                      prod( erf((x+1)/(sqrt(2)*l)) - erf((x-1)/(sqrt(2)*l)) );
      Ikmean  = (pi*l^2/8)^(d/2) * ...
                  ( sqrt(2*l^2/pi)*(exp(-2/l^2) - 1) + 2*erf(sqrt(2)/l) )^d;
    end
  end

end
