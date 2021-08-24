%%%%%%%%%
%%%
%%%  EXAMPLE 3:
%%%  Gauss-Hermite sparse grid kernel quadrature on a zero
%%%  coupon bond model. Integration in dimensions 9 to 299.
%%%
%%% Toni Karvonen, 2017
%%%
%%%%%%%%%

  %% Initializations
    addpath('../fskq')

    % Permanent parameters
    T = 5;
    kappa = 0.1817303;
    theta = 0.0825398957;
    sigma = 0.0125901;
    r0 = 0.021673;
    
    % Seed for MC reproducibility
    rng(9798)
    
  %% Select and loop over the dimensions
    dims = [10:2:100 105:5:300];
    Qs   = [];
    wces = [];
    Ns   = [];
    zcbs = [];
    MCs  = [];
  
    for d = dims
       
      D = d - 1;
    
      % The true solution
      zcb = zcb_true(kappa, theta, sigma, r0, T, d);
      zcbs = [zcbs zcb];
      
      % Generate the sparse grid (this is the most time-consuming part)
      q  = 2;
      XS = gh_seq(q);
      us = sparse_gens(XS, D);
      us = us(:,2:end); % We do not want the central point
      [Us Ls] = fss_gen(us, D);

      % Integrand evaluations
      Y = zcb_integrand(cell2mat(Us), kappa, theta, sigma, r0, T);    
      
      % Sparse grid kernel quadrature
      l  = d;
      [k kmean Ikmean] = kq_kernel('gauss', l, D, 'normal');
      [Q, wce, wr] = kq_fss(Y, Us, k, kmean, Ikmean, 'true');
      Qs = [Qs Q];
      wces = [wces wce];
      N = sum(Ls);
      Ns = [Ns N];
      
      % Standard MC
      X = randn(D, N);
      Y = zcb_integrand(X, kappa, theta, sigma, r0, T);
      MC = sum(Y)/N;
      MCs = [MCs MC];
      
      % Tell where we are at the moment
      fprintf('Dimension = %i\n', d);
      
    end
  
  %% Plot
  
    subplot(311)
    plot(dims, Qs, dims, MCs, dims, zcbs, '--')
    legend('GHSGKQ', 'Monte Carlo', 'Ground truth')
    title('Integral values')
    
    relErr = abs(Qs - zcbs)./zcbs;
    relErrMC = abs(MCs - zcbs)./zcbs;
    
    subplot(312)
    plot(dims, relErr, dims, relErrMC)
    legend('Relative error (GHSGKQ)', 'Relative error (MC)')
    title('Relative error')  

    subplot(313)
    plot(dims, wces)
    legend('Worst-case error')
    title('Worst-case error')
