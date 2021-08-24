%%%%%%%%%
%%%
%%% EXAMPLE 2:
%%% Integration example with a priori known length-scale
%%%
%%% Toni Karvonen, 2017
%%%
%%%%%%%%%

  %% Initializations
    addpath('../fskq')
    
    % Integrand and kernel
    d = 11;
    l = 0.8;
    
    [k kmean Ikmean] = kq_kernel('gauss', l, d, 'uniform');
    
    xf = linspace(0.2,0.5,d)';
    f = @(x) exp(-norm(x(:)-xf).^2/(2*l^2));
    ground_truth = kmean(xf);
    
    isotropic = 'true';

    % Seed for reproducibility
    rng(163321)

  %% Monte Carlo kernel quadrature with preassigned length-scale
    nKMC = 12000;
    KMCstep = 1000;
    KMCstepn = nKMC/KMCstep;
    KMCs = [];
    KMCwces = [];
    X = [];
    
    for i = 1:KMCstepn
      X = [X, 2*rand(d,KMCstep) - 1];
      Y = funceval(f, X);
      [Q wce] = kq(Y, X, k, kmean, Ikmean, isotropic);
      KMCs = [KMCs Q];
      KMCwces = [KMCwces wce];
    end
  
  %% Sparse grid kernel quadrature with Clenshaw-Curtis points
  % Levels 8 and 9 take some time because billions of kernel evaluations
  % are needed.
    qmax = 9;
    Qs = [];
    wces = [];
    gens = [];
    Ns = [];
    
    for q = 1:qmax
      tic
      XS = cc_seq(q);
      us = sparse_gens(XS, d);
      [Us Ls] = fss_gen(us, d);
      Y = funceval(f, cell2mat(Us));
      [Q, wce, wr] = kq_fss(Y, Us, k, kmean, Ikmean, isotropic);
      Qs = [Qs Q];
      wces = [wces wce];
      Ns = [Ns length(Y)];
      gens = [gens length(Ls)];
      % Some reporting
      fprintf('Level q = %i/%i\n', q, qmax)
      toc
    end
  
  %% Plot
    
    % Compute relative errors errors
    KMCerrs = abs(KMCs(1,:) - ground_truth)/ground_truth;
    FSKQerrs = abs(Qs(1,:) - ground_truth)/ground_truth;
    
    % Plot
    subplot(211)
    loglog(KMCstep:KMCstep:nKMC, KMCerrs, Ns, FSKQerrs, '--x')
    title('Relative erros')
    subplot(212)
    semilogx(KMCstep:KMCstep:nKMC, KMCwces, Ns, wces, '--x')
    title('Kernel quadrature standard deviations')
