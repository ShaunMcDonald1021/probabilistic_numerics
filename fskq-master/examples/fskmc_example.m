%%%%%%%%%
%%%
%%% EXAMPLE 1:
%%% Integration example and comparison with fully symmetric kernel Monte Carlo quadrature
%%%
%%% Toni Karvonen, 2017
%%%
%%%%%%%%%

  %% Initializations
    addpath('../fskq')
  
    % The dimension, integrand and basic kernel without length-scale
    d     = 3;
    f     = @(x) exp( ( sin(5*norm(x)))^2 - x(1)^2 - 0.5*x(2)^2 - 2*x(3)^4 );
    [k kmean Ikmean kwol kdl] = kq_kernel('gauss', 1, d, 'normal');
    isotropic = 'true';
    
    % Seed for reproducibility
    rng(31536)
  
  %% Monte Carlo baseline
    N = 1e7;
    ground_truth = sum(funceval(f, randn(d, N)))/N;
  
  %% Run different methods
    
    warning('off')

    Nmax_FSS = 50; % Maximum number of fully symmetric sets
    N_FSS   = fss_numel(randn(d,1),d); % Number of MC samples per FSS
    
    % For Monte Carlo
    Q_MC = [];
    X_MC = [];
    Y_MC = [];
    
    % For kernel Monte Carlo
    Q_KMC   = [];
    V_KMC   = [];
    
    % For FSS kernel Monte Carlo
    Q_FSKMC   = [];
    V_FSKMC   = [];
    Us         = {};
    X_FSKMC   = [];
    Y_FSKMC   = [];
    
    % Length-scales
    ells  = [];
    l     = 1;
    
    for n = 1:Nmax_FSS
      
      % Monte Carlo
        X = randn(d, N_FSS);
        X_MC = [X_MC X];
        Y_MC = [Y_MC funceval(f, X)];
        Q_MC = [Q_MC sum(Y_MC)/(n*N_FSS)];
      
      % Kernel Monte Carlo
        % Fit the length-scale
        l     = eq_fit(X_MC, Y_MC, kwol, isotropic, l, kdl);
        ells  = [ells l];
        % Quadrature
        [k kmean Ikmean] = kq_kernel('gauss', l, d, 'normal');
        [Q V]   = kq(funceval(f, X_MC), X_MC, k, kmean, Ikmean, isotropic);
        Q_KMC   = [Q_KMC Q];
        V_KMC   = [V_KMC V];
      
      % FSS kernel Monte Carlo
        % Generate new MC FSS
        U         = fss(randn(d, 1), d);
        Us{n}     = U;
        % Quadrature
        X_FSKMC  = [X_FSKMC U];
        Y_FSKMC  = [Y_FSKMC funceval(f, U)];
        [Q, V]    = kq_fss(funceval(f, cell2mat(Us)), Us, k, kmean, Ikmean, isotropic);
        Q_FSKMC  = [Q_FSKMC Q];
        V_FSKMC  = [V_FSKMC V];
      
      % Display progress
        fprintf('Step %i/%i\n', n, Nmax_FSS)
          
      end
      
      warning('on')

  %% Plot
    
    Ns = 1:Nmax_FSS;
    
    subplot(311)
    plot(Ns, Q_MC, Ns, Q_KMC, Ns, Q_FSKMC)
    legend('MC', 'KMC', 'FSKMC')
    title('Quadrature estimates')
    subplot(312)
    plot(Ns,  sqrt(V_KMC), Ns, sqrt(V_FSKMC))
    legend('KMC', 'FSKMC')
    title('Standard deviations')  
    subplot(313)
    plot(Ns, ells)
    legend('KMC length-scale')
    title('Fitted length-scale')
