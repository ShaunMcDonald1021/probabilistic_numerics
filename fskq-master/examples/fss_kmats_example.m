%%%%%%%%%
%%%
%%% Kernel matrix illustrations with fully symmetric sets
%%%
%%% Toni Karvonen, 2017
%%%
%%%%%%%%%

  %% Initializations
    addpath('../fskq')
    
    % Global parameters
    l = 1;
    k = @(r) exp(-r.^2/(2*l^2));
  
  %% 1st matrix
    d = 3;
    us = [0 0.5 1 2 ; 
          0 0.4 1 1 ; 
          0 0.1 1 2 ];
    Us = fss_gen(us, d);
    X = cell2mat(Us);
    K1 = kmat(X, k);
  
  %% 2nd matrix
    d = 4;
    us = [0.05 1 2 0.5 0.5;
          0 1 1 0.4 0.5; 
          0 1 2 0.1 1];
    Us = fss_gen(us, d);
    X = cell2mat(Us);
    K2 = kmat(X, k);
  
  %% 3rd matrix (random)
    rng(78)
    d = 4;
    N = 200;
    X = randn(d, N);
    K3 = kmat(X, k);
  
  %% Plot
    figure
    imagesc(K1)
    axis equal tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])  
    set(gca,'position',[0 0 1 1],'units','normalized')
    
    figure
    imagesc(K2)
    axis equal tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])    
    set(gca,'position',[0 0 1 1],'units','normalized')

    figure
    imagesc(K3)
    axis equal tight
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])    
    set(gca,'position',[0 0 1 1],'units','normalized')  

