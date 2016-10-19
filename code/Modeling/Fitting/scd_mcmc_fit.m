function [LOOP] = scd_mcmc_fit(compute_proba,param0,sigma_pert,paramlimits,NN)
%% 
    % scd_mcmc_fit(compute_proba,param0,sigma_pert,paramlimits,NN)
    %
    % Inputs : ------------------------------------------------------------------------------------------------
    % compute_proba     function    with which the likelihood will be calculated at each iteration    
    % param0            vector      initial parameters
    % sigma_pert        vector      std deviations of perturbation
    % param_limits      matrix      Limits within which the parameters should belong : [mins ; maxs]
    % NN                Number      Number of iterations
    % ---------------------------------------------------------------------------------------------------------
    %
    % Example :
    % scd_mcmc_fit(@(param) scd_compute_proba(Ax,param),[0.8 0.5 3],[0.05 0.1 0.25],[ 0 0 1 ; 1 3 15 ],2000);
        
    %% Beginning of the loop

    LOOP.welcome=0;
    LOOP.allparam   = zeros(NN+1,size(param0,2));
    LOOP.parameters = zeros(NN+1,size(param0,2));
    LOOP.logLnewvec = zeros(NN+1,1);
    LOOP.logL       = zeros(NN+1,1);
    LOOP.ratio      = zeros(NN+1,1);
    
for i = 1 : NN+1      % New Random parameters and model curve
    %
j_progress(i/(NN+1));

    if i==1
    LOOP.parameters(i,:) = param0 ;
    LOOP.allparam(i,:)   = param0 ;
    else
        % Perturbating parameters with the right noise
        LOOP.allparam(i,:) = random('Rician',LOOP.parameters(i-1,:),sigma_pert) ;
        %j_progress(i/NN)
    end
       
    %% Calculating new likelihood
    
    % Verifying that basic conditions are satisfied
    for j=1:size(param0,2)
        in_limits_param = paramlimits(1,j) <= LOOP.allparam(i,j) <= paramlimits(2,j); 
    end
    
    % Calculating new Likelihood function %%%%%%%----------------------------------------------------------------------
    if in_limits_param           
        LOOP.logLnew = compute_proba(LOOP.allparam(i,:));
    else
        LOOP.logLnew=-Inf; 
    end

    LOOP.logLnewvec(i) = LOOP.logLnew;
    
    %%  Core of the MCMC loop
    
    if i==1
        ratio = 1; 
    else
        ratio = min(1,exp((LOOP.logLnew-LOOP.logL(i-1)))); % pnew/p(i-1)
    end
    
    LOOP.ratio(i) = ratio ;
    
    if  rand < ratio ; % Then Update parameters %%%%%%%------------------------
        LOOP.parameters(i,:) = LOOP.allparam(i,:);
        LOOP.logL(i)=LOOP.logLnew;
        LOOP.welcome=LOOP.welcome+1;
        
    else
        LOOP.parameters(i,:) = LOOP.parameters(i-1,:); % Otherwise keep previous fh
        LOOP.logL(i) = LOOP.logL(i-1) ;
    end
       
    %%--------------------------------------------------------------------------------------------------------------------------------------
    
end