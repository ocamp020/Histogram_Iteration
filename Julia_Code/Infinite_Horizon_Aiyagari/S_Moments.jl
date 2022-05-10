###################################################################
###################################################################
###################################################################
## Compute moments using the histogram iteration method 


###################################################################
###################################################################
## Add Simulation Functions + 
    include("Simulation_Panel.jl")

## Run Simulations
    M_P = Model_Panel() ; 
    M_P = Simulate_Panel(M_Aiyagari,M_P) ;

    fig_sample = M_P.N_Min:1000:M_P.N_Panel ; 
    fig_N      = length(fig_sample)         ; 

###################################################################
###################################################################
## Average Wealth, Percentiles, and Top Wealth Shares 
    pct_list     = [90;95;99;99.9;99.99] ;
    av_a_S       = zeros(fig_N)   ;
    pct_S        = zeros(5,fig_N) ;
    Top_Shares_S = zeros(5,fig_N) ; 
    for i = 1:fig_N
        a_aux       = M_P.a_mat[1:fig_sample[i],end]   ;
        av_a_S[i]   = sum( a_aux )/i     ;
        pct_S[:,i]  = percentile( a_aux , pct_list ) ;
        Top_Shares_S[:,i] = [ 100*sum( a_aux[ a_aux.>=pct_S[p] ]  )/(fig_sample[i]*av_a_S[i])  for p in 1:5] ;
    end 

## Figures: Top Shares 1% and 0.1%     


###################################################################
###################################################################
## Lorenz Curve 
    Lorenz_S = zeros(100,3) ; 
    sample_vec = [10000, 50000, 100000] ; 
    for i=1:3
        a_aux       = M_P.a_mat[1:sample_vec[i],end]   ; av_a_aux = sum(a_aux) ;
        p_aux       = percentile( a_aux , 1:100 )    ;
        Lorenz_S[:,i] = [ 100*sum( a_aux[ a_aux.>=p_aux[p] ]  )/((M_P.N_Min+i)*av_a_aux)  for p in 1:100] ; 
    end 

## Figure with Lorenz Curve for 10k, 50k, 100k    



###################################################################
###################################################################
## Pareto Tail 
    # 10.000 observations 
    Pareto_a_10k = M_P.a_mat[1:10000,end]                      ; # Select first 10.000 observations 
    Pareto_a_10k = sort( Pareto_a_10k[ Pareto_a_10k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_10k = collect(length(Pareto_a_10k):-1:1)./length(Pareto_a_10k) ; # Counter CDF = 1- CDF

    # 50.000 observations 
    Pareto_a_50k = M_P.a_mat[1:50000,end]                      ; # Select first 50.000 observations 
    Pareto_a_50k = sort( Pareto_a_50k[ Pareto_a_50k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_50k = collect(length(Pareto_a_50k):-1:1)./length(Pareto_a_50k) ; # Counter CDF = 1- CDF

    # 100.000 observations 
    Pareto_a_100k = M_P.a_mat[1:100000,end]                      ; # Select first 100.000 observations 
    Pareto_a_100k = sort( Pareto_a_100k[ Pareto_a_100k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_100k = collect(length(Pareto_a_100k):-1:1)./length(Pareto_a_100k) ; # Counter CDF = 1- CDF

## Figure with all pareto tails 



###################################################################
###################################################################
## 10 Year transition rates for decile 
    deciles_a_S     = zeros(11,fig_N)    ;
    Tr_decicles_a_S = zeros(10,10,fig_N) ;
    for i=1:fig_N
        a_aux_0          = M_P.a_mat[1:fig_sample[i],1]   ;
        a_aux_T          = M_P.a_mat[1:fig_sample[i],end] ;
        deciles_a_S[:,i] = percentile( a_aux_T , 0:10 )   ; # Deciles based on end of sample (hoping for stationariety)
        for p=1:10
            ind_d = findall(x-> deciles_a_S[p]<=x<=deciles_a_S[p+1], a_aux_0 ) ; # Find "i" in each decile in t=1
            a_T = a_aux_T[ind_d] ; # Follow them to t=10
            Tr_decicles_a_S[p,:,i] = [100*sum( deciles_a_S[pT].<=a_T.<=deciles_a_S[pT+1] )/length(ind_d) for pT in 1:10] ; # Get transition rates 
        end 
    end 

## Figure with transitions of bottom decile


## Figure with transitions of top decile


###################################################################
###################################################################
## 5 Year autocorrelation of consumption first quintile
    av_c_q_S  = zeros(2,fig_N)    ; # Average Consumption in first quintile 
    sd_c_q_S  = zeros(2,fig_N)    ; # Std Dev of Consumption in first quintile 
    cor_c_q_S = zeros(fig_N)      ; # Cor of Consumption in first quintile 
    for i=1:fig_N
        # Fix current sample 
        a_aux_0       = M_P.a_mat[1:fig_sample[i],9]   ;
        c_aux_0       = M_P.c_mat[1:fig_sample[i],9]   ;
        c_aux_T       = M_P.c_mat[1:fig_sample[i],end] ;
        # Find index of first quintile 
        ind_q         = findall(x-> 0<=x<=deciles_a_S[3], a_aux_0 ) ; # Index of first quintile of assets in T-1
        # Compute moments 
        av_c_q_S[1,i] = mean( c_aux_0[ind_q] ) ; sd_c_q_S[1,i] = std( c_aux_0[ind_q] ) ;
        av_c_q_S[2,i] = mean( c_aux_T[ind_q] ) ; sd_c_q_S[2,i] = std( c_aux_T[ind_q] ) ;
        cor_c_q_S[i]  = cor( [c_aux_0[ind_q] c_aux_T[ind_q]]  ;dims=1)[2]              ;  
    end 
    
## Figure     