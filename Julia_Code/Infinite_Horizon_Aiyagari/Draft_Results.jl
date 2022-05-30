###################################################################
###################################################################
###################################################################
## Graphs and Tables for Draft 




###################################################################
###################################################################
## Histogram 


###################################################################
## Top Wealth Shares from Histogram 
function H_Moments_Top_Shares(M::Model,pct_list)
    
    # Allocate results 
        Wealth_Stats = zeros(length(pct_list)+1) ;
        Wealth_Share = zeros(length(pct_list)  ) ;

    # Compute distribution of assets 
        Γ_a   = dropdims( sum( M.Γ , dims=(3,2) ) , dims=(3,2) ) ; # Assets 
        CDF_a = cumsum(Γ_a) ;
    
    # Compute moments from distribution 
        Wealth_Stats[end] = sum( M.a_grid_fine.*Γ_a ) ; # Average 
        for j=1:length(pct_list)
            ind_top = CDF_a.>=pct_list[j]/100 ; 
            Wealth_Stats[j] = M.a_grid_fine[ind_top][1] ;
            Wealth_Share[j] = 100*sum( M.a_grid_fine[ind_top].*Γ_a[ind_top] )/Wealth_Stats[end] ; 
        end

    # Pareto Coefficient 
        ind     = M.a_grid_fine.>=1000 ;
        grid_1M = M.a_grid_fine[ind]   ;
        Γ_a_1M  = Γ_a[ind]/sum(Γ_a[ind])        ; Γ_a_1M = Γ_a_1M/sum(Γ_a_1M) ; 
        CCDF_1M = 1 .- cumsum(Γ_a_1M)           ;
        P_coeff = (log.(grid_1M[1:end-1]./1000)'*log.(grid_1M[1:end-1]./1000))\log.(grid_1M[1:end-1]./1000)'*log.(CCDF_1M[1:end-1])

    # Return Moments 
        return Wealth_Stats, Wealth_Share, P_coeff
end 



###################################################################
## Transition Rates Across Wealth Deciles
function H_Moments_Decile_Transitions(M::Model,N::Int)
    
    # Allocate results 
        Tr_deciles_a = zeros(10,10) ;

    # Compute distribution of assets 
        Γ_a   = dropdims( sum( M.Γ , dims=(3,2) ) , dims=(3,2) ) ; # Assets 
        CDF_a = cumsum(Γ_a) ;

    # Find deciles 
        deciles_a     = Array{Int64}(undef,11)   ;
        deciles_l     = Array{Float64}(undef,11) ;
        deciles_a[1]  = 0 ; deciles_l[1] = M.a_grid_fine[1] ;
        for i=2:10 
        deciles_a[i]  = collect(1:M.n_a_fine)[(100*CDF_a).>=(10*(i-1))][1] ;
        deciles_l[i]  = M.a_grid_fine[deciles_a[i]]                      ;
        end 
        deciles_a[11] = M.n_a_fine ; deciles_l[11] = M.a_grid_fine[end]  ;
        deciles_mat = [deciles_a deciles_l] ;

    # For each decile compute transitions   
        for i=1:10
            # Conditional distribution on current decile 
            Γ_0                                    = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ)          ; 
            Γ_0[deciles_a[i]+1:deciles_a[i+1],:,:] = M.Γ[deciles_a[i]+1:deciles_a[i+1],:,:] ; 
            Γ_0 = Γ_0/sum(Γ_0) ; 
            # Iterate distribution 
            Γ_N = Histogram_Iteration(M,n_H,Γ_0) ;
            # Obtain transitions (sum over new distribution within deciles)
            Tr_deciles_a[i,:] = [sum( Γ_N[deciles_a[j]+1:deciles_a[j+1],:,:] ) for j in 1:10] ;
        end 

    # Return Moments 
        return Tr_deciles_a, deciles_mat
end 



###################################################################
## 3 Year autocorrelation of consumption 
function H_Moments_C_Correlation(M::Model,N::Int,a_min,a_max)

    # Get average consumption and standard deviation 
        # First quintile distribution 
        Γ_q                     = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ)   ; 
        Γ_q[a_min:a_max,:,:]    = M.Γ[a_min:a_max,:,:]            ; 
        Γ_q                     = Γ_q/sum(Γ_q) ; 
        # First quintile consumption 
        av_c_q = sum( M.G_c_fine[a_min:a_max,:,:].*Γ_q[a_min:a_max,:,:] ) ;
        sd_c_q = sqrt( sum( ((M.G_c_fine[a_min:a_max,:,:].-av_c_q).^2).*Γ_q[a_min:a_max,:,:] ) ) ;
        # Future average and standard deviation conditional on first quintile
            # Iterate distribution 
            Γ_qN = Histogram_Iteration(M,N,Γ_q) ;
            # Obtain consumption
            av_c_q_N = sum( M.G_c_fine.*Γ_qN  )                                   ;
            sd_c_q_N = sqrt( sum( ((M.G_c_fine.-av_c_q_N).^2).*Γ_qN  ) ) ; 

    # Compute integrand of correlation for first and tenth deciles  
        cov_c = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
        for i_ζ=1:M.n_ζ # Current ζ
        for i_ϵ=1:M.n_ϵ # Current ϵ
            for i_a=a_min:a_max # Current a
                # For each state in state space (a,ϵ,ζ) get conditional distribution 
                Γ_0 = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
                Γ_0[i_a,i_ϵ,i_ζ] = 1 ; 
                # Iterate distribution
                Γ_N = Histogram_Iteration(M,N,Γ_0) ;
                # Get portion of integrand (C_0 - av_c)*(C_n - av_c)*Γ_n
                cov_c[i_a,i_ϵ,i_ζ] = sum( ( M.G_c_fine[i_a,i_ϵ,i_ζ].-av_c_q )*( M.G_c_fine.-av_c_q_N ).*Γ_N ) ;
            end 

        end 
        end 

    # Integrate with respect to initial distribution of first quintile
        cov_c_q  = sum( cov_c.*Γ_q )  ;
        cor_c_q  = cov_c_q/sqrt(sd_c_q^2*sd_c_q_N^2) ; 

    # Return Moment 
        return cor_c_q 
end 


###################################################################
## Run Histogram Simulation for Different Grids 
    H_grid_size = 100:50:1000 ; 
    n_H = length(H_grid_size) ;
    H_Γ_timed = zeros(n_H)    ; 
    H_Γ_bytes = zeros(n_H)    ;
    H_M_timed = zeros(n_H,3)  ; 
    H_M_bytes = zeros(n_H,3)  ;

    H_Wealth_Stats = zeros(n_H,6) ; 
    H_Wealth_Share = zeros(n_H,5) ; 
    H_Pareto_Coeff = zeros(n_H  ) ; 
    pct_list = [90;95;99;99.9;99.99] ; 

    H_Decile       = zeros(11,2 ,n_H) ;
    H_Decile_Tr    = zeros(10,10,n_H) ;
    N_Decile_Tr    = 9                ;

    H_Cons_Corr    = zeros(n_H)       ;
    N_Cons_Corr    = 2                ;

    for i=1:n_H 

        println("   Histogram with $(H_grid_size[i]) Grid Points")

        # Set up model structure 
        M_Hist = Model(n_a_fine=H_grid_size[i],method=1,read_flag=true) ;
        # Solve for stationary distribution and save time and allocation - Adjust grid size
        M_Hist, H_Γ_timed[i], H_Γ_bytes[i] = @timed Aiyagari_Equilibrium(M_Hist);
        
       
        ## Moments 
        
        # 1-2) Top Wealth Shares and Pareto Coefficient 
        out, time, memory = @timed H_Moments_Top_Shares(M_Hist,pct_list) ;
        H_Wealth_Stats[i,:] = out[1] ; H_Wealth_Share[i,:] = out[2] ; H_Pareto_Coeff[i] = out[3] ; H_M_timed[i,1] = time ; H_M_bytes[i,1] = memory ; 


        # 3) Decile Transitions 
        out, time, memory = @timed H_Moments_Decile_Transitions(M_Hist,N_Decile_Tr) ;
        H_Decile_Tr[:,:,i] = out[1] ; H_Decile[:,:,i] = out[2] ; H_M_timed[i,2] = time ;  H_M_bytes[i,2] = memory ; 
        
        
        # 4) Consumption Autocorrelation for first quintile
        out, time, memory = @timed H_Moments_C_Correlation(M_Hist,N_Cons_Corr,1,Int(H_Decile[3,1,i])) ;
        H_Cons_Corr[i] = out ;  H_M_timed[i,3] = time ; H_M_bytes[i,3] = memory ;

    end 



###################################################################
###################################################################
## Simulation 




###################################################################
## Run Simulation for Different Panel Size 

    # Set up model structures 
    M_Simul = Model(method=1,read_flag=true) ;
    M_Panel = Model_Panel(N_Panel=10000)   ; 

    # Set up discrete observations 
    S_sample = 1000:1000:M_Panel.N_Panel ; 
    N_S      = length(S_sample)          ;
    pct_list = [90;95;99;99.9;99.99]     ;  

    S_M_timed      = zeros(n_S,3)     ; 
    S_M_bytes      = zeros(n_S,3)     ;
    
    S_Wealth_Stats = zeros(n_S,6)     ; 
    S_Wealth_Share = zeros(n_S,5)     ; 
    H_Pareto_Coeff = zeros(n_S  )     ;  

    H_Decile       = zeros(11,n_S)    ;
    H_Decile_Tr    = zeros(10,10,n_S) ;

    H_Cons_Corr    = zeros(n_S)       ;
    
    # Solve model 
    M_Simul, S_Γ_timed, S_Γ_bytes = @timed Aiyagari_Equilibrium(M_Simul);

    # Simulate Panel 
    M_Panel = Simulate_Panel_Dynasty(M_Simul,M_Panel) ; 

    
    ## Moments 

    # 1-2) Top Wealth Shares and Pareto Coefficient 
    for i=1:N_S 
        tic()

        # Select sample 
        a_sample  = M_Panel.a_mat[1:S_sample[i],end] ;

        # Average wealth 
        S_Wealth_Stats[i,end]     = mean( a_sample )   ;

        # Percentiles 
        S_Wealth_Stats[i,1:end-1] = percentile( a_sample , pct_list ) ;

        # Top Shares 
        S_Wealth_Share[i,:] = [ 100*sum( a_sample[ a_sample.>=S_Wealth_Stats[i,p] ]  )/(S_sample[i]*S_Wealth_Stats[i,end])  for p in 1:5] ;

        # Pareto Coefficient 
        Pareto_sample = sort( a_sample[ a_sample.>=1000 ] ) ; # Select and sort observations above $1M  
        Pareto_CCDF   = collect(length(Pareto_sample):-1:1)./length(Pareto_sample) ; # Counter CDF = 1- CDF
        H_Pareto_Coeff[i] = (log.(Pareto_sample[1:end-1]./1000)'*log.(Pareto_sample[1:end-1]./1000))\log.(Pareto_sample[1:end-1]./1000)'*log.(Pareto_CCDF[1:end-1]) ;

        # Time it 
        S_M_timed[i,1] = toc() ; 
    end 


    # 3) Decile Transitions 
    for i=1:N_S 
        tic()

        a_aux_0       = M_Panel.a_mat[1:S_sample[i],1]   ;
        a_aux_T       = M_Panel.a_mat[1:S_sample[i],end] ;
        H_Decile[:,i] = percentile( a_aux_T , 0:10:100 ) ; # Deciles based on end of sample (hoping for stationariety)
        for p=1:10
            ind_d = findall(x-> H_Decile[p,i]<=x<=H_Decile[i,p+1], a_aux_0 ) ; # Find "i" in each decile in t=1
            a_T   = a_aux_T[ind_d] ; # Follow them to t=10
            H_Decile_Tr[p,:,i] = [100*sum( H_Decile[pT,i].<=a_T.<=H_Decile[pT+1,i] )/length(ind_d) for pT in 1:10] ; # Get transition rates 
        end 

        # Time it 
        S_M_timed[i,2] = toc() ; 
    end 
    
    
    # 4) Consumption Autocorrelation for first quintile
    for i=1:N_S 
        tic()

        # Fix current sample 
        a_aux_0       = M_Panel.a_mat[1:fig_sample[i],8]   ;
        c_aux_0       = M_Panel.c_mat[1:fig_sample[i],8]   ;
        c_aux_T       = M_Panel.c_mat[1:fig_sample[i],end] ;
        # Find index of first quintile 
        ind_q         = findall(x-> 0<=x<=S_Decile[3,i], a_aux_0 ) ; # Index of first quintile of assets in T-1
        # Compute moments 
        S_Cons_Corr[i]= cor( [c_aux_0[ind_q] c_aux_T[ind_q]]  ;dims=1)[2]              ;  

        # Time it 
        S_M_timed[i,3] = toc() ; 
    end 






###################################################################
###################################################################
## Graphs and Tables 
    


###################################################################
## Top 1%, 0.1%, 0.01% Shares 
    # Top 1%
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  H_M_timed[:,1] , H_Wealth_Share[:,3] ,label="Top 1% - Histogram")
    plot!( S_M_timed[:,1] , S_Wealth_Share[:,3] ,label="Top 1% - Simulation")
    ylims!(0,1)
    # xlims!(1,M_P.N_Panel/1000)
    title!("Top 1% Share of Wealth",titlefont=14)
    xlabel!("Time: Seconds",labelsize=18)
    savefig("./"*Fig_Folder*"/Draft_Top_1_Share.pdf")

    # Top 0.1%
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  H_M_timed[:,1] , H_Wealth_Share[:,4] ,label="Top 0.1% - Histogram")
    plot!( S_M_timed[:,1] , S_Wealth_Share[:,4] ,label="Top 0.1% - Simulation")
    ylims!(0,1)
    # xlims!(1,M_P.N_Panel/1000)
    title!("Top 0.1% Share of Wealth",titlefont=14)
    xlabel!("Time: Seconds",labelsize=18)
    savefig("./"*Fig_Folder*"/Draft_Top_01_Share.pdf")

    # Top 0.01%
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  H_M_timed[:,1] , H_Wealth_Share[:,5] ,label="Top 0.01% - Histogram")
    plot!( S_M_timed[:,1] , S_Wealth_Share[:,5] ,label="Top 0.01% - Simulation")
    ylims!(0,1)
    # xlims!(1,M_P.N_Panel/1000)
    title!("Top 0.01% Share of Wealth",titlefont=14)
    xlabel!("Time: Seconds",labelsize=18)
    savefig("./"*Fig_Folder*"/Draft_Top_001_Share.pdf")


###################################################################
## Pareto Tail 
    # Results from Histogram 
    ind     = M_Simul.a_grid_fine.>=1000 ;
    grid_1M = M_Simul.a_grid_fine[ind]   ;
    Γ_a_1M  = Γ_a[ind]/sum(Γ_a[ind])     ; Γ_a_1M = Γ_a_1M/sum(Γ_a_1M) ; 
    CCDF_1M = 1 .- cumsum(Γ_a_1M)        ;
    P_coeff = (log.(grid_1M[1:end-1]./1000)'*log.(grid_1M[1:end-1]./1000))\log.(grid_1M[1:end-1]./1000)'*log.(CCDF_1M[1:end-1])

    # 10k Simulation 
    Pareto_a_10k = M_Panel.a_mat[1:10000,end]                  ; # Select first 10.000 observations 
    Pareto_a_10k = sort( Pareto_a_10k[ Pareto_a_10k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_10k = collect(length(Pareto_a_10k):-1:1)./length(Pareto_a_10k) ; # Counter CDF = 1- CDF
    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α=$(round(P_coeff,digits=2))",12)
    scatter!( log.(Pareto_a_10k[1:end-1]./1000)  , log.(Pareto_p_10k[1:end-1])  , marker=(:circle ,3)    , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[10].*log.(grid_1M[1:end-1]./1000) , w=2 , label=nothing )
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Draft_Pareto_10k.pdf")

    # 50k Simulation 
    Pareto_a_50k = M_Panel.a_mat[1:50000,end]                  ; # Select first 50.000 observations 
    Pareto_a_50k = sort( Pareto_a_50k[ Pareto_a_50k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_50k = collect(length(Pareto_a_50k):-1:1)./length(Pareto_a_50k) ; # Counter CDF = 1- CDF
    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α=$(round(P_coeff,digits=2))",12)
    scatter!( log.(Pareto_a_50k[1:end-1]./1000)  , log.(Pareto_p_50k[1:end-1])  , marker=(:circle ,3)    , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[50].*log.(grid_1M[1:end-1]./1000) , w=2 , label=nothing )
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Draft_Pareto_50k.pdf")


    # 100k Simulation 
    Pareto_a_100k = M_Panel.a_mat[1:100000,end]                  ; # Select first 100.000 observations 
    Pareto_a_100k = sort( Pareto_a_100k[ Pareto_a_100k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_100k = collect(length(Pareto_a_100k):-1:1)./length(Pareto_a_100k) ; # Counter CDF = 1- CDF
    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α=$(round(P_coeff,digits=2))",12)
    scatter!( log.(Pareto_a_100k[1:end-1]./1000)  , log.(Pareto_p_100k[1:end-1])  , marker=(:circle ,3)    , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[100].*log.(grid_1M[1:end-1]./1000) , w=2 , label=nothing )
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Draft_Pareto_100k.pdf")

    # 500k Simulation 
    Pareto_a_500k = M_Panel.a_mat[1:500000,end]                  ; # Select first 500.000 observations 
    Pareto_a_500k = sort( Pareto_a_500k[ Pareto_a_500k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_500k = collect(length(Pareto_a_500k):-1:1)./length(Pareto_a_500k) ; # Counter CDF = 1- CDF
    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α=$(round(P_coeff,digits=2))",12)
    scatter!( log.(Pareto_a_500k[1:end-1]./1000)  , log.(Pareto_p_500k[1:end-1])  , marker=(:circle ,3)    , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[500].*log.(grid_1M[1:end-1]./1000) , w=2 , label=nothing )
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Draft_Pareto_500k.pdf")


    # 1M Simulation 
    Pareto_a_1M = M_Panel.a_mat[1:1000000,end]                  ; # Select first 1.000.000 observations 
    Pareto_a_1M = sort( Pareto_a_1M[ Pareto_a_1M.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_1M = collect(length(Pareto_a_1M):-1:1)./length(Pareto_a_1M) ; # Counter CDF = 1- CDF
    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α=$(round(P_coeff,digits=2))",12)
    scatter!( log.(Pareto_a_1M[1:end-1]./1000)  , log.(Pareto_p_1M[1:end-1])  , marker=(:circle ,3)    , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[1000].*log.(grid_1M[1:end-1]./1000) , w=2 , label=nothing )
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Draft_Pareto_1M.pdf")


###################################################################
## Decile Transition
    # Persistence first decile 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  H_M_timed[:,2] , 100*dropdims(H_Decile_Tr[1,1,:]) ,label="d1-d1 Transition - Histogram")
    plot!( S_M_timed[:,2] , 100*dropdims(S_Decile_Tr[1,1,:]) ,label="d1-d1 Transition - Simulation")
    ylims!(0,1)
    # xlims!(1,M_P.N_Panel/1000)
    # title!("Top 1% Share of Wealth",titlefont=14)
    xlabel!("Time: Seconds",labelsize=18)
    savefig("./"*Fig_Folder*"/Draft_Decile_Tr_1.pdf")

    # Persistence tenth decile 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  H_M_timed[:,2] , 100*dropdims(H_Decile_Tr[10,10,:]) ,label="d10-d10 Transition - Histogram")
    plot!( S_M_timed[:,2] , 100*dropdims(S_Decile_Tr[10,10,:]) ,label="d10-d10 Transition - Simulation")
    ylims!(0,1)
    # xlims!(1,M_P.N_Panel/1000)
    # title!("Top 1% Share of Wealth",titlefont=14)
    xlabel!("Time: Seconds",labelsize=18)
    savefig("./"*Fig_Folder*"/Draft_Decile_Tr_10.pdf")


###################################################################
## Consumption 3 Year Auto-Correlation 
    Mat = [ H_Cons_Corr[1]  H_Cons_Corr[4]  H_Cons_Corr[9]   S_Cons_Corr[10]                 S_Cons_Corr[100]                 S_Cons_Corr[500] ;
            H_Γ_timed[1]    H_Γ_timed[4]    H_Γ_timed[9]     S_Γ_timed+M_Panel.t_vec[10000]  S_Γ_timed+M_Panel.t_vec[100000]  S_Γ_timed+M_Panel.t_vec[500000] ;
            H_Timed[1,3]    H_Timed[4,3]    H_Timed[9,3]     S_M_timed[10,3]                 S_M_timed[100,3]                 S_M_timed[500,3]   ] ;
