###################################################################
###################################################################
###################################################################
## Compute moments using the histogram iteration method 


###################################################################
###################################################################
## Add Simulation Functions
    include("Simulation_Panel.jl")

## Run Simulations
    M_P = Model_Panel() ; 
    M_P = Simulate_Panel(M_Aiyagari,M_P) ;

    fig_sample = [1000 ; collect(5000:5000:M_P.N_Panel)] ; 
    fig_N      = length(fig_sample)         ; 

    # Save Results 
    open(File_Folder*"/S_a_mat.csv", "w") do io
        writedlm(io, M_P.a_mat , ',')
    end;
    open(File_Folder*"/S_c_mat.csv", "w") do io
        writedlm(io, M_P.c_mat , ',')
    end;
    open(File_Folder*"/S_e_mat.csv", "w") do io
        writedlm(io, M_P.ϵ_mat , ',')
    end;
    open(File_Folder*"/S_z_mat.csv", "w") do io
        writedlm(io, M_P.ζ_mat , ',')
    end;
    

###################################################################
###################################################################
## Average Wealth, Percentiles, and Top Wealth Shares 
    pct_list     = [90;95;99;99.9;99.99] ;
    av_a_S       = zeros(fig_N)   ;
    pct_S        = zeros(5,fig_N) ;
    Top_Shares_S = zeros(5,fig_N) ; 
    for i = 1:fig_N
        a_sample       = M_P.a_mat[1:fig_sample[i],end] ;
        av_a_S[i]   = mean( a_sample )                  ;
        pct_S[:,i]  = percentile( a_sample , pct_list ) ;
        Top_Shares_S[:,i] = [ 100*sum( a_sample[ a_sample.>=pct_S[p] ]  )/(fig_sample[i]*av_a_S[i])  for p in 1:5] ;
    end 

## Figures: Top percentiles 1% and 0.1%     
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(  fig_sample/1000 , pct_S[3,:] , marker=(:circle ,3,:cornflowerblue),label=nothing)
    hline!( [Top_shares[3,2]] ,c=:orange  , w=3 , label=nothing )
    plot!( fig_sample/1000 , pct_S[4,:] , marker=(:circle ,3,:cornflowerblue),label=nothing)
    hline!( [Top_shares[4,2]] ,c=:orange  , w=3 , label=nothing )
    ylims!(floor(minimum(pct_S[3,:]/500))*500,ceil(maximum(pct_S[4,:]/500))*500)
    xlims!(1,M_P.N_Panel/1000)
    title!("99th and 99.9th percentiles",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/Top_pct_Simul.pdf")


## Figures: Top Shares 1% and 0.1%     
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( fig_sample/1000 , Top_Shares_S[3,:] , marker=(:circle ,3,:cornflowerblue),label=nothing)
    hline!( [Top_shares[3,3]] ,c=:orange  , w=3 , label=nothing )
    plot!( fig_sample/1000 , Top_Shares_S[4,:] , marker=(:circle ,3,:cornflowerblue),label=nothing)
    hline!( [Top_shares[4,3]] ,c=:orange  , w=3 , label=nothing )
    ylims!(0,ceil(maximum(Top_Shares_S[3,:]/5))*5)
    xlims!(1,M_P.N_Panel/1000)
    title!("Top 1% and 0.1% Shares",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/Top_Share_Simul.pdf")


###################################################################
###################################################################
## Lorenz Curve 
    Lorenz_S = zeros(100,4) ; 
    sample_vec = [10000, 50000, 100000, 500000] ; 
    for i=1:4
        a_sample       = M_P.a_mat[1:sample_vec[i],end]   ; av_a_aux = sum(a_sample) ;
        p_aux       = percentile( a_sample , 1:100 )    ;
        Lorenz_S[:,i] = [ 100*sum( a_sample[ a_sample.<=p_aux[p] ]  )/av_a_aux  for p in 1:100] ; 
    end 

## Figure with Lorenz Curve for 10k, 50k, 100k    
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(1:100,1:100,w=1,linecolor=:gray70,label=nothing,aspect_ratio=1)
    plot!( 100*CDF_a  , 100*Lorenz_a  , w=3, c=:cornflowerblue , label=nothing ,aspect_ratio=1)   
    plot!( 1:100      , Lorenz_S[:,1] , w=2 , label=nothing ,aspect_ratio=1)   
    plot!( 1:100      , Lorenz_S[:,2] , w=2 , label=nothing ,aspect_ratio=1)   
    plot!( 1:100      , Lorenz_S[:,3] , w=2 , label=nothing ,aspect_ratio=1)   
    plot!( 1:100      , Lorenz_S[:,4] , w=2 , label=nothing ,aspect_ratio=1)   
    title!("Lorenz curve",titlefont=14)
    ylims!(0,100); xlims!(0,100)
    savefig("./"*Fig_Folder*"/Distribution_Wealth_Lorenz_Simul.pdf")


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
    Pareto_a_100k = M_P.a_mat[1:100000,end]                       ; # Select first 100.000 observations 
    Pareto_a_100k = sort( Pareto_a_100k[ Pareto_a_100k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_100k = collect(length(Pareto_a_100k):-1:1)./length(Pareto_a_100k) ; # Counter CDF = 1- CDF

    # 500.000 observations 
    Pareto_a_500k = M_P.a_mat[1:500000,end]                       ; # Select first 500.000 observations 
    Pareto_a_500k = sort( Pareto_a_500k[ Pareto_a_500k.>=1000 ] ) ; # Select and sort observations above $1M  
    Pareto_p_500k = collect(length(Pareto_a_500k):-1:1)./length(Pareto_a_500k) ; # Counter CDF = 1- CDF

## Figure with all pareto tails 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    scatter!( log.(grid_1M[1:end-1]./1000)       , log.(CCDF_1M[1:end-1])       , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    scatter!( log.(Pareto_a_10k[1:end-1]./1000)  , log.(Pareto_p_10k[1:end-1])  , marker=(:circle ,3)                 , markerstrokewidth=0 , label=nothing )   
    scatter!( log.(Pareto_a_50k[1:end-1]./1000)  , log.(Pareto_p_50k[1:end-1])  , marker=(:circle ,3)                 , markerstrokewidth=0 , label=nothing )   
    scatter!( log.(Pareto_a_100k[1:end-1]./1000) , log.(Pareto_p_100k[1:end-1]) , marker=(:circle ,3)                 , markerstrokewidth=0 , label=nothing )   
    scatter!( log.(Pareto_a_500k[1:end-1]./1000) , log.(Pareto_p_500k[1:end-1]) , marker=(:circle ,3)                 , markerstrokewidth=0 , label=nothing )   
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Distribution_Wealth_Pareto_Simul.pdf")


###################################################################
###################################################################
## 10 Year transition rates for decile 
    deciles_a_S     = zeros(11,fig_N)    ;
    Tr_decicles_a_S = zeros(10,10,fig_N) ;
    for i=1:fig_N
        a_aux_0          = M_P.a_mat[1:fig_sample[i],1]   ;
        a_aux_T          = M_P.a_mat[1:fig_sample[i],end] ;
        deciles_a_S[:,i] = percentile( a_aux_T , 0:10:100 )   ; # Deciles based on end of sample (hoping for stationariety)
        for p=1:10
            ind_d = findall(x-> deciles_a_S[p]<=x<=deciles_a_S[p+1], a_aux_0 ) ; # Find "i" in each decile in t=1
            a_T = a_aux_T[ind_d] ; # Follow them to t=10
            Tr_decicles_a_S[p,:,i] = [100*sum( deciles_a_S[pT].<=a_T.<=deciles_a_S[pT+1] )/length(ind_d) for pT in 1:10] ; # Get transition rates 
        end 
    end 

## Figure with transitions of bottom decile
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( fig_sample/1000 , Tr_decicles_a_S[1,:,:]' ,label=nothing)
    hline!( 100*[Tr_deciles_a[1,:]] , c=:orange , w=2 , label=nothing )
    ylims!(0,ceil(maximum(Tr_decicles_a_S[1,:,:]/5))*5)
    xlims!(1,M_P.N_Panel/1000)
    title!("1st Decile Transitions",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/Tr_Decile_1_Simul.pdf")

## Figure with transitions of top decile
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( fig_sample/1000 , Tr_decicles_a_S[end,:,:]' ,label=nothing)
    hline!( 100*[Tr_deciles_a[end,:]] , c=:orange , w=2 , label=nothing )
    ylims!(0,ceil(maximum(Tr_decicles_a_S[end,:,:]/5))*5)
    xlims!(1,M_P.N_Panel/1000)
    title!("10th Decile Transitions",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/Tr_Decile_10_Simul.pdf")


###################################################################
###################################################################
## 5 Year autocorrelation of consumption first quintile
    av_c_q_S  = zeros(2,fig_N)    ; # Average Consumption in first quintile 
    sd_c_q_S  = zeros(2,fig_N)    ; # Std Dev of Consumption in first quintile 
    cor_c_q_S = zeros(fig_N)      ; # Cor of Consumption in first quintile 
    for i=1:fig_N
        # Fix current sample 
        a_aux_0       = M_P.a_mat[1:fig_sample[i],8]   ;
        c_aux_0       = M_P.c_mat[1:fig_sample[i],8]   ;
        c_aux_T       = M_P.c_mat[1:fig_sample[i],end] ;
        # Find index of first quintile 
        ind_q         = findall(x-> 0<=x<=deciles_a_S[3], a_aux_0 ) ; # Index of first quintile of assets in T-1
        # Compute moments 
        av_c_q_S[1,i] = mean( c_aux_0[ind_q] ) ; sd_c_q_S[1,i] = std( c_aux_0[ind_q] ) ;
        av_c_q_S[2,i] = mean( c_aux_T[ind_q] ) ; sd_c_q_S[2,i] = std( c_aux_T[ind_q] ) ;
        cor_c_q_S[i]  = cor( [c_aux_0[ind_q] c_aux_T[ind_q]]  ;dims=1)[2]              ;  
    end 
    
## Figure   
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( fig_sample/1000 , cor_c_q_S ,label=nothing)
    hline!( [cor_c_q] ,c=:orange  , w=3 , label=nothing )
    ylims!(0,1)
    xlims!(1,M_P.N_Panel/1000)
    title!("Auto-Correlation of Consumption",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/Cons_Corr_Simul.pdf")