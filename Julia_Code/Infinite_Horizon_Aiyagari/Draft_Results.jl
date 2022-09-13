###################################################################
###################################################################
###################################################################
## Simulate to get restuls for draft's graphs and tables



#=
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
            Γ_N = Histogram_Iteration(M,N,copy(Γ_0)) ;
            # Obtain transitions (sum over new distribution within deciles)
            Tr_deciles_a[i,:] = [100*sum( Γ_N[deciles_a[j]+1:deciles_a[j+1],:,:] ) for j in 1:10] ;
        end 

    # Return Moments 
        return Tr_deciles_a, deciles_mat
end 



###################################################################
## 3 Year autocorrelation of consumption 
function H_Moments_Auto_Correlation(M::Model,N::Int,a_min,a_max)

    # Get average consumption and standard deviation 
        # Initial distribution over asset range
        Γ_q                  = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ)   ; 
        Γ_q[a_min:a_max,:,:] = M.Γ[a_min:a_max,:,:]            ; 
        Γ_q                  = Γ_q/sum(Γ_q)                    ; 
            # println("Γ_q=$(sum(Γ_q)) - Γ_q_active=$(sum(Γ_q[a_min:a_max,:,:])) ")
        # Consumption 
        av_c_q = sum( M.G_c_fine[a_min:a_max,:,:].*Γ_q[a_min:a_max,:,:] ) ;
        sd_c_q = sqrt( sum( ((M.G_c_fine[a_min:a_max,:,:].-av_c_q).^2).*Γ_q[a_min:a_max,:,:] ) ) ;
        # Assets 
        av_a_q = sum( M.a_mat_fine[a_min:a_max,:,:].*Γ_q[a_min:a_max,:,:] ) ;
        sd_a_q = sqrt( sum( ((M.a_mat_fine[a_min:a_max,:,:].-av_a_q).^2).*Γ_q[a_min:a_max,:,:] ) ) ;
        # ϵ
        av_ϵ_q = sum( log.(M.ϵ_mat_fine[a_min:a_max,:,:]).*Γ_q[a_min:a_max,:,:] ) ;
        sd_ϵ_q = sqrt( sum( ((log.(M.ϵ_mat_fine[a_min:a_max,:,:]).-av_ϵ_q).^2).*Γ_q[a_min:a_max,:,:] ) ) ;
        # ζ
        av_ζ_q = sum( log.(M.ζ_mat_fine[a_min:a_max,:,:]).*Γ_q[a_min:a_max,:,:] ) ;
        sd_ζ_q = sqrt( sum( ((log.(M.ζ_mat_fine[a_min:a_max,:,:]).-av_ζ_q).^2).*Γ_q[a_min:a_max,:,:] ) ) ;
    # Future average and standard deviation conditional on first quintile
        # Iterate distribution 
        Γ_qN = Histogram_Iteration(M,N,copy(Γ_q)) ;
            # println("Γ_qN=$(sum(Γ_qN)) - Γ_qN_active=$(sum(Γ_qN[a_min:a_max,:,:])) ")
        # Consumption
        av_c_q_N = sum( M.G_c_fine.*Γ_qN  )                            ;
        sd_c_q_N = sqrt( sum( ((M.G_c_fine.-av_c_q_N).^2).*Γ_qN  ) )   ; 
        # Assets 
        av_a_q_N = sum( M.a_mat_fine.*Γ_qN  )                          ;
        sd_a_q_N = sqrt( sum( ((M.a_mat_fine.-av_a_q_N).^2).*Γ_qN  ) ) ; 
        # ϵ 
        av_ϵ_q_N = sum( log.(M.ϵ_mat_fine).*Γ_qN  )                          ;
        sd_ϵ_q_N = sqrt( sum( ((log.(M.ϵ_mat_fine).-av_ϵ_q_N).^2).*Γ_qN  ) ) ; 
        # ζ 
        av_ζ_q_N = sum( log.(M.ζ_mat_fine).*Γ_qN  )                          ;
        sd_ζ_q_N = sqrt( sum( ((log.(M.ζ_mat_fine).-av_ζ_q_N).^2).*Γ_qN  ) ) ; 

    # Compute integrand of correlation 
        cov_c = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
        cov_a = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
        cov_ϵ = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
        cov_ζ = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
        for i_ζ=1:M.n_ζ # Current ζ
        for i_ϵ=1:M.n_ϵ # Current ϵ
            for i_a=a_min:a_max # Current a
                # For each state in state space (a,ϵ,ζ) get conditional distribution 
                Γ_0 = zeros(M.n_a_fine,M.n_ϵ,M.n_ζ) ; 
                Γ_0[i_a,i_ϵ,i_ζ] = 1 ; 
                # Iterate distribution
                Γ_N = Histogram_Iteration(M,N,copy(Γ_0)) ;
                # Get portion of integrand (C_0 - av_c)*(C_n - av_c)*Γ_n
                cov_c[i_a,i_ϵ,i_ζ] = sum( (      M.G_c_fine[i_a,i_ϵ,i_ζ]   .-av_c_q )*(      M.G_c_fine   .-av_c_q_N ).*Γ_N ) ;
                cov_a[i_a,i_ϵ,i_ζ] = sum( (      M.a_mat_fine[i_a,i_ϵ,i_ζ] .-av_a_q )*(      M.a_mat_fine .-av_a_q_N ).*Γ_N ) ;
                cov_ϵ[i_a,i_ϵ,i_ζ] = sum( ( log.(M.ϵ_mat_fine[i_a,i_ϵ,i_ζ]).-av_ϵ_q )*( log.(M.ϵ_mat_fine).-av_ϵ_q_N ).*Γ_N ) ;
                cov_ζ[i_a,i_ϵ,i_ζ] = sum( ( log.(M.ζ_mat_fine[i_a,i_ϵ,i_ζ]).-av_ζ_q )*( log.(M.ζ_mat_fine).-av_ζ_q_N ).*Γ_N ) ;
            end 
        end 
        end 

    # Integrate with respect to initial distribution
        cov_c_q  = sum( cov_c.*Γ_q )  ;     cor_c_q  = 100*cov_c_q/(sd_c_q*sd_c_q_N) ; 
        cov_a_q  = sum( cov_a.*Γ_q )  ;     cor_a_q  = 100*cov_a_q/(sd_a_q*sd_a_q_N) ; 
        cov_ϵ_q  = sum( cov_ϵ.*Γ_q )  ;     cor_ϵ_q  = 100*cov_ϵ_q/(sd_ϵ_q*sd_ϵ_q_N) ; 
        cov_ζ_q  = sum( cov_ζ.*Γ_q )  ;     cor_ζ_q  = 100*cov_ζ_q/(sd_ζ_q*sd_ζ_q_N) ; 
        # println("Γ_q=$(sum(Γ_q)) - Γ_q_active=$(sum(Γ_q[a_min:a_max,:,:])) ")
        
    # Print other moments 
        # println("Moments for autocorrelation of consumption:")
        # println("cor_c=$(round(cor_c_q,digits=2)) - av_c_0=$(round(av_c_q,digits=3)) - av_c_T=$(round(av_c_q_N,digits=3)) - sd_c_0=$(round(sd_c_q,digits=2)) - sd_c_T=$(round(sd_c_q_N,digits=2)) - cov_0N=$(round(cov_c_q,digits=3))")        
        # println("cor_a=$(round(cor_a_q,digits=2)) - av_a_0=$(round(av_a_q,digits=3)) - av_a_T=$(round(av_a_q_N,digits=3)) - sd_a_0=$(round(sd_a_q,digits=2)) - sd_a_T=$(round(sd_a_q_N,digits=2)) - cov_0N=$(round(cov_a_q,digits=3))")        
        # println("cor_ϵ=$(round(cor_ϵ_q,digits=2)) - av_ϵ_0=$(round(av_ϵ_q,digits=3)) - av_ϵ_T=$(round(av_ϵ_q_N,digits=3)) - sd_ϵ_0=$(round(sd_ϵ_q,digits=2)) - sd_ϵ_T=$(round(sd_ϵ_q_N,digits=2)) - cov_0N=$(round(cov_ϵ_q,digits=3))")        
        # println("cor_ζ=$(round(cor_ζ_q,digits=2)) - av_ζ_0=$(round(av_ζ_q,digits=3)) - av_ζ_T=$(round(av_ζ_q_N,digits=3)) - sd_ζ_0=$(round(sd_ζ_q,digits=2)) - sd_ζ_T=$(round(sd_ζ_q_N,digits=2)) - cov_0N=$(round(cov_ζ_q,digits=3))")        
        println("cor_c=$(round(cor_c_q,digits=2)) - cor_a=$(round(cor_a_q,digits=2)) - cor_ϵ=$(round(cor_ϵ_q,digits=2)) - cor_ζ=$(round(cor_ζ_q,digits=2)) - a∈[$(round(M.a_mat_fine[a_min],digits=2)),$(round(M.a_mat_fine[a_max],digits=2))]")        
    # Return Moment 
        return cor_c_q, cor_a_q, cor_ϵ_q, cor_ζ_q 
end 


###################################################################
## Run Histogram Simulation for Different Grids 
    H_grid_size = [250 250 500 1000 5000] ; 
    n_H = length(H_grid_size) ;
    H_Γ_timed = zeros(n_H)    ; 
    H_Γ_bytes = zeros(n_H)    ;
    H_M_timed = zeros(n_H,3)  ; 
    H_M_bytes = zeros(n_H,3)  ;

    H_a_grid  = zeros(H_grid_size[end],n_H) ; # Asset Grid
    H_Γ_a     = zeros(H_grid_size[end],n_H) ; # Distribution of assets

    H_Wealth_Stats = zeros(n_H,6) ; # 5 Percentiles (levels) + Average Assets
    H_Wealth_Share = zeros(n_H,5) ; # 5 Top Shares by percentiles
    H_Pareto_Coeff = zeros(n_H  ) ; # Pareto Coefficient
    pct_list = [90;95;99;99.9;99.99] ; # Percentiles to be computed

    H_Decile       = zeros(11,2 ,n_H) ;
    H_Decile_Tr    = zeros(10,10,n_H) ;
    N_Decile_Tr    = 9                ;

    H_Cons_Corr    = zeros(n_H)       ;
    H_A_Corr       = zeros(n_H)       ;
    H_ϵ_Corr       = zeros(n_H)       ;
    H_ζ_Corr       = zeros(n_H)       ;
    N_Cons_Corr    = 2                ;

    for i=1:n_H 

        println(" ")
        println("   Histogram with $(H_grid_size[i]) Grid Points")
        println(" ")

        # Set up model structure 
        M_Hist = Model(n_a_fine=H_grid_size[i],method=1,read_flag=true) ;
        # Solve for stationary distribution and save time and allocation - Adjust grid size
        M_Hist, H_Γ_timed[i], H_Γ_bytes[i] = @timed Aiyagari_Equilibrium(M_Hist);
        
        ## Grid and Distribution
        H_a_grid[1:H_grid_size[i],i] = M_Hist.a_grid_fine ;
        H_Γ_a[1:H_grid_size[i],i]    = dropdims( sum( M_Hist.Γ , dims=(3,2) ) , dims=(3,2) ) ;
       
        ## Moments 
        
        # 1-2) Top Wealth Shares and Pareto Coefficient 
        out, time, memory = @timed H_Moments_Top_Shares(M_Hist,pct_list) ;
        H_Wealth_Stats[i,:] = out[1] ; H_Wealth_Share[i,:] = out[2] ; H_Pareto_Coeff[i] = out[3] ; 
        H_M_timed[i,1] = time ; H_M_bytes[i,1] = memory ; 


        # 3) Decile Transitions 
        out, time, memory = @timed H_Moments_Decile_Transitions(M_Hist,N_Decile_Tr) ;
        H_Decile_Tr[:,:,i] = out[1] ; H_Decile[:,:,i] = out[2] ; H_M_timed[i,2] = time ;  H_M_bytes[i,2] = memory ; 
        
        
        # 4) Consumption Autocorrelation for first quintile
        out, time, memory = @timed H_Moments_Auto_Correlation(M_Hist,N_Cons_Corr,1,Int(H_Decile[3,1,i])) ;
        H_Cons_Corr[i] = out[1] ; H_A_Corr[i] = out[2] ; H_ϵ_Corr[i] = out[3] ; H_ζ_Corr[i] = out[4] ;  
        H_M_timed[i,3] = time ; H_M_bytes[i,3] = memory ;

    end 

    open(Hist_Folder*"/H_G_timed.csv", "w") do io
    writedlm(io, H_Γ_timed , ',')
    end;

    open(Hist_Folder*"/H_G_bytes.csv", "w") do io
    writedlm(io, H_Γ_bytes , ',')
    end;

    open(Hist_Folder*"/H_a_grid.csv", "w") do io
    writedlm(io, H_a_grid , ',')
    end;

    open(Hist_Folder*"/H_G_a.csv", "w") do io
    writedlm(io, H_Γ_a , ',')
    end;

    open(Hist_Folder*"/H_M_timed.csv", "w") do io
    writedlm(io, H_M_timed , ',')
    end;

    open(Hist_Folder*"/H_M_bytes.csv", "w") do io
    writedlm(io, H_M_bytes , ',')
    end;

    open(Hist_Folder*"/H_Wealth_Stats.csv", "w") do io
    writedlm(io, H_Wealth_Stats , ',')
    end;

    open(Hist_Folder*"/H_Wealth_Share.csv", "w") do io
    writedlm(io, H_Wealth_Share , ',')
    end;
    
    open(Hist_Folder*"/H_Pareto_Coeff.csv", "w") do io
    writedlm(io, H_Pareto_Coeff , ',')
    end;

    open(Hist_Folder*"/H_Decile.csv", "w") do io
    writedlm(io, H_Decile , ',')
    end;

    open(Hist_Folder*"/H_Decile_Tr.csv", "w") do io
    writedlm(io, H_Decile_Tr , ',')
    end;
 
    open(Hist_Folder*"/H_Cons_Corr.csv", "w") do io
    writedlm(io, H_Cons_Corr, ',')
    end;

    open(Hist_Folder*"/H_A_Corr.csv", "w") do io
    writedlm(io, H_A_Corr, ',')
    end;

    open(Hist_Folder*"/H_eps_Corr.csv", "w") do io
    writedlm(io, H_ϵ_Corr, ',')
    end;

    open(Hist_Folder*"/H_z_Corr.csv", "w") do io
    writedlm(io, H_ζ_Corr, ',')
    end;
    
                   
=#


#=
###################################################################
###################################################################
## Simulation 




###################################################################
## Run Simulation for Different Panel Size 

    # Set up model structures 
    M_Simul = Model(method=1,read_flag=true) ;
    # M_Panel = Model_Panel(N_Panel=1000000)   ; 

    # Set up discrete observations 
    S_sample = [50000 250000 500000 1000000 10000000]    ; 
    N_S      = length(S_sample)         ;
    pct_list = [90;95;99;99.9;99.99]    ;  

    S_M_timed      = zeros(N_S,4)       ; # 1-> Simulation 2->Top Shares 3->Decile Transition 4->Auto-correlation
    S_M_bytes      = zeros(N_S,4)       ; # 1-> Simulation 2->Top Shares 3->Decile Transition 4->Auto-correlation
    
    S_Wealth_Sample= zeros(N_S,S_sample[N_S]) ;
    S_Wealth_Stats = zeros(N_S,6)       ; 
    S_Wealth_Share = zeros(N_S,5)       ; 
    S_Pareto_Coeff = zeros(N_S  )       ;  

    S_Decile       = zeros(11,N_S)      ;
    S_Decile_Tr    = zeros(10,10,N_S)   ;

    S_Cons_Corr    = zeros(N_S)         ;
    S_A_Corr       = zeros(N_S)         ;
    S_ϵ_Corr       = zeros(N_S)         ;
    S_ζ_Corr       = zeros(N_S)         ;
    
    # Solve model 
    M_Simul, S_Γ_timed, S_Γ_bytes = @timed Aiyagari_Equilibrium(M_Simul);

    # # Simulate Panel 
    # M_Panel = Simulate_Panel_Dynasty(M_Simul,M_Panel) ; 


=#    
    ## Moments 

    for i=N_S # 1:N_S

        println(" ")
        println("Simulation with $(S_sample[i]) agents")
        println(" ")

    # Simulate Panel      
        M_Panel = Model_Panel(N_Panel=S_sample[i])   ;   
        M_Panel, S_M_timed[i,1], S_M_bytes[i,1] = @timed Simulate_Panel(M_Simul,M_Panel,false) ; 

    # 1-2) Top Wealth Shares and Pareto Coefficient 
        a, S_M_timed[i,2], S_M_bytes[i,2] = @timed begin 

        # Select sample 
        a_sample  = M_Panel.a_mat[:,end] ;# M_Panel.a_mat[1:S_sample[i],end] ;
        S_Wealth_Sample[i,1:S_sample[i]] = a_sample  ;

        # Average wealth 
        S_Wealth_Stats[i,end]     = mean( a_sample )   ;

        # Percentiles 
        S_Wealth_Stats[i,1:end-1] = percentile( a_sample , pct_list ) ;

        # Top Shares 
        S_Wealth_Share[i,:] = [ 100*sum( a_sample[ a_sample.>=S_Wealth_Stats[i,p] ]  )/(S_sample[i]*S_Wealth_Stats[i,end])  for p in 1:5] ;

        # Pareto Coefficient 
        Pareto_sample = sort( a_sample[ a_sample.>=1000 ] ) ; # Select and sort observations above $1M  
        Pareto_CCDF   = collect(length(Pareto_sample):-1:1)./length(Pareto_sample) ; # Counter CDF = 1- CDF
        S_Pareto_Coeff[i] = (log.(Pareto_sample[1:end-1]./1000)'*log.(Pareto_sample[1:end-1]./1000))\log.(Pareto_sample[1:end-1]./1000)'*log.(Pareto_CCDF[1:end-1]) ;

        # Time it 
        end


    # 3) Decile Transitions 
        a, S_M_timed[i,3], S_M_bytes[i,3] = @timed begin 

        a_aux_0       = M_Panel.a_mat[:,1]   ;
        a_aux_T       = M_Panel.a_mat[:,end] ;
        S_Decile[:,i] = percentile( a_aux_T , 0:10:100 ) ; # Deciles based on end of sample (hoping for stationariety)
        for p=1:10
            ind_d = findall(x-> S_Decile[p,i]<=x<=S_Decile[p+1,i], a_aux_0 ) ; # Find "i" in each decile in t=1
            a_T   = a_aux_T[ind_d] ; # Follow them to t=10
            S_Decile_Tr[p,:,i] = [100*sum( S_Decile[pT,i].<=a_T.<=S_Decile[pT+1,i] )/length(ind_d) for pT in 1:10] ; # Get transition rates 
        end 

        # Time it 
        end
    
    
    # 4) Consumption Autocorrelation for first quintile
        a, S_M_timed[i,4], S_M_bytes[i,4] = @timed begin 

        # Fix current sample and future sample 
        a_aux_0 = M_Panel.a_mat[:, 8 ] ; c_aux_0 = M_Panel.c_mat[:, 8 ] ; ϵ_aux_0 = M_Simul.ϵ_grid[M_Panel.ϵ_mat[:, 8 ]] ; ζ_aux_0 = M_Simul.ζ_grid[M_Panel.ζ_mat[:, 8 ]]   ;
        a_aux_T = M_Panel.a_mat[:,end] ; c_aux_T = M_Panel.c_mat[:,end] ; ϵ_aux_T = M_Simul.ϵ_grid[M_Panel.ϵ_mat[:,end]] ; ζ_aux_T = M_Simul.ζ_grid[M_Panel.ζ_mat[:,end]] ;
        # Find index of first quintile 
        ind_q         = findall(x-> 0<=x<=S_Decile[3,i], a_aux_0 ) ; # Index of first quintile of assets in T-1
        # Compute moments 
        S_Cons_Corr[i]= cor( [     c_aux_0[ind_q]       c_aux_T[ind_q] ]  ;dims=1)[2]  ;
        S_A_Corr[i]   = cor( [     a_aux_0[ind_q]       a_aux_T[ind_q] ]  ;dims=1)[2]  ;
        S_ϵ_Corr[i]   = cor( [log.(ϵ_aux_0[ind_q]) log.(ϵ_aux_T[ind_q])]  ;dims=1)[2]  ;
        S_ζ_Corr[i]   = cor( [log.(ζ_aux_0[ind_q]) log.(ζ_aux_T[ind_q])]  ;dims=1)[2]  ; 
            # Deconstruct measure: mean, var, cov 
            av_c_0   = mean(     c_aux_0[ind_q] ) ; av_c_T   = mean(     c_aux_T[ind_q] ) ; cov_c_0T = cov([c_aux_0[ind_q] c_aux_T[ind_q]]) ;  
            av_a_0   = mean(     a_aux_0[ind_q] ) ; av_a_T   = mean(     a_aux_T[ind_q] ) ; cov_a_0T = cov([a_aux_0[ind_q] a_aux_T[ind_q]]) ;  
            av_ϵ_0   = mean(log.(ϵ_aux_0[ind_q])) ; av_ϵ_T   = mean(log.(ϵ_aux_T[ind_q])) ; cov_ϵ_0T = cov([log.(ϵ_aux_0[ind_q]) log.(ϵ_aux_T[ind_q])]) ;  
            av_ζ_0   = mean(log.(ζ_aux_0[ind_q])) ; av_ζ_T   = mean(log.(ζ_aux_T[ind_q])) ; cov_ζ_0T = cov([log.(ζ_aux_0[ind_q]) log.(ζ_aux_T[ind_q])]) ;  
            println("Moments for autocorrelation of consumption:")
            println("cor_c=$(round(100*S_Cons_Corr[i],digits=2)) - av_c_0=$(round(av_c_0,digits=3)) - av_c_T=$(round(av_c_T,digits=3)) - sd_c_0=$(round(sqrt(cov_c_0T[1,1]),digits=2)) - sd_c_T=$(round(sqrt(cov_c_0T[2,2]),digits=2)) - cov_0N=$(round(cov_c_0T[1,2],digits=3))") 
            println("cor_a=$(round(100*S_A_Corr[i]   ,digits=2)) - av_a_0=$(round(av_a_0,digits=3)) - av_a_T=$(round(av_a_T,digits=3)) - sd_a_0=$(round(sqrt(cov_a_0T[1,1]),digits=2)) - sd_a_T=$(round(sqrt(cov_a_0T[2,2]),digits=2)) - cov_0N=$(round(cov_a_0T[1,2],digits=3))") 
            println("cor_ϵ=$(round(100*S_ϵ_Corr[i]   ,digits=2)) - av_ϵ_0=$(round(av_ϵ_0,digits=3)) - av_ϵ_T=$(round(av_ϵ_T,digits=3)) - sd_ϵ_0=$(round(sqrt(cov_ϵ_0T[1,1]),digits=2)) - sd_ϵ_T=$(round(sqrt(cov_ϵ_0T[2,2]),digits=2)) - cov_0N=$(round(cov_ϵ_0T[1,2],digits=3))") 
            println("cor_ζ=$(round(100*S_ζ_Corr[i]   ,digits=2)) - av_ζ_0=$(round(av_ζ_0,digits=3)) - av_ζ_T=$(round(av_ζ_T,digits=3)) - sd_ζ_0=$(round(sqrt(cov_ζ_0T[1,1]),digits=2)) - sd_ζ_T=$(round(sqrt(cov_ζ_0T[2,2]),digits=2)) - cov_0N=$(round(cov_ζ_0T[1,2],digits=3))") 
        # Time it 
        end  

    end

    open(MC_Folder*"/S_M_timed.csv", "w") do io
    writedlm(io, S_M_timed , ',')
    end;

    open(MC_Folder*"/S_M_bytes.csv", "w") do io
    writedlm(io, S_M_bytes , ',')
    end;


    open(MC_Folder*"/S_G_timed.csv", "w") do io
    writedlm(io, S_Γ_timed , ',')
    end;

    open(MC_Folder*"/S_G_bytes.csv", "w") do io
    writedlm(io, S_Γ_bytes , ',')
    end;


    S_Wealth_1 = S_Wealth_Sample[1:N_S-1,1:S_sample[N_S-1]] ; 
    S_Wealth_2 = S_Wealth_Sample[N_S,:] ; 
    open(MC_Folder*"/S_Wealth_Sample.csv", "w") do io
        writedlm(io, S_Wealth_1 , ',')
        end;
    open(MC_Folder*"/S_Wealth_Sample_Large.csv", "w") do io
        writedlm(io, round.(S_Wealth_2; digits=3) , ',')
        end;


    open(MC_Folder*"/S_Wealth_Stats.csv", "w") do io
    writedlm(io, S_Wealth_Stats , ',')
    end;

    open(MC_Folder*"/S_Wealth_Share.csv", "w") do io
    writedlm(io, S_Wealth_Share , ',')
    end;
    
    open(MC_Folder*"/S_Pareto_Coeff.csv", "w") do io
    writedlm(io, S_Pareto_Coeff , ',')
    end;

    open(MC_Folder*"/S_Decile.csv", "w") do io
    writedlm(io, S_Decile , ',')
    end;

    open(MC_Folder*"/S_Decile_Tr.csv", "w") do io
    writedlm(io, S_Decile_Tr , ',')
    end;
    
    open(MC_Folder*"/S_Cons_Corr.csv", "w") do io
    writedlm(io, S_Cons_Corr, ',')
    end;
        
    open(MC_Folder*"/S_A_Corr.csv", "w") do io
    writedlm(io, S_A_Corr, ',')
    end;

    open(MC_Folder*"/S_eps_Corr.csv", "w") do io
    writedlm(io, S_ϵ_Corr, ',')
    end;

    open(MC_Folder*"/S_z_Corr.csv", "w") do io
    writedlm(io, S_ζ_Corr, ',')
    end;
