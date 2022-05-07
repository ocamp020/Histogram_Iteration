###################################################################
###################################################################
###################################################################
## Compute moments using the histogram iteration method 



###################################################################
###################################################################
## 10 Year transition rates for decile 
println("\n===============================================")
println("Computing transition matrix for deciles")
# Set time horizon 
    n_H = 10 ;

# Find deciles 
    deciles_a    = Array{Int64}(undef,11) ;
    deciles_a[1] = 0         ;
    for i=2:10 
    deciles_a[i]  = collect(1:M_Aiyagari.n_a_fine)[(100*CDF_a).>=(10*(i-1))][1] ;
    end 
    deciles_a[11] = M_Aiyagari.n_a_fine ;

# For each decile compute transitions   
    Tr_deciles_a = zeros(10,10) ; 
    for i=1:10
        println("\n Decile $i")
        # Conditional distribution on current decile 
        Γ_0                                    = zeros(M_Aiyagari.n_a_fine,M_Aiyagari.n_ϵ,M_Aiyagari.n_ζ)   ; 
        Γ_0[deciles_a[i]+1:deciles_a[i+1],:,:] = M_Aiyagari.Γ[deciles_a[i]+1:deciles_a[i+1],:,:]            ; 
        Γ_0 = Γ_0/sum(Γ_0) ; 
        # Iterate distribution 
        Γ_N = Histogram_Iteration(M_Aiyagari,n_H,Γ_0) ;
        # Obtain transitions (sum over new distribution within deciles)
        Tr_deciles_a[i,:] = [sum( Γ_N[deciles_a[j]+1:deciles_a[j+1],:,:] ) for j in 1:10] ;
    end 

# Print transtion matrix
    println("   Wealth Deciles Transition Matrix (N=$n_H) ")
    for i=1:10
        println("   Decile $i, a<=\$$(round(M_Aiyagari.a_grid_fine[deciles_a[i+1]],digits=2))k // Transition:  $(round.(100*Tr_deciles_a[i,:],digits=2))")
    end 

# Save transition matrix to file 


println("===============================================\n")


###################################################################
###################################################################
## 5 Year autocorrelation of consumption 
println("\n===============================================")
println("Computing autocorrelation of consumption")
# Set time horizon 
    n_H = 2 ;

# Get average consumption and standard deviation 
    av_c = sum( M_Aiyagari.G_c_fine.*M_Aiyagari.Γ )                      ;
    sd_c = sqrt( sum( ((M_Aiyagari.G_c_fine.-av_c).^2).*M_Aiyagari.Γ ) ) ;

# Compute integrand of correlation   
    cov_c = zeros(M_Aiyagari.n_a_fine,M_Aiyagari.n_ϵ,M_Aiyagari.n_ζ) ; 
    for i_ζ=1:M_Aiyagari.n_ζ # Current ζ
    for i_ϵ=1:M_Aiyagari.n_ϵ # Current ϵ
    for i_a=1:M_Aiyagari.n_a_fine # Current a
        # For each state in state space (a,ϵ,ζ) get conditional distribution 
        Γ_0 = zeros(M_Aiyagari.n_a_fine,M_Aiyagari.n_ϵ,M_Aiyagari.n_ζ) ; 
        Γ_0[i_a,i_ϵ,i_ζ] = 1 ; 
        # Iterate distribution
        Γ_N = Histogram_Iteration(M_Aiyagari,n_H,Γ_0) ;
        # Get portion of integrand (C_0 - av_c)*(C_n - av_c)*Γ_n
        cov_c[i_a,i_ϵ,i_ζ] = sum( (M_Aiyagari.G_c_fine[i_a,i_ϵ,i_ζ]-av_c)*(M_Aiyagari.G_c_fine.-av_c).*Γ_N ) ;
    end 
    end 
    end 

# Integrate with respect to initial distribution 
    cov_c = sum( cov_c.*M_Aiyagari.Γ )  ;
    cor_c = cov_c/(sd_c^2)              ; 

# Print result
    println("   Average Consumption:  \$$(round(av_c,digits=2))k ")
    println("   St. Dv. Consumption:  \$$(round(sd_c,digits=2))k ")
    println("   $n_H years auto-corr:   $(round(cor_c,digits=2)) ")

# Save result to file 

println("===============================================\n")