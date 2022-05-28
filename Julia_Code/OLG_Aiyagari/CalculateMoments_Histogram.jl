###################################################################
###################################################################
###################################################################
## Compute moments using the histogram iteration method 



###################################################################
###################################################################
## Cohort evolution of 90-10 percentiles of wealth 
println("\n===============================================")
println("90-10 Percentiles for Cohorts")
## Newborns 
    pct_C_age_1 = zeros(2,M_Aiyagari.p.Max_Age) ; 
    for h=1:M_Aiyagari.p.Max_Age
        Γ_aux  = sum(M_Aiyagari.Γ[:,:,h],dims=(2,3))[:,1]   ; 
        CDF_h  = cumsum(Γ_aux/sum(Γ_aux))                   ; 
        pct_C_age_1[1,h]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_h).>=(10)][1] ] ;
        pct_C_age_1[2,h]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_h).>=(90)][1] ] ;
    end 

## 45-year olds with more than median income
    age_0 = 26 ; 
    pct_C_age_2 = zeros(2,M_Aiyagari.p.Max_Age-(age_0-1)) ;  
    # Initial value 
    Γ_h = zeros(size(M_Aiyagari.Γ))                                  ;
    Γ_h[:,med_ϵ:end,age_0]  = M_Aiyagari.Γ[:,med_ϵ:end,age_0]/sum(M_Aiyagari.Γ[:,med_ϵ:end,age_0]) ; 
    CDF_1  = cumsum( sum( Γ_h , dims=(2,3) )[:,1] ) ;
    pct_C_age_2[1,1]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_1).>=(10)][1] ] ;
    pct_C_age_2[2,1]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_1).>=(90)][1] ] ;
    # Simulate cohort forward 
    for h=2:M_Aiyagari.p.Max_Age-(age_0-1)
    global Γ_h = Histogram_Iteration(M_Aiyagari,1,Γ_h) ; 
    global Γ_h[:,:, 1:end .!= (h+(age_0-1)) ] .= 0  ; 
    global Γ_h .= Γ_h/sum(Γ_h) ; 
    CDF_h  = cumsum( sum( Γ_h , dims=(2,3) )[:,1] ) ;
    pct_C_age_2[1,h]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_h).>=(10)][1] ] ;
    pct_C_age_2[2,h]  = M_Aiyagari.a_grid_fine[ collect(1:M_Aiyagari.n_a_fine)[(100*CDF_h).>=(90)][1] ] ;
    end 
    
println("===============================================\n")


###################################################################
###################################################################
## Auto-correlation of wealth Ages 35 and 65 (conditional on survival)
println("\n===============================================")
println("Autocorr of Wealth: ages 35-65")
    age_0 = 16        ; # Start agents at age 35 
    n_H   = 65-(35-1) ; # Simulate for n_H years 
    # Get a_grid vector
    a_grid_vec = collect(M_Aiyagari.a_grid_fine) ;
    # Turn off death to ensure balance panel 
    M_Aiyagari_C = Model(M_Aiyagari; p=Par(M_Aiyagari.p; Surv_Pr = [ones(M_Aiyagari.p.Max_Age-1);0]) ) ; 
    # Initial distribution 
    Γ_h = zeros(size(M_Aiyagari.Γ))                                        ;
    Γ_h[:,:,age_0]  = M_Aiyagari.Γ[:,:,age_0]/sum(M_Aiyagari.Γ[:,:,age_0]) ; 
    # Average and standard deviation at initial distribution 
    av_a_0 =  sum( a_grid_vec.*sum(Γ_h,dims=(2,3)) ) ;
    sd_a_0 =  sqrt( sum( ((a_grid_vec.-av_a_0).^2).*sum(Γ_h,dims=(2,3)) ) );
    av_a_N =  sum( a_grid_vec.*sum(M_Aiyagari.Γ[:,:,n_H+age_0],dims=(2,3)) )/sum(M_Aiyagari.Γ[:,:,n_H+age_0]) ;
    sd_a_N =  sqrt( sum( ((a_grid_vec.-av_a_N).^2).*sum(M_Aiyagari.Γ[:,:,n_H+age_0],dims=(2,3))/sum(M_Aiyagari.Γ[:,:,n_H+age_0]) ) );
    # Follow each initial state and fill in integrand 
    cov_a_3565 = zeros(M_Aiyagari.n_a_fine,M_Aiyagari.n_ϵ) ;
    for i_ϵ=1:M_Aiyagari.n_ϵ
        println(" Iterating with i_ϵ=$i_ϵ")
    for i_a=1:M_Aiyagari.n_a_fine 
        # println(" Iterating with i_ϵ=$i_ϵ , i_a=$i_a")
        # For each state in state space (a,ϵ,ζ) get conditional distribution 
        Γ_0 = zeros(M_Aiyagari.n_a_fine,M_Aiyagari.n_ϵ,M_Aiyagari.p.Max_Age) ; 
        Γ_0[i_a,i_ϵ,age_0] = 1 ; 
        # Iterate distribution
        Γ_N = Histogram_Iteration(M_Aiyagari_C,n_H,Γ_0) ;
        # Fill in integrand (a_0 - av_a)*(a_n - av_a)*Γ_N
        cov_a_3565[i_a,i_ϵ] = sum( ( a_grid_vec[i_a].-av_a_0 )*( a_grid_vec.-av_a_N ).*sum(Γ_N,dims=(2,3))/sum(Γ_N) ) ;
    end 
    end 
    # Integrate covariance 
    cor_a_3565 = sum( cov_a_3565.*Γ_h[:,:,age_0]  )/sqrt(sd_a_0^2*sd_a_N^2) ;


println("===============================================\n")
