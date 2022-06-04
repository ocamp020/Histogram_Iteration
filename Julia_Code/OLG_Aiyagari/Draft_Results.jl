###################################################################
###################################################################
###################################################################
## Graphs and Tables for Draft 




###################################################################
###################################################################
## Histogram 



###################################################################
## Wealth Profile
function H_Moments_Wealth_Profile(M::Model,pct_list,age_0::Int)
    
    # Define marginal distributions 
    Γ_a   = dropdims( sum( M.Γ , dims=(3,2) ) , dims=(3,2) ) ; # Assets 
    Γ_ϵ   = dropdims( sum( M.Γ , dims=(1,3) ) , dims=(1,3) ) ; # Labor Efficiency 
    Γ_age = dropdims( sum( M.Γ , dims=(1,2) ) , dims=(1,2) ) ; # Returns 

    # Define index for median shocks
    med_ϵ   = convert(Int64,round(M.n_ϵ/2)) ;

    # Unconditioned Wealth Profile 
    pct_a_age_nb   = zeros(M.p.Max_Age,length(pct_list)+1) ; 
    for h=1:M.p.Max_Age
        Γ_aux  = sum(M.Γ[:,:,h],dims=(2,3))[:,1]   ; 
        CDF_h  = cumsum(Γ_aux/sum(Γ_aux))                   ; 
        pct_a_age_nb[h,1:end-1]  = [ M.a_grid_fine[ collect(1:M.n_a_fine)[(100*CDF_h).>=(pct_list[p])][1] ]  for p in 1:length(pct_list) ];
        pct_a_age_nb[h,end]      = sum(M.a_mat_fine[:,:,h].*M.Γ[:,:,h])/Γ_age[h]
    end         
    
    # 45-year olds with more than median income
    pct_a_age_45 = zeros(M.p.Max_Age-(age_0-1),length(pct_list)+1) ;  
    # Initial value 
    Γ_h = zeros(size(M.Γ))                                  ;
    Γ_h[:,med_ϵ:end,age_0]  = M.Γ[:,med_ϵ:end,age_0]/sum(M.Γ[:,med_ϵ:end,age_0]) ; 
    CDF_1 = cumsum( sum( Γ_h , dims=(2,3) )[:,1] ) ;
    pct_a_age_45[1,1:end-1] = [ M.a_grid_fine[ collect(1:M.n_a_fine)[(100*CDF_1).>=(pct_list[p])][1] ] for p in 1:length(pct_list) ];
    pct_a_age_45[1,end]     = sum(M.a_mat_fine[:,:,age_0].*Γ_h[:,:,age_0]) ; 
    # Simulate cohort forward 
    for h=2:M.p.Max_Age-(age_0-1)
    Γ_h = Histogram_Iteration(M,1,Γ_h) ; 
    Γ_h[:,:, 1:end .!= (h+(age_0-1)) ] .= 0  ; 
    Γ_h .= Γ_h/sum(Γ_h) ; 
    CDF_h  = cumsum( sum( Γ_h , dims=(2,3) )[:,1] ) ;
    pct_a_age_45[h,1:end-1] = [ M.a_grid_fine[ collect(1:M.n_a_fine)[(100*CDF_h).>=(pct_list[p])][1] ] for p in 1:length(pct_list) ];
    pct_a_age_45[h,end]     = sum(M.a_mat_fine[:,:,h].*Γ_h[:,:,h]) ; 
    end 
    
    # Return Moments 
    return pct_a_age_nb, pct_a_age_45 
end 


###################################################################
## Wealth Autocorrelation 35-55
function H_Moments_Wealth_Corr(M::Model,age_0::Int,age_T::Int)
    n_H   = age_T-age_0 ; # Simulate for n_H years 
    # Get a_grid vector
    a_grid_vec = collect(M.a_grid_fine) ;
    # Turn off death to ensure balance panel 
    M_C = Model(M; p=Par(M.p; Surv_Pr = [ones(M.p.Max_Age-1);0]) ) ; 
    # Initial distribution 
    Γ_h = zeros(size(M.Γ))                               ;
    Γ_h[:,:,age_0]  = M.Γ[:,:,age_0]/sum(M.Γ[:,:,age_0]) ; 
    # Average and standard deviation at initial distribution 
    av_a_0 =  sum( a_grid_vec.*sum(Γ_h,dims=(2,3)) ) ;
    sd_a_0 =  sqrt( sum( ((a_grid_vec.-av_a_0).^2).*sum(Γ_h,dims=(2,3)) ) );
    av_a_N =  sum( a_grid_vec.*sum(M.Γ[:,:,n_H+age_0],dims=(2,3)) )/sum(M.Γ[:,:,n_H+age_0]) ;
    sd_a_N =  sqrt( sum( ((a_grid_vec.-av_a_N).^2).*sum(M.Γ[:,:,n_H+age_0],dims=(2,3))/sum(M.Γ[:,:,n_H+age_0]) ) );
    # Follow each initial state and fill in integrand 
    cov_a = zeros(M.n_a_fine,M.n_ϵ) ;
    for i_ϵ=1:M.n_ϵ
        # println(" Iterating with i_ϵ=$i_ϵ")
    for i_a=1:M.n_a_fine 
        # println(" Iterating with i_ϵ=$i_ϵ , i_a=$i_a")
        # For each state in state space (a,ϵ,ζ) get conditional distribution 
        Γ_0 = zeros(M.n_a_fine,M.n_ϵ,M.p.Max_Age) ; 
        Γ_0[i_a,i_ϵ,age_0] = 1 ; 
        # Iterate distribution
        Γ_N = Histogram_Iteration(M_C,n_H,Γ_0) ;
        # Fill in integrand (a_0 - av_a)*(a_n - av_a)*Γ_N
        cov_a[i_a,i_ϵ] = sum( ( a_grid_vec[i_a].-av_a_0 )*( a_grid_vec.-av_a_N ).*sum(Γ_N,dims=(2,3))/sum(Γ_N) ) ;
    end 
    end 
    # Integrate covariance 
    cor_a = sum( cov_a.*Γ_h[:,:,age_0]  )/sqrt(sd_a_0^2*sd_a_N^2) ;

    # Return Moment
    return cor_a
end 




###################################################################
## Run Histogram Simulation for Different Grids 
    H_grid_size = 100:50:1000 ; 
    n_H = length(H_grid_size) ;
    H_Γ_timed = zeros(n_H)    ; 
    H_Γ_bytes = zeros(n_H)    ;
    H_M_timed = zeros(n_H,2)  ; 
    H_M_bytes = zeros(n_H,2)  ;

    pct_list = [90;95;99;99.9;99.99] ; 
    age_0_Wealth_Profile =  26       ; 
    age_0_Wealth_Corr    =  16       ; 
    age_T_Wealth_Corr    =  36       ; 

    H_Wealth_Profile_NB = zeros(M_Aiyagari.p.Max_Age,6,n_H) ; 
    H_Wealth_Profile_45 = zeros(M_Aiyagari.p.Max_Age-age_0_Wealth_Profile+1,6,n_H) ; 
    H_Wealth_Corr       = zeros(n_H)   ; 

    for i=1:n_H 

        println("   Histogram with $(H_grid_size[i]) Grid Points")

        # Set up model structure 
        M_Hist = Model(n_a_fine=H_grid_size[i],read_flag=false) ;
        # Solve for stationary distribution and save time and allocation - Adjust grid size
        M_Hist, H_Γ_timed[i], H_Γ_bytes[i] = @timed Aiyagari_Equilibrium(M_Hist);
        
        ## Moments 
        
        # 1) Wealth Profiles 
        out, time, memory = @timed H_Moments_Wealth_Profile(M_Hist,pct_list,age_0_Wealth_Profile) ;
        H_Wealth_Profile_NB[:,:,i] = out[1] ; H_Wealth_Profile_45[:,:,i] = out[2] ; H_M_timed[i,1] = time ; H_M_bytes[i,1] = memory ; 

        # 2) Wealth Autocorrelation 35-55
        out, time, memory = @timed H_Moments_Wealth_Corr(M_Hist,age_0_Wealth_Corr,age_T_Wealth_Corr) ;
        H_Wealth_Corr[i] = out ;  H_M_timed[i,2] = time ; H_M_bytes[i,2] = memory ;

    end 


    open(Hist_Folder*"/H_G_timed.csv", "w") do io
    writedlm(io, H_Γ_timed , ',')
    end;

    open(Hist_Folder*"/H_G_bytes.csv", "w") do io
    writedlm(io, H_Γ_bytes , ',')
    end;

    open(Hist_Folder*"/H_M_timed.csv", "w") do io
    writedlm(io, H_M_timed , ',')
    end;

    open(Hist_Folder*"/H_M_bytes.csv", "w") do io
    writedlm(io, H_M_bytes , ',')
    end;

    open(Hist_Folder*"/H_Wealth_Profile_NB.csv", "w") do io
    writedlm(io, H_Wealth_Profile_NB , ',')
    end;

    open(Hist_Folder*"/H_Wealth_Profile_45.csv", "w") do io
    writedlm(io, H_Wealth_Profile_45 , ',')
    end;
    
    open(Hist_Folder*"/H_Wealth_Corr.csv", "w") do io
    writedlm(io, H_Wealth_Corr , ',')
    end;

                    
    

###################################################################
###################################################################
## Simulation 




###################################################################
## Run Simulation for Different Panel Size 

    # Set up age limits 
    age_0_Wealth_Profile =  26           ; 
    age_0_Wealth_Corr    =  16           ; 
    age_T_Wealth_Corr    =  36           ; 

    # Set up model structures
    M_Simul = Model(read_flag=false) ;
    M_Panel = Model_Panel( N_Panel=500000) ;
    M_C_45  = Model_Cohort(N_Panel=500000,T_Panel=M_Aiyagari.p.Max_Age-(age_0_Wealth_Profile-1))    ; 
    M_C_35  = Model_Cohort(N_Panel=500000,T_Panel=M_Aiyagari.p.Max_Age-(age_0_Wealth_Corr-1)   )    ; 
        

    # Set up discrete observations 
    S_sample = 10000:1000:M_Panel.N_Panel ; 
    N_S      = length(S_sample)          ;
    pct_list = [90;95;99;99.9;99.99]     ; 
    med_ϵ  = convert(Int64,round(M.n_ϵ/2)) ;

    S_M_timed      = zeros(N_S,2)     ; 
    S_M_bytes      = zeros(N_S,2)     ;
    
    S_Wealth_Profile_NB = zeros(M_Aiyagari.p.Max_Age,6,N_S) ; 
    S_Wealth_Profile_45 = zeros(M_Aiyagari.p.Max_Age,6,N_S) ; 
    S_Wealth_Corr       = zeros(N_S)   ; 



    # Solve model 
    M_Simul, S_Γ_timed, S_Γ_bytes = @timed Aiyagari_Equilibrium(M_Simul);

    # Simulate Panel 
    M_Panel = Simulate_Panel_Timed(M_Simul,M_Panel) ; 

    # Simulate Cohorts
    M_C_45 = Simulate_Cohort_Timed(M_Simul,M_C_45,age_0_Wealth_Profile) ; 
    M_C_35 = Simulate_Cohort_Timed(M_Simul,M_C_35,age_0_Wealth_Corr   ) ; 
    
    ## Moments 

    # 1) Wealth Profiles 
    for i=1:N_S 
        a, S_M_timed[i,1], S_M_bytes[i,1] = @timed begin 

        ## New Borns            
        # Select sample 
        a_sample    = M_Panel.a_mat[1:S_sample[i],end]    ;
        h_sample    = M_Panel.h_mat[1:S_sample[i],end]    ;
        for age = 1:M_Simul.p.Max_Age
            
        if sum(h_sample.==age)>0
            S_Wealth_Profile_NB[age,1:end-1,i]  .= percentile( a_sample[h_sample.==age] , pct_list ) ;
            S_Wealth_Profile_NB[age,  end  ,i]   = mean(       a_sample[h_sample.==age]  ) ;
        else 
            println(" Error in sample at i=$i, sample=$(S_sample[i]) and age=$age: $(sum(h_sample.==age)) Observations")
            S_Wealth_Profile_NB[age,:,i]   .= NaN ;
        end
    
        end 
 
        ## 45-year olds with more than median income
        # Select sample 
        a_sample    = M_C_45.a_mat[1:S_sample[i],:]    ;
        h_sample    = M_C_45.h_mat[1:S_sample[i],:]    ;
        ϵ_sample    = M_C_45.ϵ_mat[1:S_sample[i],:]    ;
        ind_45      = ϵ_sample[:,1].>=med_ϵ            ;
        a_sample    = a_sample[ind_45,:]               ;
        h_sample    = h_sample[ind_45,:]               ;
        for age = 1:M_Aiyagari.p.Max_Age-(age_0_Wealth_Profile-1)
        
            ind_alive   = h_sample[:,age].==1          ; 

        if sum(ind_alive)>0 

            S_Wealth_Profile_45[age,1:end-1,i]  = percentile( a_sample[ind_alive,age] , pct_list ) ;
            S_Wealth_Profile_45[age,  end  ,i]  = mean( a_sample[ind_alive,age]  ) ;

        else 
            S_Wealth_Profile_45[age,:,i] .= NaN ; 
        end 
        
        end 
        
        # Time it 
        end
    end 


    # 2) Wealth Autocorrelation 35-65
    for i=1:N_S 
        a, S_M_timed[i,2], S_M_bytes[i,2] = @timed begin 

        # Fix current sample 
        a_aux_0       = M_C_35.a_mat[1:S_sample[i],1]   ;
        a_aux_T       = M_C_35.a_mat[1:S_sample[i],age_T_Wealth_Corr-age_0_Wealth_Corr] ;
        # Compute moments 
        S_Wealth_Corr[i]  = cor( [a_aux_0 a_aux_T]  ;dims=1)[2]              ;  
        
        # Time it 
        end
    end 
    

    open(MC_Folder*"/S_M_timed.csv", "w") do io
    writedlm(io, S_M_timed , ',')
    end;

    open(MC_Folder*"/S_M_bytes.csv", "w") do io
    writedlm(io, S_M_bytes , ',')
    end;

    open(MC_Folder*"/S_Wealth_Profile_NB.csv", "w") do io
    writedlm(io, S_Wealth_Profile_NB , ',')
    end;

    open(MC_Folder*"/S_Wealth_Profile_45.csv", "w") do io
    writedlm(io, S_Wealth_Profile_45 , ',')
    end;
    
    open(MC_Folder*"/S_Wealth_Corr.csv", "w") do io
    writedlm(io, S_Wealth_Corr , ',')
    end;

                        




###################################################################
###################################################################
## Graphs and Tables 
    

###################################################################
## Wealth Profiles
    # 10k Observations 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    # Average percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,end, 9] , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,10] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,end, 9] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,end,10] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 90th percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,3, 9] , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,3,10] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,3, 9] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,3,10] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 99th percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,4, 9] , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,4,10] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,4, 9] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,4,10] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    ylims!(0,ceil(maximum(S_Wealth_Profile_45[:,4]/10))*10)
    title!("Age Profile: Average, 90th, 99th Percentiles of Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_10k.pdf")



    # 100k Observations 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    # Average percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,end, 9]  , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,100] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,end, 9]  , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,end,100] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 90th percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,3, 9]  , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,3,100] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,3, 9]  , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,3,100] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 99th percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,4, 9]  , marker=(:circle ,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,4,100] , marker=(:diamond,2 ,:gray70),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[age_0_Wealth_Profile:end,4, 9]  , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[age_0_Wealth_Profile:end,4,100] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    ylims!(0,ceil(maximum(S_Wealth_Profile_45[:,4]/10))*10)
    title!("Age Profile: Average, 90th, 99th Percentiles of Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_100k.pdf")


###################################################################
## Wealth Auto-Correlation 
    Mat = [ H_Wealth_Corr[1]  H_Wealth_Corr[4]  H_Wealth_Corr[9]   S_Wealth_Corr[10]                 S_Wealth_Corr[100]                 S_Wealth_Corr[500] ;
            H_Γ_timed[1]      H_Γ_timed[4]      H_Γ_timed[9]       S_Γ_timed+M_Panel.t_vec[10000]    S_Γ_timed+M_Panel.t_vec[100000]    S_Γ_timed+M_Panel.t_vec[500000] ;
            H_Timed[1,2]      H_Timed[4,2]      H_Timed[9,2]       S_M_timed[10,2]                   S_M_timed[100,2]                   S_M_timed[500,2]   ] ;


 