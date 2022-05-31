###################################################################
###################################################################
###################################################################
## Compute moments using the Monte-Carlo Simulation method 


###################################################################
###################################################################

## Run Simulations - Economy 
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
    open(File_Folder*"/S_h_mat.csv", "w") do io
        writedlm(io, M_P.h_mat , ',')
    end;

    M_C_45 = Model_Cohort(T_Panel=M_Aiyagari.p.Max_Age-25) ;
    M_C_45 = Simulate_Cohort(M_Aiyagari,M_C_45,26)         ;

    M_C_35 = Model_Cohort(T_Panel=M_Aiyagari.p.Max_Age-15) ;
    M_C_35 = Simulate_Cohort(M_Aiyagari_C,M_C_35,16)       ;


###################################################################
###################################################################
## Cohort evolution of average and 90-10 percentiles of wealth 
println("\n===============================================")
println("Average and 90-10 Percentiles for Cohorts")
## Newborns 
    sample_vec = [10000, 50000, 100000, 500000]  ; 
    pct_list   = [10;50;90;99]                   ;
    av_a_S     = zeros(M_Aiyagari.p.Max_Age,4)   ;
    pct_S      = zeros(4,M_Aiyagari.p.Max_Age,4) ; 
    for i = 1:4
        a_sample    = M_P.a_mat[1:sample_vec[i],end]    ;
        h_sample    = M_P.h_mat[1:sample_vec[i],end]    ;
        for age = 1:M_Aiyagari.p.Max_Age
        if sum(h_sample.==age)>0
        av_a_S[age,i]   = mean( a_sample[h_sample.==age] )                  ;
        pct_S[:,age,i]  = percentile( a_sample[h_sample.==age] , pct_list ) ;
        else 
        println(" Error in sample at i=$i, sample=$(sample_vec[i]) and age=$age: $(sum(h_sample.==age)) Observations")
        av_a_S[age,i]   = NaN ;
        pct_S[:,age,i] .= NaN ;
        end
        end 
    end 

    # Figure Average
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter(  20:(19+M_Aiyagari.p.Max_Age) , age_profile_a , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , av_a_S[:,1]   , marker=(:circle ,2.5 ), label=nothing)
    # scatter!( 20:(19+M_Aiyagari.p.Max_Age) , av_a_S[:,2]   , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , av_a_S[:,3]   , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , av_a_S[:,4]   , marker=(:circle ,2.5 ), label=nothing)
    ylims!(0,ceil(maximum(age_profile_a/10))*10)
    title!("Age Profile: Average Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(20:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_a_Simul.pdf")

    # Figure Percentiles
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    # 10th percentile 
    scatter(  20:(19+M_Aiyagari.p.Max_Age) , pct_C_age_1[1,:] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[1,:,1]     , marker=(:circle ,2.5 ), label=nothing)
    # scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[1,:,2]     , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[1,:,3]     , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[1,:,4]     , marker=(:circle ,2.5 ), label=nothing)
    # 90th percentile 
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_C_age_1[2,:] , marker=(:diamond ,3 ,:cornflowerblue),label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[3,:,1]     , marker=(:diamond ,2.5 ), label=nothing)
    # scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[3,:,2]     , marker=(:diamond ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[3,:,3]     , marker=(:diamond ,2.5 ), label=nothing)
    scatter!( 20:(19+M_Aiyagari.p.Max_Age) , pct_S[3,:,4]     , marker=(:diamond ,2.5 ), label=nothing)
    ylims!(0,ceil(maximum(pct_C_age_1[2,:]/10))*10)
    title!("Age Profile: 10th and 90th Percentiles of Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(20:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_a_pct_Simul.pdf")

## 45-year olds with more than median income
    age_0      = 26 ; 
    sample_vec = [10000, 50000, 100000, 500000]  ; 
    pct_list   = [10;50;90;99]                   ;
    pct_2_S    = zeros(4,M_Aiyagari.p.Max_Age-(age_0-1),4) ; 
    for i = 1:4
        a_sample    = M_C_45.a_mat[1:sample_vec[i],:]    ;
        h_sample    = M_C_45.h_mat[1:sample_vec[i],:]    ;
        ϵ_sample    = M_C_45.ϵ_mat[1:sample_vec[i],:]    ;
        ind_45      = ϵ_sample[:,1].>=med_ϵ              ;
        a_sample    = a_sample[ind_45,:]                 ;
        h_sample    = h_sample[ind_45,:]                 ;
        for age = 1:M_Aiyagari.p.Max_Age-(age_0-1)
        ind_alive   = h_sample[:,age].==1                ; 
        if sum(ind_alive)>0 
        pct_2_S[:,age,i]  = percentile( a_sample[ind_alive,age] , pct_list ) ;
        else 
        pct_2_S[:,age,i] .= NaN ; 
        end 
        end 
    end 

    # Figure Percentiles for 45 year old cohort
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    # 10th percentile 
    scatter(  45:(19+M_Aiyagari.p.Max_Age) , pct_C_age_2[1,:] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[1,:,1]   , marker=(:circle ,2.5 ), label=nothing)
    # scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[1,:,2]   , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[1,:,3]   , marker=(:circle ,2.5 ), label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[1,:,4]   , marker=(:circle ,2.5 ), label=nothing)
    # 90th percentile 
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_C_age_2[2,:] , marker=(:diamond ,3 ,:cornflowerblue),label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[3,:,1]   , marker=(:diamond ,2.5 ), label=nothing)
    # scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[3,:,2]   , marker=(:diamond ,2.5 ), label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[3,:,3]   , marker=(:diamond ,2.5 ), label=nothing)
    scatter!( 45:(19+M_Aiyagari.p.Max_Age) , pct_2_S[3,:,4]   , marker=(:diamond ,2.5 ), label=nothing)
    ylims!(0,ceil(maximum(pct_C_age_2[2,:]/10))*10)
    title!("Age Profile: 10th and 90th Percentiles of Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_a_pct_45_Simul.pdf")

println("===============================================\n")


###################################################################
###################################################################
## Auto-correlation of wealth Ages 35 and 65 (conditional on survival)
println("\n===============================================")
println("Autocorr of Wealth: ages 35-65")
    age_0        = 16                ; # Start agents at age 35 
    n_H          = 65-(35-1)         ; # Simulate for n_H years 
    cor_a_3565_S = zeros(fig_N)      ; # Cor of Consumption in first quintile 
    for i=1:fig_N
        # Fix current sample 
        a_aux_0       = M_C_35.a_mat[1:fig_sample[i],1]   ;
        a_aux_T       = M_C_35.a_mat[1:fig_sample[i],n_H] ;
        # Compute moments 
        cor_a_3565_S[i]  = cor( [a_aux_0 a_aux_T]  ;dims=1)[2]              ;  
    end 


    # Figure 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot( fig_sample/1000 , cor_a_3565_S ,label=nothing)
    hline!( [cor_a_3565] ,c=:orange  , w=3 , label=nothing )
    ylims!(0,1)
    xlims!(1,M_P.N_Panel/1000)
    title!("Auto-Correlation of Wealth, ages 35-65",titlefont=14)
    xlabel!("Sample Size: Thousands",labelsize=18)
    savefig("./"*Fig_Folder*"/A_Corr_3565_Simul.pdf")


println("===============================================\n")
    