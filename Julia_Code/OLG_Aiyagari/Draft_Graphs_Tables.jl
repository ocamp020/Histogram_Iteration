###################################################################
###################################################################
## Graphs and Tables 
    
color_vec_H = [:black , :lightskyblue2 , :cornflowerblue , :royalblue4] ;# range(colorant"cornflowerblue", stop=colorant"orange", length=4) ;


###################################################################
## Load results from csv files 
    H_grid_size =  [250 250 500 1000] ; # 100:50:1000 ; 
    n_H = length(H_grid_size) ;
    pct_list = [90;95;99;99.9;99.99] ;
    S_sample = [50000 100000 250000 500000] ;  # 10000:1000:M_Panel.N_Panel ; 
    N_S      = length(S_sample)          ;
    age_0_Wealth_Profile =  26           ; 
    age_0_Wealth_Corr    =  16           ; 
    age_T_Wealth_Corr    =  36           ; 
    
    H_Γ_timed =          readdlm(Hist_Folder*"/H_G_timed.csv", ',', Float64) ;
    H_Γ_bytes =          readdlm(Hist_Folder*"/H_G_bytes.csv", ',', Float64) ;
    H_M_timed = reshape( readdlm(Hist_Folder*"/H_M_timed.csv", ',', Float64) , n_H , 3 ) ;
    H_M_bytes = reshape( readdlm(Hist_Folder*"/H_M_bytes.csv", ',', Float64) , n_H , 3 ) ;

    H_Wealth_Profile_NB  = reshape( readdlm(Hist_Folder*"/H_Wealth_Profile_NB.csv", ',', Float64) , p.Max_Age , 6 , N_S ) ;
    H_Wealth_Profile_45  = reshape( readdlm(Hist_Folder*"/H_Wealth_Profile_45.csv", ',', Float64) , p.Max_Age-age_0_Wealth_Profile+1 , 6 , N_S ) ;
    H_Wealth_Corr        = reshape( readdlm(Hist_Folder*"/H_Wealth_Corr.csv"      , ',', Float64) , n_H , 2 );    
    
    S_M_timed = reshape( readdlm(MC_Folder*"/S_M_timed.csv", ',', Float64) , N_S , 6 ) ;
    S_M_bytes = reshape( readdlm(MC_Folder*"/S_M_bytes.csv", ',', Float64) , N_S , 4 ) ;

    S_Wealth_Profile_NB  = reshape( readdlm(MC_Folder*"/S_Wealth_Profile_NB.csv", ',', Float64) , p.Max_Age , 6 , N_S ) ;
    S_Wealth_Profile_45  = reshape( readdlm(MC_Folder*"/S_Wealth_Profile_45.csv", ',', Float64) , p.Max_Age , 6 , N_S ) ;
    S_Wealth_Corr        = reshape( readdlm(MC_Folder*"/S_Wealth_Corr.csv"      , ',', Float64) , N_S , 2 ) ;    
 
    
###################################################################
## Wealth Profiles Figures 
    y_tick_pr  = [0,250,500,750,1000,1250,1500,1750] ;
    y_label_pr = ["\$0m","\$0.25m","\$0.50m","\$0.75m","\$1.00m","\$1.25m","\$1.50m","\$1.75m"]; # ["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"]
    S_label    = ["100k","250k","350k","500k"] ;
    for i=1:N_S 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing,legendfontsize=9)
    scatter( [0],[0],marker=(:dtriangle ,15,:white),label=nothing) # label="Ave.  Wealth"
    scatter!([0],[0],marker=(:circle    ,15,:white),label=nothing) # label="p99   Wealth"
    scatter!([0],[0],marker=(:diamond   ,15,:white),label=nothing) # label="p99.9 Wealth"
    scatter!([0],[0],marker=(:rect,15,:cornflowerblue),label=nothing) # label="Histogram, N=500"
    scatter!([0],[0],marker=(:rect,15,:orangered2    ),label=nothing) # label="Monte-Carlo, N="*S_label[i]
    # Average 
    scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,2] , marker=(:dtriangle , 4 , 0.75 ,:cornflowerblue), markerstrokewidth=0 , markerstrokecolor=:cornflowerblue ,label=nothing)
    scatter!(  45:(19+p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,i] , marker=(:dtriangle , 4 , 0.95 ,:orangered2    ), markerstrokewidth=0.10 , markerstrokecolor=:grey25,label=nothing)
    # 90th percentile 
    scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end, 3 ,2] , marker=(:circle , 4.25 , 0.75 ,:cornflowerblue), markerstrokewidth=0 , markerstrokecolor=:cornflowerblue,label=nothing)
    scatter!(  45:(19+p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end, 3 ,i] , marker=(:circle , 4.25 , 0.95 ,:orangered2    ), markerstrokewidth=0.10 , markerstrokecolor=:grey25,label=nothing)
    # 99th percentile 
    scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end, 4 ,2] , marker=(:diamond, 4.25 , 0.75 ,:cornflowerblue), markerstrokewidth=0 , markerstrokecolor=:cornflowerblue,label=nothing)
    scatter!(  45:(19+p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end, 4 ,i] , marker=(:diamond, 4.25 , 0.95 ,:orangered2    ), markerstrokewidth=0.10 , markerstrokecolor=:grey25,label=nothing)
    # Formatting 
    ylims!(0,ceil(maximum(S_Wealth_Profile_NB[:,4,1]/250))*250)
    xlabel!("Age",labelsize=18); 
    yticks!(y_tick_pr,y_label_pr)# ylabel!("Thousands of Dollars",labelsize=18)
    xlims!(45,100); xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_$(Int(S_sample[i]/1000))k.pdf")
    end 

    # Comparing Histograms 
    alpha_vec = range(0.85,0.95,length=4) ;
    ms_vec    = range(4.50,3.25,length=4) ;
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing,legendfontsize=9)
    scatter( [0],[0],marker=(:dtriangle ,15,:white),label=nothing) # label="Ave.  Wealth"
    scatter!([0],[0],marker=(:circle    ,15,:white),label=nothing) # label="p99   Wealth"
    scatter!([0],[0],marker=(:diamond   ,15,:white),label=nothing) # label="p99.9 Wealth"
    for i=[2 3 4] 
        scatter!([0],[0],marker=(:rect,15,color_vec_H[i]),label=nothing) # label="Histogram, N=$(H_grid_size[i])"
        # Average 
        scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,i] , marker=(:dtriangle , ms_vec[i] , alpha_vec[i] ,color_vec_H[i]), markerstrokewidth=0 , markerstrokecolor=color_vec_H[i] ,label=nothing)
        # 90th percentile 
        scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end, 3 ,i] , marker=(:circle    , ms_vec[i] , alpha_vec[i] ,color_vec_H[i]), markerstrokewidth=0 , markerstrokecolor=color_vec_H[i] ,label=nothing)
        # 99th percentile 
        scatter!(  45:(19+p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end, 4 ,2] , marker=(:diamond   , ms_vec[i] , alpha_vec[i] ,color_vec_H[i]), markerstrokewidth=0 , markerstrokecolor=color_vec_H[i] ,label=nothing)
        
    end 
    # Formatting 
    ylims!(0,ceil(maximum(S_Wealth_Profile_NB[:,4,1]/250))*250)
    yticks!(y_tick_pr,y_label_pr); # ylabel!("Thousands of Dollars",labelsize=18)
    xlabel!("Age",labelsize=18); 
    xlims!(45,100); xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_Hist.pdf")

    # Figure for legend 
    step_1 = 2.6 ; step_2 = 3.5 ;
    x_vec_m = [collect(0.50:step_1:(0.5+step_1*3))' collect((0.5+step_1*4):step_2:(0.5+step_1*4)+step_2*2)'] ;  x_vec_a = x_vec_m.+0.30 ;
    gr(axis=([], false)) #axis=([], false) grid = false
    scatter( [x_vec_m[1]],[1.0],marker=(:dtriangle , 4.5 ,:white        ),label=nothing) ; annotate!( [x_vec_a[1]] , [1.0] , text("Ave. Wealth" ,:left,5)) ;
    scatter!([x_vec_m[2]],[1.0],marker=(:circle    , 4.5 ,:white        ),label=nothing) ; annotate!( [x_vec_a[2]] , [1.0] , text("p99 Wealth"  ,:left,5)) ;
    scatter!([x_vec_m[3]],[1.0],marker=(:diamond   , 4.5 ,:white        ),label=nothing) ; annotate!( [x_vec_a[3]] , [1.0] , text("p99.9 Wealth",:left,5)) ;
    scatter!([x_vec_m[4]],[1.0],marker=(:rect      , 4.5 ,:orangered2   ), markerstrokewidth=0.10 , markerstrokecolor=:grey25       ,label=nothing) ; annotate!( [x_vec_a[4]] , [1.0] , text("Monte-Carlo" ,:left,5)) ;
    scatter!([x_vec_m[5]],[1.0],marker=(:rect      , 4.5 ,color_vec_H[2]), markerstrokewidth=0    , markerstrokecolor=color_vec_H[2],label=nothing) ; annotate!( [x_vec_a[5]] , [1.0] , text("Histogram, N=$(H_grid_size[1])",:left,5)) ;
    scatter!([x_vec_m[6]],[1.0],marker=(:rect      , 4.5 ,color_vec_H[3]), markerstrokewidth=0    , markerstrokecolor=color_vec_H[3],label=nothing) ; annotate!( [x_vec_a[6]] , [1.0] , text("Histogram, N=$(H_grid_size[2])",:left,5)) ;
    scatter!([x_vec_m[7]],[1.0],marker=(:rect      , 4.5 ,color_vec_H[4]), markerstrokewidth=0    , markerstrokecolor=color_vec_H[4],label=nothing) ; annotate!( [x_vec_a[7]] , [1.0] , text("Histogram, N=$(H_grid_size[4])",:left,5)) ;
    plot!([0.0, 0.0, x_vec_m[end]+step_2+0.10 , x_vec_m[end]+step_2+0.10 , 0.0] , [0.97, 1.03, 1.03, 0.97, 0.97],color=:gray75,w=0.5,label=nothing)
    xlims!(-0.5,x_vec_m[end]+step_2+1.0); ylims!(0.5,1.5);
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_Legend.pdf")


###################################################################
## Autocorrelation 
    Mat_Corr = [ "Wealth Autocorrelation 35-65" "" "" "" "";
                "Histogram"      "" "" "" "";
                "Grid Size"         H_grid_size               ;
                "Ave. a(35)"        H_Wealth_Profile_NB[age_0_Wealth_Corr,end,:]' ;
                "p99  a(35)"        H_Wealth_Profile_NB[age_0_Wealth_Corr, 3 ,:]' ;
                "Ave. a(65)"        H_Wealth_Profile_NB[age_T_Wealth_Corr,end,:]' ;
                "p99  a(65)"        H_Wealth_Profile_NB[age_T_Wealth_Corr, 3 ,:]' ;
                "Wealth Corr 35-40" 100*H_Wealth_Corr[:,1]'        ; 
                "Wealth Corr 35-50" 100*H_Wealth_Corr[:,2]'        ; 
                "-" "-" "-" "-" "-";
                "Model Time"        H_Γ_timed'                ;
                "Profile Time"      H_M_timed[:,1]'           ;
                "Corr Time 35-45"   H_M_timed[:,2]'           ;
                "Corr Time 35-45"   H_M_timed[:,3]'           ;
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "Simulation"     "" "" "" "";
                "Sample Size (k)"   Int.(S_sample/1000)       ;
                "Ave. a(35)"        S_Wealth_Profile_NB[age_0_Wealth_Corr,end,:]' ;
                "p99  a(35)"        S_Wealth_Profile_NB[age_0_Wealth_Corr, 3 ,:]' ;
                "Ave. a(65)"        S_Wealth_Profile_NB[age_T_Wealth_Corr,end,:]' ;
                "p99  a(65)"        S_Wealth_Profile_NB[age_T_Wealth_Corr, 3 ,:]' ;
                "Wealth Corr 35-40" 100*S_Wealth_Corr[:,1]'        ; 
                "Wealth Corr 35-50" 100*S_Wealth_Corr[:,2]'        ; 
                "-" "-" "-" "-" "-";
                "Simulation Time"   S_M_timed[:,1]'           ;
                "Cohort Simul Time" S_M_timed[:,6]'           ;
                "Profile Time"      S_M_timed[:,2]'           ;
                "Corr Time 35-45"   S_M_timed[:,3]'           ;
                "Corr Time 35-55"   S_M_timed[:,4]'           ;
                "-" "-" "-" "-" "-";];

        open("./"*Fig_Folder*"/Table_Auto_Corr.csv", "w") do io
        writedlm(io, Mat_Corr, ',')
        end;

#=

    # N_Age
    N_Age_Short = M_Aiyagari.p.Max_Age-(age_0_Wealth_Profile-1) ;
    # 100k Observations 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
    scatter(0,0,marker=(:circle ,0.1 ,:white))
    # # Average 
    # scatter(   45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,2] , marker=(:circle ,2 ,:gray70),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,end,1] , marker=(:diamond,2 ,:gray70),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[           :            ,end,2] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[       1:N_Age_Short    ,end,1] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 90th percentile 
    scatter!(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,3,2] , marker=(:circle ,3 ,:gray70),label=nothing)
    scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,3,4] , marker=(:diamond,3 ,:gray70),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[           :            ,3,2] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[       1:N_Age_Short    ,3,1] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # 99th percentile 
    scatter!(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_NB[age_0_Wealth_Profile:end,4,2] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_NB[age_0_Wealth_Profile:end,4,4] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , H_Wealth_Profile_45[           :            ,4,2] , marker=(:circle ,3 ,:cornflowerblue),label=nothing)
    # scatter!(  45:(19+M_Aiyagari.p.Max_Age) , S_Wealth_Profile_45[       1:N_Age_Short    ,4,1] , marker=(:diamond,3 ,:cornflowerblue),label=nothing)
    # title!("Age Profile: Average, 90th, 99th Percentiles of Assets",titlefont=14)
    ylims!(0,ceil(maximum(S_Wealth_Profile_45[:,4,1]/10))*10)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(45:10:100)
    savefig("./"*Fig_Folder*"/Draft_Wealth_Profile_45_10k.pdf")



    # 100k Observations 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    # Average 
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


 =#