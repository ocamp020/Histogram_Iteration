###################################################################
###################################################################

## Define marginal distributions 
    Γ_a   = dropdims( sum( M_Aiyagari.Γ , dims=(3,2) ) , dims=(3,2) ) ; # Assets 
    Γ_ϵ   = dropdims( sum( M_Aiyagari.Γ , dims=(1,3) ) , dims=(1,3) ) ; # Labor Efficiency 
    Γ_age = dropdims( sum( M_Aiyagari.Γ , dims=(1,2) ) , dims=(1,2) ) ; # Returns 

## Define index for median shocks
    med_ϵ   = convert(Int64,round(M.n_ϵ/2)) ;
    ref_age = 31                            ;

###################################################################
###################################################################



###################################################################
###################################################################
## Labor income stats
println("===============================================")
    println("\n Income Grid and Probability")
    println("   Node - Value - PDF - Γ")
    aux = [M_Aiyagari.ϵ_grid*p.w M_Aiyagari.MP_ϵ.PDF Γ_ϵ]
    for i_ϵ=1:M_Aiyagari.n_ϵ
        println("   $i_ϵ:  $(round.(aux[i_ϵ,:],digits=4))")
    end 
    av_y = sum(M_Aiyagari.y_mat_fine.*M_Aiyagari.Γ) 
    sd_y = sqrt( sum( ((M_Aiyagari.y_mat_fine .- av_y).^2).*M_Aiyagari.Γ ) )
    age_profile_y = [ sum(M_Aiyagari.y_mat_fine[:,:,i].*M_Aiyagari.Γ[:,:,i])/Γ_age[i] for i in 1:M_Aiyagari.p.Max_Age ]         
    println("\n    Expected value: $(round( av_y ,digits=3)) Thousand Dollars")
    println("    Standard Deviation: $(round( sd_y , digits=3)) Thousand Dollars \n")
println("===============================================\n")

## Asset Distribution Stats 
    println("===============================================")
    println("\n Asset Distribution Stats")
    av_a     = sum( M_Aiyagari.a_grid_fine.*Γ_a ) ;
    sd_a     = sqrt( sum( ((M_Aiyagari.a_grid_fine[1:end] .- av_a).^2).*Γ_a )  ) ;
    CDF_a    = cumsum(Γ_a) ;
    Lorenz_a = cumsum(M_Aiyagari.a_grid_fine.*Γ_a)/av_a  
    Top_shares = zeros(5,3) ; 
        Top_shares[:,1] = [90;95;99;99.9;99.99];
        for i=1:5 
            ind_top = CDF_a.>=Top_shares[i,1]/100
            Top_shares[i,2] = M_Aiyagari.a_grid_fine[ind_top][1] ;
            Top_shares[i,3] = 100*sum( M_Aiyagari.a_grid_fine[ind_top].*Γ_a[ind_top] )/av_a ; 
        end
    age_profile_a  = [ sum(M_Aiyagari.a_mat_fine[:,:,i].*M_Aiyagari.Γ[:,:,i])/Γ_age[i] for i in 1:M_Aiyagari.p.Max_Age ]         
    age_profile_ap = [ sum(M_Aiyagari.G_ap_fine[:,:,i].*M_Aiyagari.Γ[:,:,i])/Γ_age[i] for i in 1:M_Aiyagari.p.Max_Age ]         
    B2A_ratio      = sum(age_profile_ap.*Γ_age.*(1 .- M_Aiyagari.p.Surv_Pr))/av_a
    println("    Expected Value: \$$(round.( av_a ,digits=2))k ")
    println("    Standard Deviation: \$$(round.( sd_a ,digits=2))k ")
    println("    Top X%  Level   Share ")
    for i=1:5
    println("    $(round(100-Top_shares[i,1],digits=2))%  \$$(round(Top_shares[i,2],digits=2))k   $(round(Top_shares[i,3],digits=2))% ")
    end 
    println("    Bequest to wealth ratio = $(round.( 100*B2A_ratio ,digits=2))% ")
println("===============================================\n")
###################################################################
###################################################################



###################################################################
###################################################################
## Plot Grid and Fine Grid (Zoom at the bottom and the top) 
l = @layout [a  ; b  c]
    # All the grid
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,ylims=(-0.25,1.25))
        p1 =
        scatter( M_Aiyagari.a_grid     ,zeros(M_Aiyagari.n_a)     ,marker=(:circle ,3,:cornflowerblue  ),label=nothing)
        scatter!(M_Aiyagari.a_grid_fine,ones( M_Aiyagari.n_a_fine),marker=(:diamond,3,:orange          ),label=nothing)
        #title!("Asset Grids",titlefont=14)
        xlabel!("Assets (thousands of dollars)",labelsize=18)
        yticks!([0,1],["Coarse","Fine"],tickfontsize=12)

    # Zoom in to lower end of the grid  
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,ylims=(-0.25,1.25))    
        p2=
        scatter( M_Aiyagari.a_grid     ,zeros(M_Aiyagari.n_a)     ,marker=(:circle ,3,:cornflowerblue  ),label=nothing)
        scatter!(M_Aiyagari.a_grid_fine,ones( M_Aiyagari.n_a_fine),marker=(:diamond,3,:orange          ),label=nothing)
        #title!("Asset Grids",titlefont=14)
        xlabel!("Assets (thousands of dollars)",labelsize=18)
        yticks!([0,1],["Coarse","Fine"],tickfontsize=12)
        xlims!(0,10)

    # Zoom in to higher end of the grid  
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,ylims=(-0.25,1.25))
        p3 = 
        scatter( M_Aiyagari.a_grid[1:end]./1000     ,zeros(M_Aiyagari.n_a)     ,marker=(:circle ,3,:cornflowerblue  ),label=nothing)
        scatter!(M_Aiyagari.a_grid_fine[1:end]./1000,ones( M_Aiyagari.n_a_fine),marker=(:diamond,3,:orange          ),label=nothing)
        #title!("Asset Grids",titlefont=14)
        xlabel!("Assets (millions of dollars)",labelsize=18)
        yticks!([0,1],["Coarse","Fine"],tickfontsize=12)
        xlims!(1,M_Aiyagari.a_max/1000)

    #Combine plots
        pjoint = plot(p1, p2, p3, layout = l)
        savefig(pjoint,"./"*Fig_Folder*"/Asset_Grids_Layout.pdf")
        # savefig(p1,"./"*Fig_Folder*"/Asset_Grids.pdf")
        # savefig(p2,"./"*Fig_Folder*"/Asset_Grids_Low.pdf")
        # savefig(p3,"./"*Fig_Folder*"/Asset_Grids_High.pdf")
###################################################################
###################################################################


###################################################################
###################################################################
# Plot Income Distribution
    med_ϵ = convert(Int64,round(M.n_ϵ/2));
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(M_Aiyagari.ϵ_grid*p.w) , 100*Γ_ϵ , marker=(:circle ,7,:cornflowerblue),label=nothing)
    vline!( [log.(M_Aiyagari.ϵ_grid[med_ϵ]*p.w)] ,c=:gray70  ,w=1,label=nothing)
    annotate!(log.(M_Aiyagari.ϵ_grid[med_ϵ]*p.w)+0.3,2,"\$$(round(M_Aiyagari.ϵ_grid[med_ϵ]*p.w,digits=1))k",10)
    ylims!(0,ceil(maximum(100*Γ_ϵ/10))*10)
    xlims!(log(1),log(1000)); xticks!(log.([1,5,10,50,100,500,1000]),["\$1k","\$5k","\$10k","\$50k","\$100k","\$500k","\$1m"])
    # xlims!(0,ceil(maximum(M_Aiyagari.ϵ_grid*p.w/500))*500)
    title!("Labor Income Distribution",titlefont=14)
    xlabel!("(log) Labor Income",labelsize=18)
    savefig("./"*Fig_Folder*"/Distribution_Income_Shocks.pdf")
###################################################################
###################################################################



###################################################################
###################################################################
## Plot asset distribution 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(M_Aiyagari.a_grid_fine) , 100*Γ_a , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    xlabel!("Log Assets",labelsize=18)
    title!("Asset Distribution",titlefont=14)
    ylims!(0,ceil(maximum(100*Γ_a/1))*1)
    xlims!(log(0.1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1000)); 
    xticks!(log.([1,10,100,1000,10000,100000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$100m"])
    savefig("./"*Fig_Folder*"/Distribution_Wealth.pdf")

## Plot CDF
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(M_Aiyagari.a_grid_fine) , 100*CDF_a , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    xlabel!("Log Assets",labelsize=18)
    title!("Cumulative Asset Distribution",titlefont=14)
    ylims!(0,100)
    xlims!(log(0.1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1000)); 
    xticks!(log.([1,10,100,1000,10000,100000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$100m"])
    savefig("./"*Fig_Folder*"/Distribution_Wealth_CDF.pdf")

## Plot Pareto Tail (Above $1 Million)
    ind     = M_Aiyagari.a_grid_fine.>=1000 ;
    grid_1M = M_Aiyagari.a_grid_fine[ind]   ;
    Γ_a_1M  = Γ_a[ind]/sum(Γ_a[ind])        ; Γ_a_1M = Γ_a_1M/sum(Γ_a_1M) ; 
    CCDF_1M = 1 .- cumsum(Γ_a_1M)           ; 
    ind     = CCDF_1M.>=1e-12               ;
    P_coeff = (log.(grid_1M[ind]./1000)'*log.(grid_1M[ind]./1000))\log.(grid_1M[ind]./1000)'*log.(CCDF_1M[ind])
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( log.(grid_1M[ind]./1000) , log.(CCDF_1M[ind]) , marker=(:circle ,3,:cornflowerblue) , markerstrokewidth=0 , label=nothing )   
    plot!( log.(grid_1M[1:end-1]./1000) , P_coeff.*log.(grid_1M[1:end-1]./1000) , w=2, c=:orange , label=nothing )
    annotate!(-log(1.3)+mean(log.(grid_1M[ind]./1000)),mean(log.(CCDF_1M[ind])),"α=$(round(P_coeff,digits=2))",12)
    xlabel!("Log Assets",labelsize=18)
    title!("Distribution Tail",titlefont=14)
    ylims!( floor(minimum(log.(CCDF_1M[ind]))/4)*4 , 0 )
    xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
    xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
    savefig("./"*Fig_Folder*"/Distribution_Wealth_Pareto.pdf")

## Plot Lorenz Curve
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(1:100,1:100,w=1,linecolor=:gray70,label=nothing,aspect_ratio=1)
    plot!( 100*CDF_a  , 100*Lorenz_a , w=3, c=:cornflowerblue , label=nothing ,aspect_ratio=1)   
    title!("Lorenz curve",titlefont=14)
    ylims!(0,100); xlims!(0,100)
    savefig("./"*Fig_Folder*"/Distribution_Wealth_Lorenz.pdf")

###################################################################
###################################################################



###################################################################
###################################################################
## Age Profiles 
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( 20:(19+M_Aiyagari.p.Max_Age) , age_profile_y , marker=(:circle ,5,:cornflowerblue),label=nothing)
    ylims!(0,ceil(maximum(age_profile_y/10))*10)
    title!("Age Profile: Labor Income",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(20:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_Y.pdf")

    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( 20:(19+M_Aiyagari.p.Max_Age) , age_profile_a , marker=(:circle ,5,:cornflowerblue),label=nothing)
    ylims!(0,ceil(maximum(age_profile_a/10))*10)
    title!("Age Profile: Assets",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Thousands of Dollars",labelsize=18)
    xticks!(20:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_a.pdf")

    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    scatter( 20:(19+M_Aiyagari.p.Max_Age) , 100*(age_profile_ap./age_profile_a.-1) , marker=(:circle ,5,:cornflowerblue),label=nothing)
    ylims!(floor(minimum(100*(age_profile_ap./age_profile_a.-1)/5))*5,ceil(maximum(100*(age_profile_ap./age_profile_a.-1)/5))*5)
    title!("Age Profile: (agg) Savings Rate",titlefont=14)
    xlabel!("Age",labelsize=18); ylabel!("Percentage Points",labelsize=18)
    xticks!(20:10:100)
    savefig("./"*Fig_Folder*"/Age_Profile_ap.pdf")
###################################################################
###################################################################



###################################################################
###################################################################
## Plot Saving Functions (median labor efficiency and interest rate)
    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    plot(M_Aiyagari.a_grid_fine,M_Aiyagari.a_grid_fine,w=1,linecolor=:gray70,label=nothing,aspect_ratio=1,xlims=(M_Aiyagari.a_grid[1],M_Aiyagari.a_grid[end]))
    plot!(M_Aiyagari.a_grid_fine,M_Aiyagari.G_ap_fine[:,med_ϵ,ref_age],w=2,linecolor=:cornflowerblue,label=nothing,aspect_ratio=1)
    title!("Savings",titlefont=14)
    xlabel!("Assets (thousands of dollars)",labelsize=18)
    ylabel!("Assets (thousands of dollars)",labelsize=18)
    xlims!(0,M_Aiyagari.a_max)
    ylims!(0,M_Aiyagari.a_max)
    savefig("./"*Fig_Folder*"/Policy_Function_Savings_Level.pdf")

    gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
    hline( [0] ,c=:gray70  ,w=1,label=nothing) 
    plot!(log.(M_Aiyagari.a_grid_fine) , 100*(M_Aiyagari.G_ap_fine[:,med_ϵ,ref_age]./M_Aiyagari.a_grid_fine.-1),w=2,linecolor=:cornflowerblue,label=nothing,xlims=(M_Aiyagari.a_grid[1],M_Aiyagari.a_grid[end]))
    title!("Savings Rate",titlefont=14)
    xlabel!("(log) Assets",labelsize=18)
    ylabel!("Saving Rate (%)",labelsize=18)
    xlims!(log(0.8),log(M_Aiyagari.a_max)); xticks!(log.([1,10,100,1000,10000,50000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$50m"])
    ylims!(-15,50)
    savefig("./"*Fig_Folder*"/Policy_Function_Savings_Rate.pdf")
###################################################################
###################################################################


###################################################################
###################################################################
## Euler Errors plots
    Euler_perc_error = zeros(M_Aiyagari.n_a_fine,3)
    for i_a=1:M_Aiyagari.n_a_fine
        # Euler percentage error (lowest points in ϵ)
        Euler_perc_error[i_a,1] = Euler_Error(1             ,ref_age,M_Aiyagari.a_grid_fine[i_a],M_Aiyagari)
        # Euler percentage error (median points in ϵ) 
        Euler_perc_error[i_a,2] = Euler_Error(med_ϵ         ,ref_age,M_Aiyagari.a_grid_fine[i_a],M_Aiyagari)
        # Euler percentage error (higest points in ϵ)
        Euler_perc_error[i_a,3] = Euler_Error(M_Aiyagari.n_ϵ,ref_age,M_Aiyagari.a_grid_fine[i_a],M_Aiyagari)
    end
    
    # Plot: lowest ϵ
    ll = @layout [a ; b ; c]
    gr(ytickfontsize=10,xtickfontsize=10,xtick_direction=:out)
    p1_euler=
    scatter(log.(M_Aiyagari.a_grid_fine), Euler_perc_error[:,1],marker=(:circle ,3,:cornflowerblue),label=nothing)
    ylabel!("Percentage")
    xlims!(log(0.8),log(M_Aiyagari.a_max)+0.2); xticks!(log.([1,10,100,1000,10000,50000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$50m"])
    ylims!(findmin(Euler_perc_error[:,1])[1],findmax(Euler_perc_error[:,1])[1])
    
    # Plot: median ϵ
    p2_euler=
    scatter(log.(M_Aiyagari.a_grid_fine), Euler_perc_error[:,2],marker=(:circle ,3,:cornflowerblue),label=nothing)
    ylabel!("Percentage")
    xlims!(log(0.8),log(M_Aiyagari.a_max)+0.2); xticks!(log.([1,10,100,1000,10000,50000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$50m"])
    ylims!(findmin(Euler_perc_error[:,2])[1],findmax(Euler_perc_error[:,2])[1])
    
    # Plot: highest ϵ
    p3_euler=
    scatter(log.(M_Aiyagari.a_grid_fine), Euler_perc_error[:,3],marker=(:circle ,3,:cornflowerblue),label=nothing)
    ylabel!("Percentage")
    xlims!(log(0.8),log(M_Aiyagari.a_max)+0.2); xticks!(log.([1,10,100,1000,10000,50000]),["\$1k","\$10k","\$100k","\$1m","\$10m","\$50m"])
    ylims!(findmin(Euler_perc_error[:,3])[1],findmax(Euler_perc_error[:,3])[1])
    
    # Combine plots
    pjoint_euler = plot(p1_euler, p2_euler, p3_euler, layout = ll)
    savefig(pjoint_euler,"./"*Fig_Folder*"/Euler_Error_Layout.pdf")
    # savefig(p1_euler,"./"*Fig_Folder*"/Euler_Error_low.pdf")
    # savefig(p2_euler,"./"*Fig_Folder*"/Euler_Error_median.pdf")
    # savefig(p3_euler,"./"*Fig_Folder*"/Euler_Error_high.pdf")
###################################################################
###################################################################