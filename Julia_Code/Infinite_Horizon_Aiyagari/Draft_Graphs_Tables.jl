###################################################################
###################################################################
## Graphs and Tables 

###################################################################
## Set up colors
color_vec_H = range(colorant"cornflowerblue", stop=colorant"orange", length=4) ;
color_vec_S = range(colorant"darkgoldenrod", stop=colorant"green", length=4) ;

###################################################################
## Set up model structures 
M_Simul = Model(method=1,read_flag=true) ;
M_Panel = Model_Panel(N_Panel=1000000)   ; 
    
###################################################################
## Load results from csv files 
H_grid_size = [250 500 750 1000] ;  
n_H = length(H_grid_size) ;
pct_list = [90;95;99;99.9;99.99] ;
S_sample = 250000:250000:1000000 ; 
N_S      = length(S_sample)          ;

H_Γ_timed =          readdlm(Hist_Folder*"/H_G_timed.csv", ',', Float64) ;
H_Γ_bytes =          readdlm(Hist_Folder*"/H_G_bytes.csv", ',', Float64) ;
H_M_timed = reshape( readdlm(Hist_Folder*"/H_M_timed.csv", ',', Float64) , n_H , 3 ) ;
H_M_bytes = reshape( readdlm(Hist_Folder*"/H_M_bytes.csv", ',', Float64) , n_H , 3 ) ;

H_a_grid  = reshape( readdlm(Hist_Folder*"/H_a_grid.csv", ',', Float64) , H_grid_size[end] , n_H ) ;
H_Γ_a     = reshape( readdlm(Hist_Folder*"/H_G_a.csv"   , ',', Float64) , H_grid_size[end] , n_H ) ;

H_Wealth_Stats  = reshape( readdlm(Hist_Folder*"/H_Wealth_Stats.csv", ',', Float64) , n_H , 6 ) ;
H_Wealth_Share  = reshape( readdlm(Hist_Folder*"/H_Wealth_Share.csv", ',', Float64) , n_H , 5 ) ;
H_Pareto_Coeff  =          readdlm(Hist_Folder*"/H_Pareto_Coeff.csv", ',', Float64) ;    
H_Decile        = reshape( readdlm(Hist_Folder*"/H_Decile.csv"   , ',', Float64) , 11 , 2  , n_H ) ;
H_Decile_Tr     = reshape( readdlm(Hist_Folder*"/H_Decile_Tr.csv", ',', Float64) , 10 , 10 , n_H ) ;
H_Cons_Corr     =          readdlm(Hist_Folder*"/H_Cons_Corr.csv", ',', Float64) ;
H_A_Corr        =          readdlm(Hist_Folder*"/H_A_Corr.csv"   , ',', Float64) ;
H_ϵ_Corr        =          readdlm(Hist_Folder*"/H_eps_Corr.csv" , ',', Float64) ;
H_ζ_Corr        =          readdlm(Hist_Folder*"/H_z_Corr.csv"   , ',', Float64) ;

S_M_timed = reshape( readdlm(MC_Folder*"/S_M_timed.csv", ',', Float64) , N_S , 4 ) ;
S_M_bytes = reshape( readdlm(MC_Folder*"/S_M_bytes.csv", ',', Float64) , N_S , 4 ) ;
S_Γ_timed =          readdlm(MC_Folder*"/S_G_timed.csv", ',', Float64)

S_Wealth_Sample = reshape( readdlm(MC_Folder*"/S_Wealth_Sample.csv", ',', Float64) , N_S , 1000000 ) ;
S_Wealth_Stats  = reshape( readdlm(MC_Folder*"/S_Wealth_Stats.csv" , ',', Float64) , N_S , 6 ) ;
S_Wealth_Share  = reshape( readdlm(MC_Folder*"/S_Wealth_Share.csv" , ',', Float64) , N_S , 5 ) ;
S_Pareto_Coeff  =          readdlm(MC_Folder*"/S_Pareto_Coeff.csv" , ',', Float64) ;    
S_Decile        = reshape( readdlm(MC_Folder*"/S_Decile.csv"   , ',', Float64) , 11 , N_S ) ;
S_Decile_Tr     = reshape( readdlm(MC_Folder*"/S_Decile_Tr.csv", ',', Float64) , 10 , 10 , N_S ) ;
S_Cons_Corr     =          readdlm(MC_Folder*"/S_Cons_Corr.csv", ',', Float64) ;
S_A_Corr        =          readdlm(MC_Folder*"/S_A_Corr.csv"   , ',', Float64) ;
S_ϵ_Corr        =          readdlm(MC_Folder*"/S_eps_Corr.csv" , ',', Float64) ;
S_ζ_Corr        =          readdlm(MC_Folder*"/S_z_Corr.csv"   , ',', Float64) ;



###################################################################
# Top Shares + Pareto Coefficient Table      
Mat_Top_Stats = [ "Top Wealth Shares + Pareto Tail" "" "" "" "";
                "Histogram"     "" "" "" "";
                "Grid Size"     H_grid_size                 ;
                "Top 0.1%"      H_Wealth_Share[:,4]'        ;
                "Top 1%  "      H_Wealth_Share[:,3]'        ;
                "Top 10% "      H_Wealth_Share[:,1]'        ;
                "Pareto  "      H_Pareto_Coeff'             ;
                "Av. Assets"    H_Wealth_Stats[:,end]'      ;
                "-" "-" "-" "-" "-";
                "Total Time"    H_Γ_timed'.+H_M_timed[:,1]' ;
                "Model Time"    H_Γ_timed'                  ;
                "Moment Time"   H_M_timed[:,1]'             ;
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "Simulation"     "" "" "" "";
                "Grid Size"     "250k" "500k" "750k" "1M"   ;
                "Top 0.1%"      S_Wealth_Share[:,4]'        ;
                "Top 1%  "      S_Wealth_Share[:,3]'        ;
                "Top 10% "      S_Wealth_Share[:,1]'        ;
                "Pareto  "      S_Pareto_Coeff'             ;
                "Av. Assets"    S_Wealth_Stats[:,end]'      ;
                "-" "-" "-" "-" "-";
                "Total Time"    S_Γ_timed.+S_M_timed[:,1]'.+S_M_timed[:,2]'  ;
                "Model Time"    repeat([S_Γ_timed],1,N_S)   ;
                "Simul Time"    S_M_timed[:,1]'             ;
                "Moment Time"   S_M_timed[:,2]'             ;
                "-" "-" "-" "-" "-";];

        open("./"*Fig_Folder*"/Table_Top_Stats.csv", "w") do io
        writedlm(io, Mat_Top_Stats, ',')
        end;


###################################################################
## Top 10%, 1%, 0.1% Shares 
#=
# Top 1%
gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
plot(  H_M_timed[:,1] , H_Wealth_Share[:,3] ,label="Top 1% - Histogram")
plot!( S_M_timed[:,1] , S_Wealth_Share[:,3] ,label="Top 1% - Simulation")
ylims!(0,100)
# xlims!(1,M_P.N_Panel/1000)
title!("Top 1% Share of Wealth",titlefont=14)
xlabel!("Time: Seconds",labelsize=18)
savefig("./"*Fig_Folder*"/Draft_Top_1_Share.pdf")

# Top 0.1%
gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
plot(  H_M_timed[:,1] , H_Wealth_Share[:,4] ,label="Top 0.1% - Histogram")
plot!( S_M_timed[:,1] , S_Wealth_Share[:,4] ,label="Top 0.1% - Simulation")
ylims!(0,100)
# xlims!(1,M_P.N_Panel/1000)
title!("Top 0.1% Share of Wealth",titlefont=14)
xlabel!("Time: Seconds",labelsize=18)
savefig("./"*Fig_Folder*"/Draft_Top_01_Share.pdf")

# Top 0.01%
gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out)
plot(  H_M_timed[:,1] , H_Wealth_Share[:,5] ,label="Top 0.01% - Histogram")
plot!( S_M_timed[:,1] , S_Wealth_Share[:,5] ,label="Top 0.01% - Simulation")
ylims!(0,100)
# xlims!(1,M_P.N_Panel/1000)
title!("Top 0.01% Share of Wealth",titlefont=14)
xlabel!("Time: Seconds",labelsize=18)
savefig("./"*Fig_Folder*"/Draft_Top_001_Share.pdf")
=#

###################################################################
## Pareto Tail 
# Results from Histogram (comparind different grid sizes)
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
        plot(log.(M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000]),M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000] , label=nothing )
        for i=1:n_H 
                ind     = H_a_grid[1:H_grid_size[i],i].>=1000 ;
                grid_1M = H_a_grid[1:H_grid_size[i],i][ind]   ;
                Γ_a     = H_Γ_a[1:H_grid_size[i],i] ; # Assets 
                Γ_a_1M  = Γ_a[ind]/sum(Γ_a[ind])     ; Γ_a_1M = Γ_a_1M/sum(Γ_a_1M) ; 
                CCDF_1M = 1 .- cumsum(Γ_a_1M)        ;
                scatter!( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:diamond ,3,0.75,color_vec_H[i]) , markerstrokewidth=0 , label=L"\alpha_{H}=%$(round(H_Pareto_Coeff[i],digits=2)),\,N=%$(H_grid_size[i])" )   
                plot!(    log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[i].*log.(grid_1M[1:end-1]./1000) , w=1.1, c=color_vec_H[i] , label=nothing )
                # annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α_H=$(round(H_Pareto_Coeff[i],digits=2))",11)
                ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
        end 
        xlabel!("Log Assets",labelsize=18)
        ylabel!("Log Counter CDF",labelsize=18)
        title!("Distribution Tail - Histogram",titlefont=14)
        xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
        xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
        savefig("./"*Fig_Folder*"/Draft_Pareto_Hist.pdf")

# Results from Monte Carlo (comparind different simulation sizes)
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
        plot(log.(M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000]),M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000] , label=nothing )
        for i=1:N_S 
                ind     = S_Wealth_Sample[i,1:S_sample[i]].>=1000 ;
                grid_1M = sort( S_Wealth_Sample[i,1:S_sample[i]][ind] )  ;
                CCDF_1M = collect(length(grid_1M):-1:1)./length(grid_1M) ; # Counter CDF = 1- CDF
                scatter!( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,0.75,color_vec_S[i]) , markerstrokewidth=0 , label=L"\alpha_{S}=%$(round(S_Pareto_Coeff[i],digits=2)),\,N=%$(Int(S_sample[i]/1000))k" )   
                plot!(    log.(grid_1M[1:end-1]./1000) , S_Pareto_Coeff[i].*log.(grid_1M[1:end-1]./1000) , w=1.1, c=color_vec_S[i] , label=nothing )
                # annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α_H=$(round(S_Pareto_Coeff[i],digits=2))",11)
                ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
        end 
        xlabel!("Log Assets",labelsize=18)
        ylabel!("Log Counter CDF",labelsize=18)
        title!("Distribution Tail - Simulation",titlefont=14)
        xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
        xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
        savefig("./"*Fig_Folder*"/Draft_Pareto_Simul.pdf")
        # Add Histogram with 500 grid points 
        ind_H     = H_a_grid[1:H_grid_size[2],2].>=1000 ;
        grid_1M_H = H_a_grid[1:H_grid_size[2],2][ind_H]   ;
        Γ_a_H     = H_Γ_a[1:H_grid_size[2],2] ; # Assets 
        Γ_a_1M_H  = Γ_a_H[ind_H]/sum(Γ_a_H[ind_H])     ; Γ_a_1M_H = Γ_a_1M_H/sum(Γ_a_1M_H) ; 
        CCDF_1M_H = 1 .- cumsum(Γ_a_1M_H)        ;
        scatter!( log.(grid_1M_H[1:end-1]./1000) , log.(CCDF_1M_H[1:end-1]) , marker=(:diamond ,3,0.5,color_vec_H[1]) , markerstrokewidth=0 , label=L"\alpha_{H}=%$(round(H_Pareto_Coeff[2],digits=2)),\,N=%$(H_grid_size[2])" )   
        plot!(    log.(grid_1M_H[1:end-1]./1000) , H_Pareto_Coeff[2].*log.(grid_1M_H[1:end-1]./1000) , w=1.1, c=color_vec_H[1] , label=nothing )
        savefig("./"*Fig_Folder*"/Draft_Pareto_Simul_vs_Hist.pdf")

        # Separate plots 
        for i=1:N_S 
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
        plot(log.(M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000]),M_Aiyagari.a_grid[M_Aiyagari.a_grid.<1000] , label=nothing )
                ind     = S_Wealth_Sample[i,1:S_sample[i]].>=1000 ;
                grid_1M = sort( S_Wealth_Sample[i,1:S_sample[i]][ind] )  ;
                CCDF_1M = collect(length(grid_1M):-1:1)./length(grid_1M) ; # Counter CDF = 1- CDF
                scatter!( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:circle ,3,0.75,color_vec_S[i]) , markerstrokewidth=0 , label=L"\alpha_{S}=%$(round(S_Pareto_Coeff[i],digits=2)),\,N=%$(Int(S_sample[i]/1000))k" )   
                plot!(    log.(grid_1M[1:end-1]./1000) , S_Pareto_Coeff[i].*log.(grid_1M[1:end-1]./1000) , w=1.1, c=color_vec_S[i] , label=nothing )
                # annotate!(-log(1.3)+mean(log.(grid_1M[1:end-1]./1000)),mean(log.(CCDF_1M[1:end-1])),"α_H=$(round(S_Pareto_Coeff[i],digits=2))",11)
                ylims!( floor(log(CCDF_1M[end-1])/4)*4 , 0 )
        
                # Add Histogram with 500 grid points 
                ind     = H_a_grid[1:H_grid_size[2],2].>=1000 ;
                grid_1M = H_a_grid[1:H_grid_size[2],2][ind]   ;
                Γ_a     = H_Γ_a[1:H_grid_size[2],2] ; # Assets 
                Γ_a_1M  = Γ_a[ind]/sum(Γ_a[ind])     ; Γ_a_1M = Γ_a_1M/sum(Γ_a_1M) ; 
                CCDF_1M = 1 .- cumsum(Γ_a_1M)        ;
                scatter!( log.(grid_1M[1:end-1]./1000) , log.(CCDF_1M[1:end-1]) , marker=(:diamond ,3,0.5,color_vec_H[1]) , markerstrokewidth=0 , label=L"\alpha_{H}=%$(round(H_Pareto_Coeff[2],digits=2)),\,N=%$(H_grid_size[2])" )   
                plot!(    log.(grid_1M[1:end-1]./1000) , H_Pareto_Coeff[2].*log.(grid_1M[1:end-1]./1000) , w=1.1, c=color_vec_H[1] , label=nothing )
        xlabel!("Log Assets",labelsize=18)
        ylabel!("Log Counter CDF",labelsize=18)
        title!("Distribution Tail - Simulation",titlefont=14)
        xlims!(log(1),log(ceil(M_Aiyagari.a_grid[end]/1000)*1)); 
        xticks!(log.([1,2,4,8,20,40,80]),["\$1m","\$2m","\$4m","\$8m","\$20m","\$40m","\$80m"])
        savefig("./"*Fig_Folder*"/Draft_Pareto_Simul_$(Int(S_sample[i]/1000))k.pdf")
        end 


###################################################################
## Decile Transition

Mat_Decile_TR = [ "Decile Transition" "" "" "" "";
                "Histogram"     "" "" "" "";
                "Grid Size"     H_grid_size                 ;
                "D1-D1"         H_Decile_Tr[1,1,:]'         ;
                "D1-D2"         H_Decile_Tr[1,2,:]'         ;
                "D2-D1"         H_Decile_Tr[2,1,:]'         ;
                "D2-D2"         H_Decile_Tr[2,2,:]'         ;
                "D10-D10"       H_Decile_Tr[10,10,:]'       ;
                "-" "-" "-" "-" "-";
                "Total Time"    H_Γ_timed'.+H_M_timed[:,2]' ;
                "Model Time"    H_Γ_timed'                  ;
                "Moment Time"   H_M_timed[:,2]'             ;
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "Simulation"     "" "" "" "";
                "Grid Size"     "250k" "500k" "750k" "1M"   ;
                "D1-D1"         S_Decile_Tr[1,1,:]'         ;
                "D1-D2"         S_Decile_Tr[1,2,:]'         ;
                "D2-D1"         S_Decile_Tr[2,1,:]'         ;
                "D2-D2"         S_Decile_Tr[2,2,:]'         ;
                "D10-D10"       S_Decile_Tr[10,10,:]'       ;
                "-" "-" "-" "-" "-";
                "Total Time"    S_Γ_timed.+S_M_timed[:,1]'.+S_M_timed[:,3]'  ;
                "Model Time"    repeat([S_Γ_timed],1,N_S);
                "Simul Time"    S_M_timed[:,1]'             ;
                "Moment Time"   S_M_timed[:,3]'             ;
                "-" "-" "-" "-" "-";];

        open("./"*Fig_Folder*"/Table_Decile_Tr.csv", "w") do io
        writedlm(io, Mat_Decile_TR, ',')
        end;

###################################################################
## Decile Transition
Mat_Decile = [ "Deciles"        "" "" "" "";
                "Histogram"     "" "" "" "";
                "Grid Size"     H_grid_size                 ;
                "D0"            H_Decile[1,2,:]'            ;
                "D1"            H_Decile[2,2,:]'            ;
                "D2"            H_Decile[3,2,:]'            ;
                "D3"            H_Decile[4,2,:]'            ;
                "D4"            H_Decile[5,2,:]'            ;
                "D5"            H_Decile[6,2,:]'            ;
                "D6"            H_Decile[7,2,:]'            ;
                "D7"            H_Decile[8,2,:]'            ;
                "D8"            H_Decile[9,2,:]'            ;
                "D9"            H_Decile[10,2,:]'           ;
                "-" "-" "-" "-" "-";
                "Total Time"    H_Γ_timed'.+H_M_timed[:,2]' ;
                "Model Time"    H_Γ_timed'                  ;
                "Moment Time"   H_M_timed[:,2]'             ;
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "Simulation"     "" "" "" "";
                "Grid Size"     "250k" "500k" "750k" "1M"   ;
                "D0"            S_Decile[1,:]'              ;
                "D1"            S_Decile[2,:]'              ;
                "D2"            S_Decile[3,:]'              ;
                "D3"            S_Decile[4,:]'              ;
                "D4"            S_Decile[5,:]'              ;
                "D5"            S_Decile[6,:]'              ;
                "D6"            S_Decile[7,:]'              ;
                "D7"            S_Decile[8,:]'              ;
                "D8"            S_Decile[9,:]'              ;
                "D9"            S_Decile[10,:]'             ;
                "-" "-" "-" "-" "-";
                "Total Time"    S_Γ_timed.+S_M_timed[:,1]'.+S_M_timed[:,3]'  ;
                "Model Time"    repeat([S_Γ_timed],1,N_S);
                "Simul Time"    S_M_timed[:,1]'             ;
                "Moment Time"   S_M_timed[:,3]'             ;
                "-" "-" "-" "-" "-";];

        open("./"*Fig_Folder*"/Table_Decile.csv", "w") do io
        writedlm(io, Mat_Decile, ',')
        end;


#= 
        # Persistence first decile 
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
        plot(  H_M_timed[:,2] , 100*H_Decile_Tr[1,1,:] ,label="d1-d1 Transition - Histogram")
        plot!( S_M_timed[:,2] , 100*S_Decile_Tr[1,1,:] ,label="d1-d1 Transition - Simulation")
        ylims!(0,1)
        # xlims!(1,M_P.N_Panel/1000)
        # title!("Top 1% Share of Wealth",titlefont=14)
        xlabel!("Time: Seconds",labelsize=18)
        savefig("./"*Fig_Folder*"/Draft_Decile_Tr_1.pdf")

        # Persistence tenth decile 
        gr(ytickfontsize=12,xtickfontsize=12,xtick_direction=:out,foreground_color_legend = nothing,background_color_legend = nothing)
        plot(  H_M_timed[:,2] , 100*dropdims(H_Decile_Tr[10,10,:]) ,label="d10-d10 Transition - Histogram")
        plot!( S_M_timed[:,2] , 100*dropdims(S_Decile_Tr[10,10,:]) ,label="d10-d10 Transition - Simulation")
        ylims!(0,1)
        # xlims!(1,M_P.N_Panel/1000)
        # title!("Top 1% Share of Wealth",titlefont=14)
        xlabel!("Time: Seconds",labelsize=18)
        savefig("./"*Fig_Folder*"/Draft_Decile_Tr_10.pdf")
=#


###################################################################
## Consumption 3 Year Auto-Correlation 

Mat_Cons_Corr = [ "Autocorrelations" "" "" "" "";
                "Histogram"      "" "" "" "";
                "Grid Size"      H_grid_size                 ;
                "Cons    Auto-Corr" H_Cons_Corr'             ;
                "Assets  Auto-Corr" H_A_Corr'                ;
                "Epsilon Auto-Corr" H_ϵ_Corr'                ;
                "Z       Auto-Corr" H_ζ_Corr'                ;
                "-" "-" "-" "-" "-";
                "Total Time"     H_Γ_timed'.+H_M_timed[:,3]' ;
                "Model Time"     H_Γ_timed'                  ;
                "Moment Time"    H_M_timed[:,3]'             ;
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "-" "-" "-" "-" "-";
                "Simulation"     "" "" "" "";
                "Grid Size"      "250k" "500k" "750k" "1M"   ;
                "Cons    Auto-Corr" S_Cons_Corr'             ;
                "Assets  Auto-Corr" S_A_Corr'                ;
                "Epsilon Auto-Corr" S_ϵ_Corr'                ;
                "Z       Auto-Corr" S_ζ_Corr'                ;
                "-" "-" "-" "-" "-";
                "Total Time"     S_Γ_timed.+S_M_timed[:,1]'.+S_M_timed[:,4]'  ;
                "Model Time"     repeat([S_Γ_timed],1,N_S);
                "Simul Time"     S_M_timed[:,1]'             ;
                "Moment Time"    S_M_timed[:,4]'             ;
                "-" "-" "-" "-" "-";];

        open("./"*Fig_Folder*"/Table_Auto_Corr.csv", "w") do io
        writedlm(io, Mat_Cons_Corr, ',')
        end;

