###################################################################
###################################################################
###################################################################
## Functions for Monte-Carlo simulation 


###################################################################
###################################################################
## Generate structure for the panel using Parameters module
@with_kw struct Model_Panel
    # Size 
    N_Panel::Int64 = 500000 ; # Number of dynasties to be simulated 
    T_Panel::Int64 = 10      ; # Number of periods to be saved of the dynasties 
    T_Simul::Int64 = 1000    ; # Number of periods to be simulated 
    N_Min::Int64   = 1000    ; # Minimum numberof dynasties for moments
    # Panel output 
    a_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Assets
    c_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Consumption 
    ϵ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Labor efficiency index
    ζ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Return index
    # Time 
    t_vec = Array{Float64}(undef,N_Panel)         ; # Simulation Time 
    # Tolerance for cross-sectional moments 
    Simul_tol = 1E-6  ;
    # Random number seed 
    rng_seed::Int64 = 3489398  ;
end





###################################################################
###################################################################
## Simulate the panel 
function Simulate_Panel(M::Model,M_P::Model_Panel,Seed_Flag::bool)
    ## Unpack relevant objects from model 
    @unpack p, MP_ϵ, MP_ζ, a_grid, a_grid_fine, a_max, Γ, G_ap, G_ap_fine, G_c, n_a_fine, θ_a_f = M
    @unpack a_min = p
    @unpack N_Panel, T_Panel, T_Simul, Simul_tol = M_P

    ## Initialize Seed 
    Random.seed!(M_P.rng_seed) ;

    ## Initialize cross-section vectors
    a_vec = Array{Float64}(undef,N_Panel) ; 
    c_vec = Array{Float64}(undef,N_Panel) ; 
    ϵ_vec = Array{Int64}(undef,N_Panel)   ; 
    ζ_vec = Array{Int64}(undef,N_Panel)  ; 

    ## PDFs for ϵ and ζ
    Γ_ϵ = MP_ϵ.Π ; 
    Γ_ζ = MP_ζ.Π ; 

    ## Censor savings 
    # G_ap      .= max.( min.(G_ap,a_max)      , a_min ) ; 
    # G_ap_fine .= max.( min.(G_ap_fine,a_max) , a_min ) ; 

    ## Draw initial conditions from stationary distribution (cross-section)
        println(" Initializing Simulation")
        # Draw ϵ[i] and ζ[i] from (stationary) marginal distribution 
        ϵ_vec .= rand(Categorical(M_Aiyagari.MP_ϵ.PDF),N_Panel) ; 
        ζ_vec .= rand(Categorical(M_Aiyagari.MP_ζ.PDF),N_Panel) ; 
        # Draw a[i] 
        if Seed_Flag # Use conditional distribution in Γ from histogram method 
            for i = 1:N_Panel
            Γ_aux     = Γ[:,ϵ_vec[i],ζ_vec[i]]./sum(Γ[:,ϵ_vec[i],ζ_vec[i]])    ; 
            a_vec[i]  = a_grid_fine[ rand(Categorical(Γ_aux)) ]                ;     
            end 
        else # Use a uniform between 0 and some upper limit (same as in Model) 
            a_vec    .= rand(Uniform(a_min,a_grid_fine[M.n_cut_fine]),N_Panel) ; 
        end 
        println("   E[a_0]=\$$(round(mean(a_vec),digits=2))k  E[ϵ_0]=$(round(mean(M.ϵ_grid[ϵ_vec]),digits=3))  E[ζ_0]=$(round(mean(M.ζ_grid[ζ_vec]),digits=3)) // max(a_0)=$(round(maximum(a_vec),digits=2))")

    ## Iterate forward 
    println(" Iterating Panel")
    for t = 1:T_Simul
        if mod(t,50)==0
            println("   Simulation Period $t")
            println("       E[a_t]=\$$(round(mean(a_vec),digits=2))k  E[ϵ_t]=$(round(mean(M.ϵ_grid[ϵ_vec]),digits=3))  E[ζ_t]=$(round(mean(M.ζ_grid[ζ_vec]),digits=3)) // max(a_t)=$(round(maximum(a_vec),digits=2))")
            # println("       Γ_ϵ = $(round.(100*[sum(ϵ_vec.==j)/N_Panel for j in 1:15],digits=3))")
            # println("       Γ_ζ = $(round.(100*[sum(ζ_vec.==j)/N_Panel for j in 1:7],digits=3))")
        end
    # Simulate each dinsaty            
    for i = 1:N_Panel 
        
        # Compute future assets: Linear interpolation (manually)
        if a_min<a_vec[i]<a_max 
            # if (a_vec[i]>100000)&(mod(t,50)==0)
            #     println("       Asset test - a_vec=$(a_vec[i]) ϵ=$(ϵ_vec[i]) ζ=$(ζ_vec[i]) - ")
            # end 
            i_lo      = min(n_a_fine-1, Grid_Inv(a_vec[i],n_a_fine,θ_a_f,a_min,a_max) ) ;
            ω         = 1 - min(1,max(0,(a_vec[i]-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
            a_vec[i]  = ω*G_ap_fine[i_lo,ϵ_vec[i],ζ_vec[i]] + (1-ω)*G_ap_fine[i_lo+1,ϵ_vec[i],ζ_vec[i]]    ;
            # # G_ap_ip   = ScaledInterpolations( a_grid , G_ap[:,ϵ_vec[i],ζ_vec[i]] , BSpline(Cubic(Line(OnGrid())))) ;
            # G_ap_ip   = ScaledInterpolations( a_grid_fine , G_ap_fine[:,ϵ_vec[i],ζ_vec[i]] , BSpline(Linear()) ) ;
            # a_vec[i]  = G_ap_ip( a_vec[i] ) ;
            # if (a_vec[i]>100000)&(mod(t,50)==0)
            #     println("                   ind=$i_lo, [$(a_grid_fine[i_lo]),$(a_grid_fine[i_lo+1])],  ω=$ω, G(lo)=[$(G_ap_fine[i_lo,ϵ_vec[i],ζ_vec[i]]),$(G_ap_fine[i_lo+1,ϵ_vec[i],ζ_vec[i]])], a'=$(a_vec[i])")
            # end 
        elseif  a_vec[i]==a_max 
            a_vec[i]  = G_ap[end,ϵ_vec[i],ζ_vec[i]] ;
        elseif  a_vec[i]==a_min 
            a_vec[i]  = G_ap[ 1 ,ϵ_vec[i],ζ_vec[i]] ;
        else 
            error("Error in simulation: assets not working")
        end 

        # Censor assets to be between a_min and a_max
        a_vec[i] = min(a_max,max(a_min,a_vec[i])) ;
    
        # Compute future ϵ[i] and ζ[i]
        ϵ_vec[i]  = rand(Categorical( Γ_ϵ[ϵ_vec[i],:] )) ;
        ζ_vec[i]  = rand(Categorical( Γ_ζ[ζ_vec[i],:] )) ; 
        
        # Compute consumption only for relevant periods 
        if t>=T_Simul-(T_Panel-1)
            if a_min<a_vec[i]<a_max 
            G_c_ip    = ScaledInterpolations( a_grid , G_c[:,ϵ_vec[i],ζ_vec[i]] , BSpline(Cubic(Line(OnGrid())))) ; 
            c_vec[i]  = G_c_ip( a_vec[i] ) ;
            elseif  a_min==a_vec[i]
            c_vec[i]  = G_c[1,ϵ_vec[i],ζ_vec[i]] ; 
            else 
            c_vec[i]  = G_c[end,ϵ_vec[i],ζ_vec[i]] ;
            end 
        end
    end
    
    ## Save results in panel  
    if t>=T_Simul-(T_Panel-1)
        M_P.a_mat[:,t-(T_Simul-T_Panel)] .= a_vec ; # Assets 
        M_P.ϵ_mat[:,t-(T_Simul-T_Panel)] .= ϵ_vec ; # Labor efficiency
        M_P.ζ_mat[:,t-(T_Simul-T_Panel)] .= ζ_vec ; # Returns
        M_P.c_mat[:,t-(T_Simul-T_Panel)] .= c_vec ; # Consumption 
    end 

    end 

    return M_P ;
end 





###################################################################
###################################################################
## Simulate the panel - Dynasty by Dynasty 
function Simulate_Panel_Dynasty(M::Model,M_P::Model_Panel)
    ## Unpack relevant objects from model 
    @unpack p, MP_ϵ, MP_ζ, a_grid, a_grid_fine, a_max, Γ, G_ap, G_ap_fine, G_c, n_a_fine, θ_a_f = M
    @unpack a_min = p
    @unpack N_Panel, T_Panel, T_Simul = M_P

    ## Initialize Seed 
    Random.seed!(M_P.rng_seed) ;

    ## PDFs for ϵ and ζ
    Γ_ϵ = MP_ϵ.Π ; 
    Γ_ζ = MP_ζ.Π ; 

    ## Censor savings 
    G_ap      .= max.( min.(G_ap,a_max)      , a_min ) ; 
    G_ap_fine .= max.( min.(G_ap_fine,a_max) , a_min ) ; 

    
    ## Iterate forward 
    println(" Iterating Dynasties")          
    for i = 1:N_Panel 
        aa, M_P.t_vec[i] = @timed begin 
        
        if mod(i,1000)==0
            println("   Dynasty $i")
        end

        ## Draw initial conditions from stationary distribution (cross-section)
        # Draw ϵ[i] and ζ[i] from (stationary) marginal distribution 
        ϵ_vec = rand(Categorical(M_Aiyagari.MP_ϵ.PDF)) ; 
        ζ_vec = rand(Categorical(M_Aiyagari.MP_ζ.PDF)) ; 
        # Draw a[i] from conditional distribution 
        Γ_aux = Γ[:,ϵ_vec,ζ_vec]./sum(Γ[:,ϵ_vec,ζ_vec])    ; 
        a_vec = a_grid_fine[ rand(Categorical(Γ_aux)) ]    ;     
        

    for t = 1:T_Simul
        # if mod(t,100)==0
        #     println("   Simulation Period $t")
        # end
        
        # Compute future assets: Linear interpolation (manually)
        if a_min<a_vec<a_max 
            i_lo   = min(n_a_fine-1, Grid_Inv(a_vec,n_a_fine,θ_a_f,a_min,a_max) ) ;
            ω      = min(1,max(0,(a_vec-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
            a_vec  = ω*G_ap_fine[i_lo,ϵ_vec,ζ_vec] + (1-ω)*G_ap_fine[i_lo+1,ϵ_vec,ζ_vec]    ;
        elseif  a_vec==a_max 
            a_vec  = G_ap[end,ϵ_vec,ζ_vec] ;
        elseif  a_vec==a_min 
            a_vec  = G_ap[ 1 ,ϵ_vec,ζ_vec] ;
        else 
            error("Error in simulation: assets not working")
        end 
        a_vec = min(a_max,max(a_min,a_vec)) ;
    
        # Compute future ϵ[i] and ζ[i]
        ϵ_vec  = rand(Categorical( Γ_ϵ[ϵ_vec,:] )) ;
        ζ_vec  = rand(Categorical( Γ_ζ[ζ_vec,:] )) ; 
        
        # Compute consumption only for relevant periods 
        if t>=T_Simul-(T_Panel-1)
            if a_min<a_vec<a_max 
            G_c_ip = ScaledInterpolations( a_grid , G_c[:,ϵ_vec,ζ_vec] , BSpline(Cubic(Line(OnGrid())))) ; 
            c_vec  = G_c_ip( a_vec ) ;
            elseif  a_min==a_vec
            c_vec  = G_c[1,ϵ_vec,ζ_vec] ; 
            else 
            c_vec  = G_c[end,ϵ_vec,ζ_vec] ;
            end 
        end

        ## Save results in panel  
        if t>=T_Simul-(T_Panel-1)
            M_P.a_mat[i,t-(T_Simul-T_Panel)] = a_vec ; # Assets 
            M_P.ϵ_mat[i,t-(T_Simul-T_Panel)] = ϵ_vec ; # Labor efficiency
            M_P.ζ_mat[i,t-(T_Simul-T_Panel)] = ζ_vec ; # Returns
            M_P.c_mat[i,t-(T_Simul-T_Panel)] = c_vec ; # Consumption 
        end 
    end

    # Save Time 
    end 
    end 

    return M_P ;
end 