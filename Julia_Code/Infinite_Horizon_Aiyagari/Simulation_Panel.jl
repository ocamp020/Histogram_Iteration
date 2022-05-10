###################################################################
###################################################################
###################################################################
## Functions for simulation 


###################################################################
###################################################################
## Generate structure for the panel using Parameters module
@with_kw struct Model_Panel
    # Size 
    N_Panel::Int64 = 100000  ; # Number of dynasties to be simulated 
    T_Panel::Int64 = 10      ; # Number of periods to be saved of the dynasties 
    T_Simul::Int64 = 2000    ; # Number of periods to be simulated 
    N_Min::Int64   = 10000   ; # Minimum numberof dynasties for moments
    # Panel output 
    a_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Assets
    c_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Consumption 
    ϵ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Labor efficiency index
    ζ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Return index
    # Random number seed 
    rng_seed::Int64 = 3489398  ;
end





###################################################################
###################################################################
## Simulate the panel 
function Simulate_Panel(M::Model,M_P::Model_Panel)
    ## Unpack relevant objects from model 
    @unpack p, MP_ϵ, MP_ζ, n_ϵ, n_ζ, ϵ_grid, ζ_grid, a_grid, a_grid_fine, a_max, Γ, G_ap, G_c = M
    @unpack a_min = p
    @unpack N_Panel, T_Panel, T_Simul = M_P

    ## Initialize Seed 
    Random.seed!(M_P.rng_seed) ;

    ## Initialize cross-section vectors
    a_vec = Array{Float64}(undef,N_Panel) ; 
    c_vec = Array{Float64}(undef,N_Panel) ; 
    ϵ_vec = Array{Int64}(undef,N_Panel)   ; 
    ζ_vec = Array{Int64}(undef,N_Panel)  ; 

    ## PDFs for ϵ and ζ
    Γ_ϵ = MP_ϵ.Π
    Γ_ζ = MP_ζ.Π

    ## Censor savings 
    G_ap .= max.(G_ap,a_max)

    ## Draw initial conditions from stationary distribution (cross-section)
        println(" Initializing Simulation")
        # Draw ϵ[i] and ζ[i] from (stationary) marginal distribution 
        ϵ_vec .= rand(Categorical(M_Aiyagari.MP_ϵ.PDF),N_Panel) ; 
        ζ_vec .= rand(Categorical(M_Aiyagari.MP_ζ.PDF),N_Panel) ; 
        # Draw a[i] from conditional distribution 
        for i = 1:N_Panel
        Γ_aux     = Γ[:,ϵ_vec[i],ζ_vec[i]]./sum(Γ[:,ϵ_vec[i],ζ_vec[i]])    ; 
        a_vec[i]  = a_grid_fine[ rand(Categorical(Γ_aux)) ]                ;     
        end 

    ## Iterate forward 
    for t = 1:T_Simul
        if mod(t,25)==0
            println(" Simulation Period $t")
        end
    # Simulate each dinsaty            
    for i = 1:N_Panel 
        
        # Compute future assets
        if a_min<a_vec[i]<a_max 
        G_ap_ip   = ScaledInterpolations( a_grid , G_ap[:,ϵ_vec[i],ζ_vec[i]] , BSpline(Cubic(Line(OnGrid())))) ;
        a_vec[i]  = G_ap_ip.( a_vec[i] ) ;
        elseif  a_min==a_vec[i]
        a_vec[i]  = G_ap[1,ϵ_vec[i],ζ_vec[i]] ; 
        else 
        a_vec[i]  = G_ap[end,ϵ_vec[i],ζ_vec[i]] ;
        end 
    
        # Compute future ϵ[i] and ζ[i]
        ϵ_vec[i]  = rand(Categorical( Γ_ϵ[ϵ_vec[i],:] )) ;
        ζ_vec[i]  = rand(Categorical( Γ_ζ[ζ_vec[i],:] )) ; 
        
        # Compute consumption only for relevant periods 
        if t>=T_Simul-(T_Panel-1)
            if a_min<a_vec[i]<a_max 
            G_c_ip    = ScaledInterpolations( a_grid , G_c[:,ϵ_vec[i],ζ_vec[i]] , BSpline(Cubic(Line(OnGrid())))) ; 
            c_vec[i]  = G_c_ip.( a_vec[i] ) ;
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