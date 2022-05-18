###################################################################
###################################################################
###################################################################
## Functions for simulation 


###################################################################
###################################################################
## Generate structure for the panel using Parameters module
@with_kw struct Model_Panel
    # Size 
    N_Panel::Int64 = 500000 ; # Number of dynasties to be simulated 
    T_Panel::Int64 = 10      ; # Number of periods to be saved of the dynasties 
    T_Simul::Int64 = 2000    ; # Number of periods to be simulated 
    N_Min::Int64   = 1000    ; # Minimum numberof dynasties for moments
    # Panel output 
    a_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Assets
    c_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Consumption 
    ϵ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Labor efficiency index
    h_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Age
    # Random number seed 
    rng_seed::Int64 = 3489398  ;
end





###################################################################
###################################################################
## Simulate the panel 
function Simulate_Panel(M::Model,M_P::Model_Panel)
    ## Unpack relevant objects from model 
    @unpack p, MP_ϵ, a_grid, a_grid_fine, a_max, n_ϵ, Γ, G_ap_fine, G_c_fine, n_a_fine, θ_a_f = M
    @unpack a_min, Surv_Pr, Max_Age = p
    @unpack N_Panel, T_Panel, T_Simul = M_P

    ## Initialize Seed 
    Random.seed!(M_P.rng_seed) ;

    ## Initialize cross-section vectors
    a_vec = Array{Float64}(undef,N_Panel) ; 
    c_vec = Array{Float64}(undef,N_Panel) ; 
    ϵ_vec = Array{Int64}(undef,N_Panel)   ; 
    h_vec = Array{Int64}(undef,N_Panel)   ; 

    ## PDFs for ϵ and ζ
    Γ_ϵ = MP_ϵ.Π ; 

    # Set median ϵ for newborns:
    med_ϵ = convert(Int64,round(n_ϵ/2));

    # Set initial wealth to (close to) $1k 
    b_ind = collect(1:n_a_fine)[a_grid_fine.>=1][1] ; 
    b     = a_grid_fine[b_ind]                      ;

    ## Censor savings 
    G_ap_fine .= max.( min.(G_ap_fine,a_max) , a_min ) ; 

    ## Draw initial conditions from stationary distribution (cross-section)
        println(" Initializing Simulation")
        # Draw ϵ[i] and h[i] from (stationary) marginal distribution 
        ϵ_vec .= rand(Categorical( sum(Γ,dims=(1,3))[:] ),N_Panel) ; 
        h_vec .= rand(Categorical( sum(Γ,dims=(1,2))[:] ),N_Panel) ; 
        # Draw a[i] from conditional distribution 
        for i = 1:N_Panel
        if h_vec[i]==1 
            a_vec[i] = b ; 
        else             
            a_vec[i]  = a_grid_fine[ rand( Categorical( Γ[:,ϵ_vec[i],h_vec[i]]/sum(Γ[:,ϵ_vec[i],h_vec[i]]) ) ) ]            ;     
        end 
        end
        println("   E[a_0]=\$$(round(mean(a_vec),digits=2))k  E[ϵ_0]=$(round(mean(M.ϵ_grid[ϵ_vec]),digits=3))  E[h_0]=$(round(mean(h_vec),digits=1)) // max(a_0)=$(round(maximum(a_vec),digits=2))")

    ## Iterate forward 
    println(" Iterating Panel")
    for t = 1:T_Simul
        if mod(t,50)==0
            println("   Simulation Period $t")
            println("       E[a_t]=\$$(round(mean(a_vec),digits=2))k  E[ϵ_t]=$(round(mean(M.ϵ_grid[ϵ_vec]),digits=3))  E[h_t]=$(round(mean(h_vec),digits=1)) // max(a_t)=$(round(maximum(a_vec),digits=2))")
        end

    # Simulate each dinsaty            
    for i = 1:N_Panel 

        Death_Draw = rand(Categorical( [ Surv_Pr[h_vec[i]] , 1-Surv_Pr[h_vec[i]] ] )) ; 
        
        if Death_Draw==2 # Agent dies

            a_vec[i] = b     ; 
            ϵ_vec[i] = med_ϵ ;
            h_vec[i] = 1     ;

        else # Agent lives 

            # Compute future assets: Linear interpolation (manually)
            if a_min<a_vec[i]<a_max 
                i_lo      = min(n_a_fine-1, Grid_Inv(a_vec[i],n_a_fine,θ_a_f,a_min,a_max) ) ;
                ω         = min(1,max(0,(a_vec[i]-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
                a_vec[i]  = ω*G_ap_fine[i_lo,ϵ_vec[i],h_vec[i]] + (1-ω)*G_ap_fine[i_lo+1,ϵ_vec[i],h_vec[i]]    ;
            elseif  a_vec[i]==a_max 
                a_vec[i]  = G_ap_fine[end,ϵ_vec[i],h_vec[i]] ;
            elseif  a_vec[i]==a_min 
                a_vec[i]  = G_ap_fine[ 1 ,ϵ_vec[i],h_vec[i]] ;
            else 
                error("Error in simulation: assets not working")
            end 
            a_vec[i] = min(a_max,max(a_min,a_vec[i])) ;
        
            # Compute future ϵ[i]
            ϵ_vec[i] = rand(Categorical( Γ_ϵ[ϵ_vec[i],:] )) ;
            
            # Update age 
            h_vec[i] = h_vec[i] + 1  ;

        end 
        
        # Compute consumption only for relevant periods 
        if t>=T_Simul-(T_Panel-1)
            if a_min<a_vec[i]<a_max 
            i_lo      = min(n_a_fine-1, Grid_Inv(a_vec[i],n_a_fine,θ_a_f,a_min,a_max) ) ;
            ω         = min(1,max(0,(a_vec[i]-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
            c_vec[i]  = ω*G_c_fine[i_lo,ϵ_vec[i],h_vec[i]] + (1-ω)*G_c_fine[i_lo+1,ϵ_vec[i],h_vec[i]] ;
            elseif  a_min==a_vec[i]
            c_vec[i]  = G_c_fine[1,ϵ_vec[i],h_vec[i]] ; 
            else 
            c_vec[i]  = G_c_fine[end,ϵ_vec[i],h_vec[i]] ;
            end 
        end
    end
    
    ## Save results in panel  
    if t>=T_Simul-(T_Panel-1)
        M_P.a_mat[:,t-(T_Simul-T_Panel)] .= a_vec ; # Assets 
        M_P.ϵ_mat[:,t-(T_Simul-T_Panel)] .= ϵ_vec ; # Labor efficiency
        M_P.h_mat[:,t-(T_Simul-T_Panel)] .= h_vec ; # Age
        M_P.c_mat[:,t-(T_Simul-T_Panel)] .= c_vec ; # Consumption 
    end 

    end 

    return M_P ;
end 









###################################################################
###################################################################
## Generate structure for the panel using Parameters module
@with_kw struct Model_Cohort
    # Size 
    N_Panel::Int64 = 500000  ; # Number of individuals to be simulated 
    T_Panel::Int64 = 100     ; # Provisional number of periods 
    # Panel output 
    a_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Assets
    c_mat = Array{Float32}(undef,N_Panel,T_Panel) ; # Consumption 
    ϵ_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Labor efficiency index
    h_mat = Array{Int32}(undef,N_Panel,T_Panel)   ; # Alive indicator (1->Alive, 0->Dead)
    # Random number seed 
    rng_seed::Int64 = 297835398  ;
end





###################################################################
###################################################################
## Simulate the cohort  
function Simulate_Cohort(M::Model,M_C::Model_Cohort,age_0::Int)
    ## Unpack relevant objects from model 
    @unpack p, MP_ϵ, a_grid, a_grid_fine, a_max, n_ϵ, Γ, G_ap_fine, G_c_fine, n_a_fine, θ_a_f = M
    @unpack a_min, Surv_Pr, Max_Age = p
    @unpack N_Panel = M_C

    ## Max_Cohort_Age 
    Max_C_Age = Max_Age - (age_0-1) ; 

    ## Initialize Seed 
    Random.seed!(M_C.rng_seed) ;

    ## Initialize cross-section vectors
    a_mat = zeros(N_Panel,Max_C_Age) ; 
    c_mat = zeros(N_Panel,Max_C_Age) ; 
    ϵ_mat = Array{Int64}(undef,N_Panel,Max_C_Age) ; 
    h_mat = Array{Int64}(undef,N_Panel,Max_C_Age) ; 
    
    ## PDFs for ϵ and ζ
    Γ_ϵ = MP_ϵ.Π ; 
    
    
    # Set median ϵ for newborns:
    med_ϵ = convert(Int64,round(n_ϵ/2));

    # Set initial wealth to (close to) $1k 
    b_ind = collect(1:n_a_fine)[a_grid_fine.>=1][1] ; 
    b     = a_grid_fine[b_ind]                      ;

    ## Censor savings 
    G_ap_fine .= max.( min.(G_ap_fine,a_max) , a_min ) ; 

    ## Draw initial conditions from stationary distribution (cross-section)
        println(" Initializing Simulation")
        # Draw ϵ[i] and ζ[i] from (stationary) marginal distribution 
        ϵ_mat[:,1] .= rand(Categorical( sum(Γ[:,:,age_0]/sum(Γ[:,:,age_0]),dims=(1,3))[:] ),N_Panel) ; 
        h_mat[:,1] .= 1 ; 
        # Draw a[i] from conditional distribution 
        if age_0==1 
                a_mat[:,1] .= b ; 
                c_mat[:,1] .= G_c_fine[b_ind,med_ϵ,1] ;
        else 
            for i = 1:N_Panel
                Γ_aux       = Γ[:,ϵ_mat[i,1],age_0]./sum(Γ[:,ϵ_mat[i,1],age_0])   ; 
                a_ind       = rand(Categorical(Γ_aux))              ;
                a_mat[i,1]  = a_grid_fine[ a_ind ]                  ;     
                c_mat[i,1]  = G_c_fine[a_ind,ϵ_mat[i,1],age_0]      ;
            end 
        end 
        println("   E[a_0]=\$$(round(mean(a_mat[:,1]),digits=2))k  E[ϵ_0]=$(round(mean(M.ϵ_grid[ϵ_mat[:,1]]),digits=3))  E[h_0]=$(round(mean(h_mat[:,1]),digits=1)) // max(a_0)=$(round(maximum(a_mat[:,1]),digits=2))")

    ## Iterate forward 
    println(" Iterating Panel")
    for t = 2:Max_C_Age

        if mod(t,5)==0
            ind = h_mat[:,t-1].>0 ; # Select only alive agents
            println("   Simulation Period $t")
            println("       E[a_t]=\$$(round(mean(a_mat[ind,t-1]),digits=2))k  E[ϵ_t]=$(round(mean(M.ϵ_grid[ϵ_mat[ind,t-1]]),digits=3))  E[h_t]=$(round(mean(h_mat[:,t-1]),digits=1)) // max(a_t)=$(round(maximum(a_mat[ind,t-1]),digits=2))")
        end

        # Set current age of agents in cohort 
        age = t+(age_0-1) ; 
        
    # Simulate each dinsaty            
    for i = 1:N_Panel 

        # Check if alive // if not just don't update and leave entries as zeros
        if h_mat[i,t-1]==1 

        Death_Draw = rand(Categorical( [ Surv_Pr[age] , 1-Surv_Pr[age] ] )) ; 
        
        if Death_Draw==2 # Agent dies, no need to update 

            h_mat[i,t] = 0 ;

        else # Agent lives 

            # Compute future assets: Linear interpolation (manually)
            if a_min<a_mat[i,t-1]<a_max 
                i_lo        = min(n_a_fine-1, Grid_Inv(a_mat[i,t-1],n_a_fine,θ_a_f,a_min,a_max) ) ;
                ω           = min(1,max(0,(a_mat[i,t-1]-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
                a_mat[i,t]  = ω*G_ap_fine[i_lo,ϵ_mat[i,t-1],age-1] + (1-ω)*G_ap_fine[i_lo+1,ϵ_mat[i,t-1],age-1]      ;
            elseif  a_mat[i,t-1]==a_max 
                a_mat[i,t]  = G_ap_fine[end,ϵ_mat[i,t-1],age-1] ;
            elseif  a_mat[i,t-1]==a_min 
                a_mat[i,t]  = G_ap_fine[ 1 ,ϵ_mat[i,t-1],age-1] ;
            else 
                error("Error in simulation: assets not working")
            end 
            a_mat[i,t] = min(a_max,max(a_min,a_mat[i,t])) ;
        
            # Compute future ϵ[i]
            ϵ_mat[i,t] = rand(Categorical( Γ_ϵ[ϵ_mat[i,t-1],:] )) ;
            
            # Compute consumption 
            if a_min<a_mat[i,t]<a_max 
                i_lo        = min(n_a_fine-1, Grid_Inv(a_mat[i,t],n_a_fine,θ_a_f,a_min,a_max) ) ;
                ω           = min(1,max(0,(a_mat[i,t]-a_grid_fine[i_lo])/(a_grid_fine[i_lo+1]-a_grid_fine[i_lo]))) ; 
                c_mat[i,t]  = ω*G_c_fine[i_lo,ϵ_mat[i,t],age] + (1-ω)*G_c_fine[i_lo+1,ϵ_mat[i,t],age] ;
            elseif  a_min==a_mat[i,t]
                c_mat[i,t]  = G_c_fine[1,ϵ_mat[i,t],age] ; 
            else 
                c_mat[i,t]  = G_c_fine[end,ϵ_mat[i,t],age] ;
            end 

            # Update age 
            h_mat[i,t] = 1  ;

        end 
        end 
        
    end
    end 
 

    M_C = Model_Cohort(M_C; a_mat=a_mat, c_mat=c_mat, ϵ_mat=ϵ_mat, h_mat=h_mat) ;
    return M_C
end 