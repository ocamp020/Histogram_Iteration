#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(c,p::Par)
    if p.γ>1
    return (c).^(1-p.γ)/(1-p.γ)
    else
    return log.(c)
    end
end

function d_utility(c,p::Par)
    return (c).^(-p.γ)
end

function d_utility_inv(x,p::Par)
    return x.^(-1/p.γ)
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Warm Glow function
function warm_glow(ap,p::Par)
    if p.γ_b>1
    return p.χ*(ap).^(1-p.γ_b)/(1-p.γ_b)
    else
    return p.χ*log.(c)
    end
end

function d_warm_glow(ap,p::Par)
    return p.χ*(ap).^(-p.γ_b)
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Policy function: PFI Fixed Point 
function Solve_Policy_Functions(T::Function,M::Model)
    # Unpack model structure
    @unpack p, n_ϵ, n_a, n_a_fine, θ_a, a_grid, a_grid_fine,method = M
    # PFI paramters
    @unpack max_iter, dist_tol, dist_tol_Δ, r, w, η, Max_Age = p
    
    println(" ")
    println("------------------------")
    println("Backwards Induction - n_ϵ=$n_ϵ, Max_Age=$Max_Age, n_a=$n_a - θ_a=$θ_a - r=$r")
        # Update value function
        G_ap, G_c = T(M)
         
        # Check borrowing constraint
        if any( G_ap.<p.a_min )
            error("Borrowing Constraint Violated")
        end

        # Interpolate to fine grid
        G_ap_fine = zeros(n_a_fine,n_ϵ,Max_Age)
        G_c_fine  = zeros(n_a_fine,n_ϵ,Max_Age)
        for age=1:Max_Age
        for i_ϵ=1:n_ϵ
        G_ap_ip = ScaledInterpolations(a_grid,G_ap[:,i_ϵ,age] , BSpline(Cubic(Line(OnGrid()))))
            G_ap_fine[:,i_ϵ,age].= G_ap_ip.(collect(a_grid_fine))
        G_c_ip  = ScaledInterpolations(a_grid,G_c[:,i_ϵ,age]  , BSpline(Cubic(Line(OnGrid()))))
            G_c_fine[:,i_ϵ,age] .= G_c_ip.(collect(a_grid_fine))
        end
        end
        
        # Update model
        M = Model(M; G_ap=G_ap,G_c=G_c,G_ap_fine=G_ap_fine,G_c_fine=G_c_fine)

    println("\n------------------------\n")
    
    return M

end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - EGM - Backwards Induction
function T_BI_G(M::Model)
    @unpack p, n_ϵ, MP_ϵ, n_a, y_mat, a_grid, a_mat_aϵ, G_ap, method = M
    @unpack β, a_min, r, w, Max_Age, Surv_Pr, Age_Π = p

    ## Allocate policy functioins 
    G_c  = Array{Float64}(undef,n_a,n_ϵ,Max_Age)
    G_ap = Array{Float64}(undef,n_a,n_ϵ,Max_Age)

    ## Warm Glow vector and matrix 
        warm_glow_vec = d_warm_glow(a_grid,p)       ;
        warm_glow_mat = repeat(warm_glow_vec,1,n_ϵ) ;

    ## Final Period 
        if p.χ>0 # Check for warm glow bequest motive 

        # Endogenous consumption for all levels of bequests (ap)
        C_endo = d_utility_inv( β*warm_glow_vec , p)
        # Define endogenous grid on assets
        A_endo = ( repeat( C_endo.+ a_grid ,1,n_ϵ) - y_mat[:,:,Max_Age] )./(1+r)
        # Interpolate functions on exogenous grid
        for i_ϵ=1:n_ϵ
            # Sort A_endo for interpolation
            sort_ind = sortperm(A_endo[:,i_ϵ])
            A_aux    = A_endo[:,i_ϵ][sort_ind]
            C_aux    = C_endo[sort_ind]
            # Check boundary condition
            if minimum(A_aux)>a_min
                a_vec = M.a_grid[a_grid.<minimum(A_aux)]
                A_aux = [a_vec ; A_aux]
                C_aux = [((1+r).*a_vec.+y_mat[1,i_ϵ,Max_Age].-a_min) ; C_aux]
            end
            C_ip                 = Spline1D(A_aux,C_aux)
            G_c[:,i_ϵ,Max_Age]  .= C_ip.(M.a_grid)
            G_ap[:,i_ϵ,Max_Age] .= (1+r).*collect(a_grid) .+ y_mat[:,i_ϵ,Max_Age] .- G_c[:,i_ϵ,Max_Age]
        end

        else 
            
            # println("   age = $Max_Age - No bequest motive, setting savings to zero")
            G_c[:,:,Max_Age]  .= (1+r).*a_mat_aϵ .+ y_mat[:,:,Max_Age] .- a_min
            G_ap[:,:,Max_Age] .= a_min 

        end 


    ## Backwards Induction Loop 
    for age in (Max_Age-1):-1:1
        # Define RHS of Euler equation for each (a',ϵ,age)
            # The matrix d_utility is (a',ϵ') for a given age, the transition matrix is transposed so that it gives (ϵ',ϵ), result is (a,ϵ)
        Euler_RHS = β*( Surv_Pr[age]*(1+r).* d_utility( (1+r).*a_mat_aϵ + y_mat[:,:,age+1] - G_ap[:,:,age+1] , p )*(MP_ϵ.Π)'  .+ (1 .- Surv_Pr[age]).*warm_glow_mat  )
        # Check Monotonicity
        if any( Euler_RHS.<0 )
            error("RHS must be monotone for EGM to work")
        end
        # Endogenous consumption for all levels of bequests (ap)
        C_endo = d_utility_inv( Euler_RHS , p)
        # Define endogenous grid on assets
        A_endo = ( C_endo .+ a_mat_aϵ - y_mat[:,:,age] )./(1+r)
        # Interpolate functions on exogenous grid
        for i_ϵ=1:n_ϵ
            # Sort A_endo for interpolation
            sort_ind = sortperm(A_endo[:,i_ϵ])
            A_aux    = A_endo[:,i_ϵ][sort_ind]
            C_aux    = C_endo[:,i_ϵ][sort_ind]
            # Check boundary condition
            if minimum(A_aux)>a_min
                a_vec = M.a_grid[a_grid.<minimum(A_aux)]
                A_aux = [a_vec ; A_aux]
                C_aux = [((1+r).*a_vec.+y_mat[1,i_ϵ,age].-a_min) ; C_aux]
            end
            C_ip             = Spline1D(A_aux,C_aux)
            G_c[:,i_ϵ,age]  .= C_ip.(M.a_grid)
            G_ap[:,i_ϵ,age] .= (1+r).*collect(a_grid) .+ y_mat[:,i_ϵ,age] .- G_c[:,i_ϵ,age]
        end
    end 
    
    
    # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (1+r)*M.a_mat[ind] + y_mat[ind] - a_min
        end
        # # Check for borrowing constraint
        # if any( G_ap.<a_min )
        #    error("Borrowing Constraint Violated")
        # end
    
    
    # Return Results
    return G_ap, G_c
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method
function Histogram_Method_Loop(M::Model,Γ_0=nothing)
    @unpack p, n_ϵ, MP_ϵ, n_a_fine, a_grid_fine, θ_a, θ_a_f, G_ap_fine = M
    @unpack a_min, Hist_max_iter, Hist_tol, Hist_η, Max_Age, Surv_Pr = p

    println("\n--------------------------------\nBegining Histogram Method with Loops")
    
    # Initial distribution
    if Γ_0==nothing
        Γ_0 = M.Γ
    end

    # Set median ϵ for newborns:
    med_ϵ = convert(Int64,round(n_ϵ/2));

    # Discretize distribution
    println("   Discretizing Choices and Computing Transition Probabilities")
    H_ind    = Array{Int64}(undef,n_a_fine,n_ϵ,Max_Age)
    H_ω_lo_s   = Array{Float64}(undef,n_a_fine,n_ϵ,Max_Age,n_ϵ)
    H_ω_hi_s   = Array{Float64}(undef,n_a_fine,n_ϵ,Max_Age,n_ϵ)
    H_ω_lo_d   = Array{Float64}(undef,n_a_fine,n_ϵ,Max_Age)
    H_ω_hi_d   = Array{Float64}(undef,n_a_fine,n_ϵ,Max_Age)
    for age=1:Max_Age # Current age
        Pr_age = Surv_Pr[age] ; # Transitions of future age conditional on current age
    for i_ϵ=1:n_ϵ # Current ϵ
        Pr_ϵp = MP_ϵ.Π[i_ϵ,:] ; # Transitions of future ϵ conditional on current ϵ
    for i_a=1:n_a_fine # Current a

        # Get index and weight of lower bound on approximation interval
        if     ( G_ap_fine[i_a,i_ϵ,age] ≥ M.a_max)
            H_ind[i_a,i_ϵ,age] = n_a_fine - 1 ;
            ω                  = 0            ;
        elseif ( G_ap_fine[i_a,i_ϵ,age] ≤ a_min)
            H_ind[i_a,i_ϵ,age] = 1 ; 
            ω                  = 1 ;
        else
            H_ind[i_a,i_ϵ,age] = Grid_Inv(G_ap_fine[i_a,i_ϵ,age],n_a_fine,θ_a_f,a_min,M.a_max) ;
            ω                  = min(1,max(0,(G_ap_fine[i_a,i_ϵ,age]-a_grid_fine[H_ind[i_a,i_ϵ,age]])/(a_grid_fine[H_ind[i_a,i_ϵ,age]+1]-a_grid_fine[H_ind[i_a,i_ϵ,age]]))) ;
        end
        
        # Transition matrices for survival 
        for i_ϵp=1:n_ϵ
        H_ω_lo_s[i_a,i_ϵ,age,i_ϵp] = Pr_age*Pr_ϵp[i_ϵp]*ω     ;
        H_ω_hi_s[i_a,i_ϵ,age,i_ϵp] = Pr_age*Pr_ϵp[i_ϵp]*(1-ω) ;
        end 

        # Transition matrices for death  
        H_ω_lo_d[i_a,i_ϵ,age] = (1-Pr_age)*ω     ;
        H_ω_hi_d[i_a,i_ϵ,age] = (1-Pr_age)*(1-ω) ;

        # Check for conservation of mass
        H_res = sum( H_ω_lo_s[i_a,i_ϵ,age,:]+H_ω_hi_s[i_a,i_ϵ,age,:] ) + H_ω_lo_d[i_a,i_ϵ,age] 
        if abs( (1-H_res) - H_ω_hi_d[i_a,i_ϵ,age])<1e-5 
            H_ω_hi_d[i_a,i_ϵ,age] = 1 - H_res 
        else 
            println("Error in transition probabilities: i_a=$i_a, i_ϵ=$i_ϵ, age=$age, residual=$(abs( (1-H_res) - H_ω_hi_d[i_a,i_ϵ,age]))")
        end 
    end
    end
    end
        # Check transitions 
        H_aux = sum(H_ω_lo_s+H_ω_hi_s,dims=4)+H_ω_lo_d+H_ω_hi_d ;
        @printf("\n   Check transition functions: maximum = %.4e, minimum = %.4e ", maximum(H_aux)-1 , minimum(H_aux)-1  )
        H_ω_hi_d = 1 .- (sum(H_ω_lo_s+H_ω_hi_s,dims=4)+H_ω_lo_d)
        @printf("\n   Check initial distribution: sum(Γ_0) = %.4e \n", sum(Γ_0)  )


    # Loop for updating histogram
    println("\n     Iterating on the Distribution")
    H_dist = 1
    for i_H=1:Hist_max_iter 
        # Update histogram
        Γ = zeros(n_a_fine,n_ϵ,Max_Age)

        # Update until second to last period of life
        for age=1:Max_Age-1  # Current age
        for i_ϵ=1:n_ϵ        # Current ϵ
        for i_a=1:n_a_fine   # Current a

            i_ap = H_ind[i_a,i_ϵ,age]    ;
            
            # If agents survive
            for i_ϵp=1:n_ϵ # Future ϵ
                Γ[i_ap  ,i_ϵp,age+1] = Γ[i_ap  ,i_ϵp,age+1] + H_ω_lo_s[i_a,i_ϵ,age,i_ϵp]*Γ_0[i_a,i_ϵ,age]
                Γ[i_ap+1,i_ϵp,age+1] = Γ[i_ap+1,i_ϵp,age+1] + H_ω_hi_s[i_a,i_ϵ,age,i_ϵp]*Γ_0[i_a,i_ϵ,age]
                # Γ[i_ap  ,i_ϵp,age+1] = Γ[i_ap  ,i_ϵp,age+1] + Surv_Pr[age]*Pr_ϵp[i_ϵp]*(  H_ω[i_a,i_ϵ,age])*Γ_0[i_a,i_ϵ,age]
                # Γ[i_ap+1,i_ϵp,age+1] = Γ[i_ap+1,i_ϵp,age+1] + Surv_Pr[age]*Pr_ϵp[i_ϵp]*(1-H_ω[i_a,i_ϵ,age])*Γ_0[i_a,i_ϵ,age]
            end
            
            # If agents die: age=1 and ϵ=median(ϵ)
                Γ[i_ap  ,med_ϵ,1] = Γ[i_ap  ,med_ϵ,1] + H_ω_lo_d[i_a,i_ϵ,age]*Γ_0[i_a,i_ϵ,age]
                Γ[i_ap+1,med_ϵ,1] = Γ[i_ap+1,med_ϵ,1] + H_ω_hi_d[i_a,i_ϵ,age]*Γ_0[i_a,i_ϵ,age]
                # Γ[i_ap  ,med_ϵ,1] = Γ[i_ap  ,med_ϵ,1] + (1-Surv_Pr[age])*(  H_ω[i_a,i_ϵ,age])*Γ_0[i_a,i_ϵ,age]
                # Γ[i_ap+1,med_ϵ,1] = Γ[i_ap+1,med_ϵ,1] + (1-Surv_Pr[age])*(1-H_ω[i_a,i_ϵ,age])*Γ_0[i_a,i_ϵ,age]
        end
        end
        end

        # Final period of life (all agents die)
        for i_ϵ=1:n_ϵ      # Current ϵ
        for i_a=1:n_a_fine # Current a
            
            i_ap = H_ind[i_a,i_ϵ,Max_Age]    ;

            Γ[i_ap  ,med_ϵ,1] = Γ[i_ap  ,med_ϵ,1] + H_ω_lo_d[i_a,i_ϵ,Max_Age]*Γ_0[i_a,i_ϵ,Max_Age]
            Γ[i_ap+1,med_ϵ,1] = Γ[i_ap+1,med_ϵ,1] + H_ω_hi_d[i_a,i_ϵ,Max_Age]*Γ_0[i_a,i_ϵ,Max_Age]
            # Γ[i_ap  ,med_ϵ,1] = Γ[i_ap  ,med_ϵ,1] + (  H_ω[i_a,i_ϵ,Max_Age])*Γ_0[i_a,i_ϵ,Max_Age]
            # Γ[i_ap+1,med_ϵ,1] = Γ[i_ap+1,med_ϵ,1] + (1-H_ω[i_a,i_ϵ,Max_Age])*Γ_0[i_a,i_ϵ,Max_Age]

        end
        end

        # Update distance
        H_dist = maximum(abs.(Γ-Γ_0))
        
        # Update initial distribution
        Γ_0 .= (1-Hist_η)*Γ .+ Hist_η*Γ_0
        
        # Report progress
        if mod(i_H,10)==0
            @printf("\n         Histogram Loop: iter = %d, dist = %.4e, check = %.3f",i_H,H_dist,sum(Γ))
        end
        # Check convergence
        if H_dist<Hist_tol
            @printf("\n     Histogram iteartion converged in iteration %d. H_dist=%.4e \n--------------------------------\n",i_H,H_dist)
            M = Model(M; Γ=Γ, H_ind=H_ind, H_ω_lo_s=H_ω_lo_s, H_ω_hi_s=H_ω_hi_s, H_ω_lo_d=H_ω_lo_d, H_ω_hi_d=H_ω_hi_d)
            return M
        end
    end
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Euler Error Function

# G_ap interpolation
function G_ap_ϵa(age::Int64,i_ϵ::Int64,a,M::Model)
    itp = ScaledInterpolations(M.a_grid,M.G_ap[:,i_ϵ,age], BSpline(Cubic(Line(OnGrid()))))
    return itp(a)
end

# Euler Equation Percentage Error
function Euler_Error(i_ϵ::Int64,age::Int64,a,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_ϵ, n_ϵ, ϵ_grid, y_mat = M
    @unpack β, r, w, Max_Age, Surv_Pr = p

    # Iterpolate G_ap at current ϵ
    ap  = min(M.a_grid[end],max(M.a_grid[1],G_ap_ϵa(age,i_ϵ,a,M)))

    # Current consumption
    c   = (1+r)*a +y_mat[1,i_ϵ,age] - ap

    # Compute left hand side of Euler equation
    LHS = d_utility(c,p)

    # Compute right hand side of Euler equation
    if age<M.p.Max_Age 
    RHS_prob = zeros(n_ϵ)
        # Marginal utility at age+1,ϵ',a' ,G_ap(a',ϵ',age+1)
        for i_ϵp=1:n_ϵ
        cp  = (1+r)*ap + y_mat[1,i_ϵp,age+1] - G_ap_ϵa(age+1,i_ϵp,ap,M)
        up  = d_utility(cp,p)
        RHS_prob[i_ϵp] = β*( Suvr_Pr[age]*MP_ϵ.Π[i_ϵ,i_ϵp]*(1+r)*up + (1-Surv_Pr[age])*d_warm_glow(ap,p) )
        end
    RHS = sum(RHS_prob)
    else 
    RHS = β*( d_warm_glow(ap,p) )
    end 
    # Return percentage errror in Euler equation
    return (RHS/LHS-1)*100
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Aiyagari Equilibrium
function Aiyagari_Equilibrium(M_in::Model)

    # Read files if required 
    if M_in.read_flag
        Γ_0    = reshape( readdlm(File_Folder*"/Distribution.csv"   , ',', Float64) ,M_in.n_a_fine, M.n_ϵ, M.p.Max_Age );
    else 
        Γ_0    = nothing ;
    end 

    # Compute Policy Functions
    M_in   = Solve_Policy_Functions(T_BI_G,M_in) ;

    # Compute Distribution
    M_in  = Histogram_Method_Loop(M_in,Γ_0) ;

    # Save Results 
    open(File_Folder*"/Policy_Function.csv", "w") do io
        writedlm(io, M_in.G_ap , ',')
    end;
    open(File_Folder*"/Distribution.csv", "w") do io
        writedlm(io, M_in.Γ , ',')
    end;

    return M_in
end