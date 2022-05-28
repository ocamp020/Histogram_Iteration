#-----------------------------------------------------------
#-----------------------------------------------------------
# Utility function
function utility(c,p::Par)
    if p.σ>1
    return (c).^(1-p.σ)/(1-p.σ)
    else
    return log.(c)
    end
end

function d_utility(c,p::Par)
    return (c).^(-p.σ)
end

function d_utility_inv(x,p::Par)
    return x.^(-1/p.σ)
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Policy function: PFI Fixed Point 
function PFI_Fixed_Point(T::Function,M::Model,G_ap_old=nothing)
    # Unpack model structure
    @unpack p, n_ϵ, n_ζ, n_a, n_a_fine, θ_a, a_grid, a_grid_fine,method = M
    # PFI paramters
    @unpack max_iter, dist_tol, dist_tol_Δ, r, w, η = p
    # Initialize variables for loop
    if G_ap_old==nothing
        G_ap_old  = (r*M.ζ_mat.+1).*M.a_mat
    end
    G_dist_new    = 100          ; # Initialize distance
    G_dist_old    = 1            ; # Initialize distance
    G_dist_change = 1            ; # Initialize change in distance
    println(" ")
    println("------------------------")
    println("PFI - n_ϵ=$n_ϵ, n_ζ=$n_ζ, n_a=$n_a - θ_a=$θ_a - r=$r")
    for iter=1:max_iter
        # Update old distance and iterations
        G_dist_old = G_dist_new
        # Update value function
        G_ap_new, G_c = T(Model(M,G_ap=copy(G_ap_old)))
        # Update new distance and iterations
        G_dist_new = sqrt(norm(G_ap_new-G_ap_old,2))
        # Update change in distance
        G_dist_change = abs(G_dist_new-G_dist_old)
        # Update old function
        G_ap_old  = (1-η)*G_ap_new .+ η*G_ap_old
        # Update change in convergence criteria
        # Report progress
        if mod(iter,250)==0
            @printf("\n   PFI Loop: iter = %d, dist = %.4e",iter,G_dist_new)
        end
        # Check convergence and return results
        if G_dist_new<=dist_tol ||  ((G_dist_change<=dist_tol_Δ)&&(G_dist_new<=dist_tol*10))
            
            if G_dist_new<=dist_tol
            @printf("\n\nDistance converged: Iterations = %d, Distance = %.4e",iter,G_dist_new)
            elseif ((G_dist_change<=dist_tol_Δ)&&(G_dist_new<=dist_tol*10))
            @printf("\n\nDistance fluctuating: Iterations = %d, Distance = %.4e",iter,G_dist_new)
            @printf("\nChange in distance converged: Iterations = %d, Δ Distance = %.4e",iter,G_dist_change)
            end
            
            println("\n------------------------\n")
            # Check borrowing constraint
            if any( G_ap_new.<p.a_min )
                error("Borrowing Constraint Violated")
            end
            # Interpolate to fine grid
            G_ap_fine = zeros(n_a_fine,n_ϵ,n_ζ)
            G_c_fine  = zeros(n_a_fine,n_ϵ,n_ζ)
            for i_ϵ=1:n_ϵ
            for i_ζ=1:n_ζ
            G_ap_ip = ScaledInterpolations( a_grid , G_ap_new[:,i_ϵ,i_ζ] , BSpline(Cubic(Line(OnGrid()))))
                G_ap_fine[:,i_ϵ,i_ζ].= G_ap_ip.(collect(a_grid_fine))
            G_c_ip  = ScaledInterpolations( a_grid , G_c[:,i_ϵ,i_ζ] , BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[:,i_ϵ,i_ζ] .= G_c_ip.(collect(a_grid_fine))
            end
            end
            # Update model
            M = Model(M; G_ap=G_ap_new,G_c=G_c,G_ap_fine=G_ap_fine,G_c_fine=G_c_fine)
            return M
        end
end
# If loop ends there was no convergence -> Error!
error("Error in PFI - Solution not found")
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Bellman operator - EGM - Iterate on Policy Functions
function T_EGM_G(M::Model)
    @unpack p, n_ζ, n_ϵ, MP_ϵ, MP_ζ, n_a, G_ap, method = M
    @unpack β, a_min, r, w = p
    # Define RHS of Euler equation for each (a',ϵ,ζ)
    # First dim is tomorrow's a in a fixed grid, second dim is present ϵ, third dim is present ζ
    if method == 1
        # Kronecker product of  Π_ζ ⊗ Π_ϵ
        Π_joint = kronecker(MP_ζ.Π, MP_ϵ.Π)
        Euler_RHS = β*Π_joint*(reshape((r*M.ζ_mat.+1).*d_utility( (r*M.ζ_mat.+1).*M.a_mat + w*M.ϵ_mat - G_ap , p ),(n_a,n_ϵ*n_ζ))')
        Euler_RHS = reshape(Euler_RHS',(n_a,n_ϵ,n_ζ))
    elseif method == 2
        # Loop calculation  
        aux_1          = (r*M.ζ_mat.+1).*d_utility( (r*M.ζ_mat.+1).*M.a_mat + w*M.ϵ_mat - G_ap , p );
        aux_2          = Array{Float64}(undef,n_ϵ,n_ζ); 
        Euler_RHS = Array{Float64}(undef,n_a,n_ϵ,n_ζ); 
        for i_ζ=1:n_ζ # Today's ζ
            Pr_ζp = MP_ζ.Π[i_ζ,:] ;
        for i_ϵ=1:n_ϵ # Today's ϵ
            Pr_ϵp = MP_ϵ.Π[i_ϵ,:] ;
        for i_ap=1:n_a # Tomorrow's Assets 
            for i_ζp=1:n_ζ
            for i_ϵp=1:n_ϵ
                aux_2[i_ϵp,i_ζp] = β*aux_1[i_ap,i_ϵp,i_ζp]*Pr_ϵp[i_ϵp]*Pr_ζp[i_ζp]
            end  
            end 
            Euler_RHS[i_ap,i_ϵ,i_ζ] = sum(aux_2[:,:])
        end 
        end 
        end
    end
    #println("Maximum = ",maximum(abs.(Euler_RHS-Euler_RHS_Loop)))
    # Check Monotonicity
    if any( Euler_RHS.<0 )
        error("RHS must be monotone for EGM to work")
    end
    # Define consumption from Euler equation
    C_endo = d_utility_inv(Euler_RHS,p)
    # Define endogenous grid on assets
    A_endo = (C_endo .+ M.a_mat - w*M.ϵ_mat)./(r*M.ζ_mat.+1)
    # Interpolate functions on exogenous grid
    G_c = Array{Float64}(undef,n_a,n_ϵ,n_ζ)
    for i_ζ=1:n_ζ
    for i_ϵ=1:n_ϵ
        # Sort A_endo for interpolation
        sort_ind = sortperm(A_endo[:,i_ϵ,i_ζ])
        A_aux    = A_endo[:,i_ϵ,i_ζ][sort_ind]
        C_aux    = C_endo[:,i_ϵ,i_ζ][sort_ind]
        # Check boundary condition
        if minimum(A_aux)>a_min
            a_vec = M.a_grid[M.a_grid.<minimum(A_aux)]
            A_aux = [a_vec ; A_aux]
            C_aux = [((r*M.ζ_mat[i_ζ].+1).*a_vec.+w*M.ϵ_grid[i_ϵ].-a_min) ; C_aux]
        end
        C_ip            = Spline1D(A_aux,C_aux)
        G_c[:,i_ϵ,i_ζ] .= C_ip.(M.a_grid)
        Ap_aux          = (r*M.ζ_mat[i_ζ].+1).*collect(M.a_grid) .+ w*M.ϵ_grid[i_ϵ] .- G_c[:,i_ϵ,i_ζ]
    end
    end
    # Update policy function
    G_ap .= (r*M.ζ_mat.+1).*M.a_mat .+ w*M.ϵ_mat .- G_c
        # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (r*M.ζ_mat[ind].+1)*M.a_mat[ind] + w*M.ϵ_mat[ind] - a_min
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
function Histogram_Method_Loop(M::Model,N_H=nothing,Γ_0=nothing)
    @unpack p, n_ϵ, n_ζ, MP_ϵ, MP_ζ, n_a_fine, a_grid_fine, θ_a, θ_a_f, G_ap_fine = M
    @unpack a_min, Hist_max_iter, Hist_tol, Hist_η = p

    println("\n--------------------------------\nBegining Histogram Method with Loops")

    # Change max iter
    if N_H==nothing
        N_H = Hist_max_iter
    end

    # Initial distribution
    if Γ_0==nothing
        Γ_0 = M.Γ
    end

    # Discretize distribution
    println("Discretizing Choices and Computing Transition Probabilities")
    H_ind    = Array{Int64}(undef,n_a_fine,n_ϵ,n_ζ)
    H_ω_lo   = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)
    H_ω_hi   = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)
    a_max    = maximum(a_grid_fine)
    for i_ζ=1:n_ζ # Current ζ
        Pr_ζp = MP_ζ.Π[i_ζ,:] ; # Transitions of future ζ conditional on current ζ
    for i_ϵ=1:n_ϵ # Current ϵ
        Pr_ϵp = MP_ϵ.Π[i_ϵ,:] ; # Transitions of future ϵ conditional on current ϵ
    for i_a=1:n_a_fine # Current a

        # Get index and weight of lower bound on approximation interval
        if     ( G_ap_fine[i_a,i_ϵ,i_ζ] ≥ a_max)
            H_ind[i_a,i_ϵ,i_ζ] = n_a_fine - 1
            ω_lo               = 0  
        elseif ( G_ap_fine[i_a,i_ϵ,i_ζ] ≤ a_min)
            H_ind[i_a,i_ϵ,i_ζ] = 1
            ω_lo               = 1 
        else
            H_ind[i_a,i_ϵ,i_ζ] = Grid_Inv(G_ap_fine[i_a,i_ϵ,i_ζ],n_a_fine,θ_a_f,a_min,a_max)
            ω_lo               = min(1,max(0,(G_ap_fine[i_a,i_ϵ,i_ζ]-a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]])/(a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]+1]-a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]])))
        end

        # Store weights for lower and upper bounds on approximation interval, including transition to future states 
        for i_ζp=1:n_ζ # Future ζ
        for i_ϵp=1:n_ϵ # Future ϵ
            H_ω_lo[i_a,i_ϵ,i_ζ,i_ϵp,i_ζp] = (  ω_lo)*Pr_ϵp[i_ϵp]*Pr_ζp[i_ζp]
            H_ω_hi[i_a,i_ϵ,i_ζ,i_ϵp,i_ζp] = (1-ω_lo)*Pr_ϵp[i_ϵp]*Pr_ζp[i_ζp]
        end 
        end 
    end
    end
    end
        # # Correct corner solutions above
        # H_weight[H_ind.==n_a_fine] .= 0
        # H_ind[H_ind.==n_a_fine]    .= n_a_fine-1
        # # Check bounds for weights
        # H_weight = min.(1,max.(0,H_weight))


    # Loop for updating histogram
    println("Iterating on the Distribution")
    H_dist = 1
    for i_H=1:N_H 
        # Update histogram
        Γ = zeros(n_a_fine,n_ϵ,n_ζ)
        for i_ζ=1:n_ζ # Current ζ
        for i_ϵ=1:n_ϵ # Current ϵ
        for i_a=1:n_a_fine # Current a
            i_ap = H_ind[i_a,i_ϵ,i_ζ]    ;
            for i_ζp=1:n_ζ # Future ζ
            for i_ϵp=1:n_ϵ # Future ϵ
                # Update is the product of probabilities by independence of F(ϵ) and F(ζ)
                Γ[i_ap,i_ϵp,i_ζp]   = Γ[i_ap  ,i_ϵp,i_ζp] + H_ω_lo[i_a,i_ϵ,i_ζ,i_ϵp,i_ζp]*Γ_0[i_a,i_ϵ,i_ζ]
                Γ[i_ap+1,i_ϵp,i_ζp] = Γ[i_ap+1,i_ϵp,i_ζp] + H_ω_hi[i_a,i_ϵ,i_ζ,i_ϵp,i_ζp]*Γ_0[i_a,i_ϵ,i_ζ]
            end
            end
        end
        end
        end
        # Update distance
        H_dist = maximum(abs.(Γ-Γ_0))
        # Update initial distribution
        Γ_0 .= (1-Hist_η)*Γ .+ Hist_η*Γ_0
        # Report progress
        if mod(i_H,10)==0
            @printf("\n   Histogram Loop: iter = %d, dist = %.4e",i_H,H_dist)
        end
        # Check convergence
        if H_dist<Hist_tol
            @printf("Histogram iteration converged in iteration %d. H_dist=%.4e \n--------------------------------\n",i_H,H_dist)
            M = Model(M; Γ=Γ, H_ind=H_ind, H_ω_lo=H_ω_lo, H_ω_hi=H_ω_hi)
            return M
        end
    end
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram iteration method
function Histogram_Iteration(M::Model,N_H,Γ_0)
    @unpack n_ϵ, n_ζ, n_a_fine, H_ind,H_ω_lo,H_ω_hi = M

    # println("\n--------------------------------\nBegining Histogram Iteration (N=$N_H)")

    for i_H=1:N_H 
        # Update histogram - Loop only though points with enough mass
        Γ = zeros(n_a_fine,n_ϵ,n_ζ)
        for ind in findall(>=(1e-12), Γ_0 )
            i_ap = H_ind[ind]    ;
            for i_ζp=1:n_ζ # Future ζ
            for i_ϵp=1:n_ϵ # Future ϵ
                Γ[i_ap,i_ϵp,i_ζp]   = Γ[i_ap  ,i_ϵp,i_ζp] + H_ω_lo[ind,i_ϵp,i_ζp]*Γ_0[ind]
                Γ[i_ap+1,i_ϵp,i_ζp] = Γ[i_ap+1,i_ϵp,i_ζp] + H_ω_hi[ind,i_ϵp,i_ζp]*Γ_0[ind]
            end
            end
        end 
        Γ_0 .= Γ 
        # println("   Iteration $i_H")
    end
    # println("End of Histogram Iteration \n--------------------------------\n")
    return Γ_0 
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Euler Error Function

# G_ap interpolation
function G_ap_ϵa(i_ζ::Int64,i_ϵ::Int64,a,M::Model)
    itp = ScaledInterpolations(M.a_grid,M.G_ap[:,i_ϵ,i_ζ], BSpline(Cubic(Line(OnGrid()))))
    return itp(a)
end

# Euler Equation Percentage Error
function Euler_Error(i_ϵ::Int64,i_ζ::Int64,a,M::Model)
    # Return percentage error in Euler equation
    @unpack p, MP_ϵ, MP_ζ, n_ϵ, n_ζ, ϵ_grid, ζ_grid = M
    @unpack β, r, w = p

    # Iterpolate G_ap at current ϵ
    ap  = min(M.a_grid[end],max(M.a_grid[1],G_ap_ϵa(i_ζ,i_ϵ,a,M)))

    # Current consumption
    c   = (1+r*ζ_grid[i_ζ])*a + w*ϵ_grid[i_ϵ] - ap

    # Compute left hand side of Euler equation
    LHS = d_utility(c,p)

    # Compute right hand side of Euler equation
    RHS_prob = zeros(n_ϵ,n_ζ)
        # Marginal utility at ϵ',ζ',a' ,G_ap(a',ϵ',ζ')
        for i_ϵp=1:n_ϵ
        for i_ζp=1:n_ζ
        cp  = (1+r*ζ_grid[i_ζp])*ap + w*ϵ_grid[i_ϵp] - G_ap_ϵa(i_ζp,i_ϵp,ap,M)
        up  = d_utility(cp,p)
        RHS_prob[i_ϵp,i_ζp] = β*(MP_ϵ.Π[i_ϵ,i_ϵp]*MP_ζ.Π[i_ζ,i_ζp])*((1+r*ζ_grid[i_ζp])*up)
        end
        end
    RHS = sum(RHS_prob)
    # Return percentage errror in Euler equation
    return (RHS/LHS-1)*100
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Aiyagari Equilibrium
function Aiyagari_Equilibrium(M_in::Model)

    # Read files if required 
    if M_in.read_flag
        G_ap_0 = reshape( readdlm(File_Folder*"/Policy_Function.csv", ',', Float64) ,M_in.n_a     , M.n_ϵ, M.n_ζ );
        Γ_0    = reshape( readdlm(File_Folder*"/Distribution.csv"   , ',', Float64) ,M_in.n_a_fine, M.n_ϵ, M.n_ζ );
    else 
        G_ap_0 = nothing ;
        Γ_0    = nothing ;
    end 

    # Compute Policy Functions
    M_in   = PFI_Fixed_Point(T_EGM_G,M_in,G_ap_0) ;

    # Compute Distribution
    M_in  = Histogram_Method_Loop(M_in,nothing,Γ_0) ;

    # Save Results 
    open(File_Folder*"/Policy_Function.csv", "w") do io
        writedlm(io, M_in.G_ap , ',')
    end;
    open(File_Folder*"/Distribution.csv", "w") do io
        writedlm(io, M_in.Γ , ',')
    end;

    return M_in
end