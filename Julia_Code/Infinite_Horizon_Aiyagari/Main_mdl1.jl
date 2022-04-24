# Computing Longitudinal Moments for Heterogeneous Agent Models
# Sergio Ocampo, Baxter Robinson, and Emmanuel Murray Leclair
# April 2022
# Aiyagari economy with inelastic labor supply
# 3 partial equilibrium models:
#       1. Aiyagari economy with infinitely lived agents and constant rate of returns
#       2. Aiyagari economy with infinitely lived agents and stochastic rate of returns
#       3. Overlaping generation Aiyagari economy with constant rate of returns
# This scripts computes longitudinal moments for model 1 (Aiyagari economy with infinitely lived agents and constant rate of returns)

## Make auxiliary directores 
mkpath("Figures")

# Load packages
using SparseArrays
using Plots
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl

# Load functions in VFI Toolbox 
include("../VFI_Toolbox.jl")

println(" ")
println("------------------------")
println("Aiyagari in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, Roots, Parameters, ScaledInterpolation")
println("Solve Aiyagari model with EGM and the histogram method")
println("------------------------")
println(" ")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Parameters and Model Structure
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par 
        # Model Parameters
        α::Float64 = 0.36 ; # Production function 
        β::Float64 = 0.96 ; # Discount factor
        γ::Float64 = 2.0  ; # Relative risk aversion (utility) parameter
        δ::Float64 = 0.05 ; # Depreciation rate
        ρ_ϵ::Float64 = 0.963 ; # Persistence of labor efficiency process
        σ_ϵ::Float64 = sqrt(0.162) ; # Standard deviation of labor efficiency innovation
        #ϵ̄::Float64 = exp(-σ_ϵ^2/(2*(1-ρ_ϵ^2))); # Reference level for labor efficiency
        # Model prices (partial equilibrium) and aggregates
        r::Float64 = 0.0379 ; # Average real return on net-worth (Fagereng et al. 2020)
        w::Float64 = 53.624 ; # U.S. (2019) - thousand $
        # Borrowing constraint
        a_min::Float64 = 0; # Borrowing constraint
        # VFI Parameters
        max_iter::Int64   = 100000; # Maximum number of iterations
        dist_tol::Float64 = 1E-6  ; # Tolerance for distance
        # Howard's Policy Iterations
        H_tol::Float64    = 1E-9  ; # Tolerance for policy function iteration
        N_H::Int64        = 20    ; # Maximum number of policy iterations
        # Histogram iteration parameters
        Hist_max_iter     = 10000 ;
        Hist_tol          = 1E-7  ;
        # Histogram iteration parameters
        N_eq              = 1000  ;
        tol_eq            = 1E-7  ;
        η                 = 0.3   ; # Dampen factor for updating capital
        # Minimum consumption for numerical optimization
        c_min::Float64    = 1E-16
    end

# Allocate paramters to object p for future calling
p = Par();

# Generate structure of model objects
    # Model 1 : Aiyagari economy with infinitely lived agents and constant rate of returns
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model parameters in their own structure
        # Capital Grid
        a_max::Float64  = 50                         # Max node of a_grid
        θ_a::Float64    = 2.5                        # Curvature of a_grid
        n_a::Int64      = 200                        # Size of a_grid
        n_a_fine::Int64 = 1000                       # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,1  ,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        # Productivity process
        n_ϵ       = 15                               # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ_ϵ,p.σ_ϵ,n_ϵ)   # Markov Process for ϵ
        ϵ_ref     = n_ϵ/sum(exp.(MP_ϵ.grid))         # Reference level for labor efficiency
        ϵ_grid    = ϵ_ref*exp.(MP_ϵ.grid)            # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid,1,n_ϵ)
        a_mat_fine= repeat(a_grid_fine,1,n_ϵ)
        ϵ_mat     = repeat(ϵ_grid',n_a,1)
        # Value and policy functions
        V         = Array{Float64}(undef,n_a,n_ϵ)       # Value Function
        G_ap      = Array{Float64}(undef,n_a,n_ϵ)       # Policy Function for capital
        G_c       = Array{Float64}(undef,n_a,n_ϵ)       # Policy Function
        V_fine    = Array{Float64}(undef,n_a_fine,n_ϵ)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_a_fine,n_ϵ)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_a_fine,n_ϵ)  # Policy Function on fine grid
        # Distribution
        Γ         = 1/(n_ϵ*n_a_fine)*ones(n_a_fine,n_ϵ) # Distribution (initiliazed to uniform)
        # Solver
        Solver    = "PFI"
    end

M = Model();

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
# Policy function: PFI Fixed Point 
function PFI_Fixed_Point(T::Function,M::Model,G_ap_old=nothing)
    # Unpack model structure
    @unpack p, n_ϵ, n_a, n_a_fine, θ_a, a_grid, a_grid_fine = M
    # PFI paramters
    @unpack max_iter, dist_tol, r, w = p
    # Initialize variables for loop
    if G_ap_old==nothing
        G_ap_old  = (1+r)*M.a_mat
    end
    G_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("PFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
    for iter=1:max_iter
        # Update value function
        G_ap_new, G_c = T(Model(M,G_ap=copy(G_ap_old)))
        # Update distance and iterations
        G_dist = sqrt(norm(G_ap_new-G_ap_old,2))
        # Update old function
        G_ap_old  = G_ap_new
        # Report progress
        if mod(iter,250)==0
            println("   PFI Loop: iter=$iter, dist=",G_dist)
        end
        # Check convergence and return results
        if G_dist<=dist_tol
            println("PFI - n_ϵ=$n_ϵ, n_a=$n_a - θ_a=$θ_a - r=$r")
            println("Iterations = $iter and Distance = ",G_dist)
            println("------------------------")
            println(" ")
            # Interpolate to fine grid
            G_ap_fine = zeros(n_a_fine,n_ϵ)
            G_c_fine  = zeros(n_a_fine,n_ϵ)
            for i_ϵ=1:n_ϵ
            G_ap_ip = ScaledInterpolations(a_grid,G_ap_new[:,i_ϵ] , BSpline(Cubic(Line(OnGrid()))))
                G_ap_fine[:,i_ϵ].= G_ap_ip.(collect(a_grid_fine))
            G_c_ip  = ScaledInterpolations(a_grid,G_c[:,i_ϵ]  , BSpline(Cubic(Line(OnGrid()))))
                G_c_fine[:,i_ϵ] .= G_c_ip.(collect(a_grid_fine))
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
    @unpack p, n_ϵ, MP_ϵ, n_a, G_ap = M
    @unpack β, a_min, r, w = p
    # Define RHS of Euler equation for each (ϵ,a')
    # Rows are present ϵ and columns are tomorrow's a in fixed grid
    Euler_RHS = (β*(1+r)*MP_ϵ.Π*d_utility( (1+r)*M.a_mat + w*M.ϵ_mat - G_ap , p )')'
    # Check Monotonicity
    if any( Euler_RHS.<0 )
        error("RHS must be monotone for EGM to work")
    end
    # Define consumption from Euler equation
    C_endo = d_utility_inv(Euler_RHS,p)
    # Define endogenous grid on assets
    A_endo = (C_endo .+ M.a_mat - w*M.ϵ_mat)/(1+r)
    # Interpolate functions on exogenous grid
    G_c = Array{Float64}(undef,n_a,n_ϵ)
    for i_ϵ=1:n_ϵ
        # Sort A_endo for interpolation
        sort_ind = sortperm(A_endo[:,i_ϵ])
        A_aux    = A_endo[:,i_ϵ][sort_ind]
        C_aux    = C_endo[:,i_ϵ][sort_ind]
        # Check boundary condition
        if minimum(A_aux)>a_min
            a_vec = M.a_grid[M.a_grid.<minimum(A_aux)]
            A_aux = [a_vec ; A_aux]
            C_aux = [((1+r)*a_vec.+w*M.ϵ_grid[i_ϵ].-a_min) ; C_aux]
        end
        C_ip        = Spline1D(A_aux,C_aux)
        G_c[:,i_ϵ] .= C_ip.(M.a_grid)
        Ap_aux      = (1+r)*collect(M.a_grid) .+ w*M.ϵ_grid[i_ϵ] .- G_c[:,i_ϵ]
    end
    # Update policy function
    G_ap .= (1+r)*M.a_mat .+ w*M.ϵ_mat .- G_c
        # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (1+r)*M.a_mat[ind] + w*M.ϵ_mat[ind] - a_min
        end
        # Check for borrowing constraint
        if any( G_ap.<a_min )
            error("Borrowing Constraint Violated")
        end
    # Return Results
    return G_ap, G_c
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Indicator variable function (returns 1 if x is between upper and lower bound)
function Indicator(lb,ub,x)
    if x > lb && x < ub
        return 1
    else
        return 0
    end
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method
function Histogram_Method_Loop(M::Model,N_H=nothing,Γ_0=nothing)
    @unpack p, n_ϵ, MP_ϵ, n_a_fine, a_grid_fine, G_ap_fine = M
    @unpack a_min, Hist_max_iter, Hist_tol = p

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
    H_ind    = Array{Int64}(undef,n_a_fine,n_ϵ)
    H_weight = Array{Float64}(undef,n_a_fine,n_ϵ)
    a_max    = maximum(a_grid_fine)
    for i_ϵ=1:n_ϵ
    for i_a=1:n_a_fine
        H_ind[i_a,i_ϵ]    = Grid_Inv(G_ap_fine[i_a,i_ϵ],n_a_fine,1,a_min,a_max)
        H_weight[i_a,i_ϵ] = (G_ap_fine[i_a,i_ϵ]-a_grid_fine[H_ind[i_a,i_ϵ]])/(a_max-a_min)
    end
    end
        # Correct corner solutions above
        H_weight[H_ind.==n_a_fine] .= 0
        H_ind[H_ind.==n_a_fine]    .= n_a_fine-1
        # Check bounds for weights
        H_weight = min.(1,max.(0,H_weight))

    # Loop for updating histogram
    H_dist = 1
    for i_H=1:N_H
        #println("iteration = $i_H")
        # Update histogram
        Γ = zeros(n_a_fine,n_ϵ)
        for i_ϵ=1:n_ϵ # Current ϵ
        for i_a=1:n_a_fine # Current a
            i_ap = H_ind[i_a,i_ϵ]
            ω_ap = H_weight[i_a,i_ϵ]
            for i_ϵp=1:n_ϵ # Future ϵ
                Γ[i_ap,i_ϵp]   = Γ[i_ap,i_ϵp]   +    ω_ap *MP_ϵ.Π[i_ϵ,i_ϵp]*Γ_0[i_a,i_ϵ]
                Γ[i_ap+1,i_ϵp] = Γ[i_ap+1,i_ϵp] + (1-ω_ap)*MP_ϵ.Π[i_ϵ,i_ϵp]*Γ_0[i_a,i_ϵ]
            end
        end
        end
        # Update distance
        H_dist = maximum(abs.(Γ-Γ_0))
        # Update initial distribution
        Γ_0 .= Γ
        # Report progress
        if mod(i_H,250)==0
            println("   Histogram Loop: iter=$i_H, dist=$H_dist")
        end
        # Check convergence
        if H_dist<Hist_tol
            println("Histogram iteartion converged in iteration $i_H. H_dist=$H_dist\n--------------------------------\n")
            M = Model(M; Γ=Γ)
            return M
        end
    end
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Cross-Sectional Moments
function CrossSection_Moments(M::Model,X,s_upper,s_lower)
    @unpack Γ = M
    # Vectorize the stationary distribution and variable of interest X
    Γ_v = vec(Γ)
    X_v = vec(X)
    # Check that dimensions match
    if size(X_v,1) != size(Γ_v)
        error("Dimensions do not match")
    end
    # Impose condition
    for i = 1:size(Γ_v)
        if Γ_v[i] > s_upper
            Γ_v[i] = 0
            X_v[i] = 0
        elseif Γ_v[i] < s_lower
            Γ_v[i] = 0
            X_v[i] = 0
        end
    end
    # Compute marginal distribution
    Γ_v_marg = Γ_v_marg./(sum(Γ_v))
    # Compute moment
    return sum(X_v.*Γ_v_marg)
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Aiyagari Equilibrium
function Aiyagari_Equilibrium(M::Model)
    @unpack p, Solver = M
    @unpack β, α, δ = p

    # Compute Policy Functions
    if Solver=="PFI"
    M   = PFI_Fixed_Point(T_EGM_G,M)
    end

    # Compute Distribution
    M  = Histogram_Method_Loop(M)

    # Return Model
    return M

    # No convergence, Display error
    error("No convervence to equilibrium after $N_eq iterations. Current distance of capital: $(100*K_dist)%")
end

println("===============================================\n Solving Aiyagari with EGM-Histogram(loop)")
@time M_Aiyagari = Aiyagari_Equilibrium(Model(Solver="PFI"))

