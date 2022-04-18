# Computing Longitudinal Moments for Heterogeneous Agent Models
# Sergio Ocampo, Baxter Robinson, and Emmanuel Murray Leclair
# April 2022
# Aiyagari economy with inelastic labor supply
# 3 partial equilibrium models:
#       1. Aiyagari economy with infinitely lived agents and constant rate of returns
#       2. Aiyagari economy with infinitely lived agents and stochastic rate of returns
#       3. Overlaping generation Aiyagari economy with constant rate of returns
# This scripts computes longitudinal moments for model 2 (Aiyagari economy with infinitely lived agents and stochastic rate of returns)

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
        ρ_ζ::Float64 = 0.80 ; # Persistence of interest rate
        σ_ζ::Float64 = sqrt(0.086) ; # Standard deviation of interest rate innovation
        z_bar::Float64 = 1; # Reference level for productivity
        ϵ̄::Float64 = exp(-σ_ϵ^2/(2*(1-ρ_ϵ^2))); # Reference level for labor efficiency
        ζ::Float64 = exp(-σ_ζ^2/(2*(1-ρ_ζ^2))); # Reference level for interest rate
        # Model prices (partial equilibrium) and aggregates
        r::Float64 = 0.0379 ; # Average real return on net-worth (Fagereng et al. 2020)
        LK::Float64 = (r/(α*z_bar))^(1/(1-α)) ; # Aggregate capital to labor ratio
        w::Float64 = (1-α)*((LK)^(-α)) ; # Implied wage
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
p = Par()

# Generate structure of model objects
    # Model 1 : Aiyagari economy with infinitely lived agents and constant rate of returns
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model parameters in their own structure
        # Capital Grid
        a_max::Float64  = 50                         # Max node of a_grid
        θ_a::Float64    = 2.5                        # Curvature of a_grid
        n_a::Int64      = 200                        # Size of a_grid
        n_a_fine::Int64 = 1000                        # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,1  ,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        # Interest rate process
        n_ζ       = 5                                # Size of ζ_grid
        MP_ζ      = Rouwenhorst95(p.ρ_ζ,p.σ_ζ,n_ζ)   # Markov Process for ζ
        ζ_grid    = p.ζ*exp.(MP_ζ.grid)              # Grid in levels
        # Productivity process
        n_ϵ       = 15                               # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ_ϵ,p.σ_ϵ,n_ϵ)   # Markov Process for ϵ
        ϵ_grid    = p.ϵ̄*exp.(MP_ϵ.grid)              # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid',n_ϵ,1,n_ζ)
        a_mat_fine= repeat(a_grid_fine',n_ϵ,1,n_ζ)
        ϵ_mat     = p.ϵ̄*exp.(repeat(MP_ϵ.grid,1,n_a,n_ζ))
        ζ_mat     = repeat(reshape(MP_ζ.grid,1,1,3),n_ϵ,n_a,1)
        # Value and policy functions
        V         = Array{Float64}(undef,n_ϵ,n_a,n_ζ)       # Value Function
        G_ap      = Array{Float64}(undef,n_ϵ,n_a,n_ζ)       # Policy Function for capital
        G_c       = Array{Float64}(undef,n_ϵ,n_a,n_ζ)       # Policy Function
        V_fine    = Array{Float64}(undef,n_ϵ,n_a_fine,n_ζ)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_ϵ,n_a_fine,n_ζ)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_ϵ,n_a_fine,n_ζ)  # Policy Function on fine grid
        # Distribution
        Γ         = 1/(n_ϵ*n_a_fine*n_ζ)*ones(n_ϵ,n_a_fine,n_ζ) # Distribution (initiliazed to uniform)
        # Solver
        Solver    = "PFI"
    end

M = Model()

# Outside the model (works)
ζ_mat     = zeros(M.n_ϵ,M.n_a,M.n_ζ)
for i_ζ=1:M.n_ζ
    ζ_mat[:,:,i_ζ] = p.ζ*exp.(fill(M.MP_ζ.grid[i_ζ],M.n_ϵ,M.n_a))
end
