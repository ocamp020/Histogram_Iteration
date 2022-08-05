# Computing Longitudinal Moments for Heterogeneous Agent Models
# Sergio Ocampo, Baxter Robinson, and Emmanuel Murray Leclair
# July 2022
# Aiyagari economy: 
#       1. Infinitely lived agents
#       2. Inelastic labor supply
#       3. Stochastic rate of returns
# This scripts computes longitudinal moments for the model
# 
# Solve:   V(ζ,ϵ,a) = max{ ((1+r(ζ))a+wϵ̄ϵ-a')^(1-σ)/(1-σ) +beta*E[V(ζ',ϵ',a')|ζ,ϵ] }
#           log(ϵ') = ρ_ϵ*log(ϵ) + η_ϵ; η_ϵ~N(0,σ_ϵ);
#           r(ζ)    = exp(ζ)r⋆    
#           log(ζ') = ρ_ζ*log(ζ) + η_ζ; η_ζ~N(0,σ_ζ); 
# The constant ϵ̄ guarantees that E[ϵ]=1 and so aggregate labor L=E[ϵ]=1

## Change to your home directory 
# # Sergio's Computer 
#    cd()
#    cd("./Dropbox/Research/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/")
# Emmanuel's Computer
    # cd()
    # cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/") # Laptop
    # cd("D:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/") # Desktop
    # cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/")
# Baxter's Computer
     cd("D:/Dropbox/Files/Economics-Research/Project-09_SIM/Code/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/")
# Compute Canada Server
#    cd("/scratch/robin370/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/")


## Make auxiliary directores
    Fig_Folder  = "Figures" ; mkpath(Fig_Folder)  ;
    File_Folder = "Files"   ; mkpath(File_Folder) ;
    Hist_Folder = "Files/Histogram"  ; mkpath(Hist_Folder) ;
    MC_Folder   = "Files/MonteCarlo" ; mkpath(MC_Folder)   ;

# Load packages
using SparseArrays
using Plots
using Colors
using LaTeXStrings
using Interpolations # Pkg.add("Interpolations") # https://github.com/JuliaMath/Interpolations.jl
using Dierckx # Pkg.add("Dierckx") # https://github.com/kbarbary/Dierckx.jl
using ForwardDiff # Pkg.add("ForwardDiff") # https://github.com/JuliaDiff/ForwardDiff.jl
using Optim # Pkg.add("Optim") # https://julianlsolvers.github.io/Optim.jl/stable/
using Optim: converged, maximum, maximizer, minimizer, iterations
using Roots # Pkg.add("Roots") # https://github.com/JuliaMath/Roots.jl
using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
using Kronecker # Pkg.add("Kronecker") # https://michielstock.github.io/Kronecker.jl
using DelimitedFiles
using Printf
using StatsBase
using Distributions
using TimerOutputs # Pkg.add("TimerOutputs")

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
        β::Float64   = 0.94  ; # Discount factor
        σ::Float64   = 2.0   ; # Relative risk aversion (utility) parameter
        ρ_ϵ::Float64 = 0.963 ; # Persistence of labor efficiency process
        σ_ϵ::Float64 = 0.162 ; # Standard deviation of labor efficiency innovation
        ρ_ζ::Float64 = 0.70  ; # Persistence of interest rate target pareto coefficient of 1.8
        σ_ζ::Float64 = 1.30  ; # Standard deviation of interest rate target top 1% share of 20%
        # Model prices (partial equilibrium)
        r::Float64 = 0.0320 ; #  Target wealth weighted 3.79% average real return on net-worth (Fagereng et al. 2020)
        w::Float64 = 53.624 ; # U.S. (2019) - thousands of dollars $
        # Borrowing constraint
        a_min::Float64 = 0.010 ; # Borrowing constraint (10 dollars)
        # VFI Parameters
        max_iter::Int64     = 20000 ; # Maximum number of iterations
        dist_tol::Float64   = 5E-6  ; # Tolerance for distance
        dist_tol_Δ::Float64 = 1E-10 ; # Tolerance for change in distance 
        η                   = 0.10  ; # Dampen factor
        # Histogram iteration parameters
        Hist_max_iter       = 1500  ; # Maximum number of iterations
        Hist_tol            = 1E-8  ; # Tolerance for distance
        Hist_η              = 0.1   ; # Dampen factor
        # Minimum consumption for numerical optimization
        c_min::Float64      = 1E-16
    end

# Allocate paramters to object p for future calling
p = Par();

# Generate structure of model objects
    # Model 1 : Aiyagari economy with infinitely lived agents and constant rate of returns
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model parameters in their own structure
        # Assets Grid
        a_max::Float64  = 100000                     # Max node of a_grid
        θ_a::Float64    = 4.5                        # Curvature of a_grid
        θ_a_f::Float64  = 4.5                        # Curvature of a_grid_fine
        n_a::Int64      = 250                        # Size of a_grid
        n_a_fine::Int64 = 500                        # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a  ,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,θ_a_f,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        n_cut_fine      = Grid_Inv(25000,n_a_fine,θ_a_f,p.a_min,a_max) # Index just below 1000
        # Interest rate process
        n_ζ       = 7                                  # Size of ζ_grid
        MP_ζ      = Tauchen86(p.ρ_ζ,p.σ_ζ,n_ζ,1.96)      # Markov Process for ζ
        ζ_ref     = 1/sum(exp.(MP_ζ.grid).*MP_ζ.PDF)   # Reference level for interest rate
        ζ_grid    = ζ_ref*exp.(MP_ζ.grid)              # Grid in levels
        # Labor productivity process
        n_ϵ       = 15                                 # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ_ϵ,p.σ_ϵ,n_ϵ)     # Markov Process for ϵ
        ϵ_ref     = 1/sum(exp.(MP_ϵ.grid).*MP_ϵ.PDF)   # Reference level for labor efficiency 
        ϵ_grid    = ϵ_ref*exp.(MP_ϵ.grid)              # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid,1,n_ϵ,n_ζ)
        ϵ_mat     = repeat(ϵ_grid',n_a,1,n_ζ)
        ζ_mat     = repeat(reshape(ζ_grid,(1,1,n_ζ)),n_a,n_ϵ,1)
        a_mat_fine= repeat(a_grid_fine,1,n_ϵ,n_ζ)
        ϵ_mat_fine= repeat(ϵ_grid',n_a_fine,1,n_ζ)
        ζ_mat_fine= repeat(reshape(ζ_grid,(1,1,n_ζ)),n_a_fine,n_ϵ,1)
        # Value and policy functions
        V         = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Value Function
        G_ap      = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Policy Function for capital
        G_c       = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Policy Function
        V_fine    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Policy Function on fine grid
        # Distribution
        Γ         = 1/(n_cut_fine*n_ϵ*n_ζ)*[ones(n_cut_fine,n_ϵ,n_ζ) ; zeros(n_a_fine-n_cut_fine,n_ϵ,n_ζ)]     # Distribution (initiliazed to uniform)
        H_ind     = Array{Int64}(undef,n_a_fine,n_ϵ,n_ζ)            # Index for discretization of savings choice 
        H_ω_lo    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)  # Transition probabilities to future states (lower bound)
        H_ω_hi    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)  # Transition probabilities to future states (lower bound)
        # Misc
        method = 1 # 1 for Kronecker and 2 for loops in expectation of PFI
        read_flag = false # Boolean for reading results from file 
    end

M = Model();

# Load functions in Functions_ModelSolution (solve the model and find stationary distribution)
include("Functions_ModelSolution.jl")

# Load functions in Functions_Montecarlo (Simulate panels of individual agents)
include("Functions_MonteCarlo.jl")


# Execute model solution 
println("\n===============================================\n Solving Aiyagari with EGM-Histogram(loop)")
    
    @time M_Aiyagari = Aiyagari_Equilibrium(Model(method=1,read_flag=false));

println("===============================================\n")

 
# # Get stats and graphs for the solution of the model 
# include("PrintStats_MakeGraphs.jl")


# # Get moments from histogram method
# include("CalculateMoments_Histogram.jl")


# # Get moments from simulation
# include("CalculateMoments_MonteCarlo.jl")


# # Run Draft Moments for Graphs and Tables 
include("Draft_Results.jl")

# Make Draft Graphs and Tables
include("Draft_Graphs_Tables.jl")


println("\n===============================================\n\n    End of Script \n\n===============================================")

