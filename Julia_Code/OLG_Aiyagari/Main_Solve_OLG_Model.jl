# Computing Longitudinal Moments for Heterogeneous Agent Models
# Sergio Ocampo, Baxter Robinson, and Emmanuel Murray Leclair
# April 2022
# Aiyagari economy: 
#       1. Finitely lived agents and overlapping generations
#       2. Inelastic labor supply
# This scripts computes longitudinal moments for the model
# 
# Solve:   V(h,ϵ,a) = max{ ((1+r)a+wϵ̄f(h,ϵ)-a')^(1-γ)/(1-γ) + beta*E[ s(h)*V(h+1,ϵ',a') + (1-s(h))*ν(a')|ϵ] }
#           log(ϵ') = ρ_ϵ*log(ϵ) + η_ϵ; η_ϵ~N(0,σ_ϵ); 
# The constant ϵ̄ guarantees that E[f(h,ϵ)]=1 and so aggregate labor L=E[ϵ]=1

## Change to your home directory 
# Sergio's Computer 
#    cd()
#    cd("./Dropbox/Research/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
# Emmanuel's Computer
    # cd()
    # cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/") # Laptop
    # cd("D:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/") # Desktop
    # cd("C:/Users/Emmanuel/Dropbox/RA_Sergio/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
# Baxter's Computer
     cd("D:/Dropbox/Files/Economics-Research/Project-09_SIM/Code/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")
# Compute Canada Server
    # cd("/scratch/robin370/Histogram_Iteration/Julia_Code/OLG_Aiyagari/")

## Make auxiliary directores
    Fig_Folder  = "Figures" ; mkpath(Fig_Folder)  ;
    File_Folder = "Files"   ; mkpath(File_Folder) ;
    Hist_Folder = "Files/Histogram"  ; mkpath(Hist_Folder) ;
    MC_Folder   = "Files/MonteCarlo" ; mkpath(MC_Folder)   ;


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
using Kronecker # Pkg.add("Kronecker") # https://michielstock.github.io/Kronecker.jl
using DelimitedFiles
using Printf
using StatsBase

# Load functions in VFI Toolbox
include("../VFI_Toolbox.jl")

println(" ")
println("------------------------")
println("Aiyagari in Julia")
println("PWD: ",pwd())
println("This code uses Plots, Interpolations, Dierckx, ForwardDiff, Optim, Roots, Parameters, ScaledInterpolation")
println("Solve OLG Aiyagari model with EGM and the histogram method")
println("------------------------")
println(" ")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Load function with population and survival probabilities
include("Setup_Demographics.jl")

#-----------------------------------------------------------
#-----------------------------------------------------------
# Parameters and Model Structure
    # Generate structure for parameters using Parameters module
    # We can set default values for our parameters
    @with_kw struct Par
        # Model Parameters
        β::Float64   = 0.94   ; # Discount factor
        σ::Float64   = 2.0    ; # Relative risk aversion (utility) parameter
        ρ_ϵ::Float64 = 0.963  ; # Persistence of labor efficiency process
        σ_ϵ::Float64 = 0.162  ; # Standard deviation of labor efficiency innovation
        # Model prices (partial equilibrium)
        r::Float64 = 0.0379   ; # Target wealth weighted 3.79% average real return on net-worth (Fagereng et al. 2020)
        w::Float64 = 53.624   ; # U.S. (2019) - tens of thousands $
        # Borrowing constraint
        a_min::Float64 = 1E-4 ; # Borrowing constraint
        # Histogram iteration parameters
        Hist_max_iter       = 1000  ; # Maximum number of iterations
        Hist_tol            = 1E-6  ; # Tolerance for distance
        Hist_η              = 0.00  ; # Dampen factor
        # Minimum consumption for numerical optimization
        c_min::Float64      = 1E-16 ; 
        # Life cycle parameters 
        Max_Age             = 81    ; # Corresponds to 100 years old 
        Surv_Pr             = Survival_Probabilities_Bell_Miller(Max_Age) ; 
        Age_Π               = Age_Transition(Max_Age,Surv_Pr)             ;
        Age_PDF             = Age_Distribution(Age_Π)                     ;
    end

# Allocate paramters to object p for future calling
p = Par();

# Generate structure of model objects
    # Model 1 : Aiyagari economy with infinitely lived agents and constant rate of returns
    @with_kw struct Model
        # Parameters
        p::Par = Par() # Model parameters in their own structure
        # Assets Grid
        a_max::Float64  = 10000                      # Max node of a_grid
        θ_a::Float64    = 3.5                        # Curvature of a_grid
        θ_a_f::Float64  = 3.5                        # Curvature of a_grid_fine
        n_a::Int64      = 250                        # Size of a_grid
        n_a_fine::Int64 = 500                        # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a  ,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,θ_a_f,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        # Labor productivity process - Transitory 
        n_ϵ       = 15                                 # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ_ϵ,p.σ_ϵ,n_ϵ)     # Markov Process for ϵ
        ϵ_ref     = 1.038479216975849/sum(exp.(MP_ϵ.grid).*MP_ϵ.PDF)   # Reference level for labor efficiency 
        ϵ_grid    = ϵ_ref*exp.(MP_ϵ.grid)              # Grid in levels
        # Labor productivity process - Life Cycle 
        age_vec   = collect(1:p.Max_Age)
        log_ξ_grid= (60*(age_vec.-1).-(age_vec.-1).^2)./1800  # Process peaks at age 50, and by age 80 gives the same income as when newborn
        ξ_ref     = 1/sum(exp.(log_ξ_grid).*p.Age_PDF)   # Reference level for labor efficiency 
        ξ_grid    = ξ_ref*exp.(log_ξ_grid)             # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid,1,n_ϵ,p.Max_Age)
        ϵ_mat     = repeat(ϵ_grid',n_a,1,p.Max_Age)
        ξ_mat     = repeat(reshape(ξ_grid,(1,1,p.Max_Age)),n_a,n_ϵ,1)
        a_mat_fine= repeat(a_grid_fine,1,n_ϵ,p.Max_Age)
        ϵ_mat_fine= repeat(ϵ_grid',n_a_fine,1,p.Max_Age)
        ξ_mat_fine= repeat(reshape(ξ_grid,(1,1,p.Max_Age)),n_a_fine,n_ϵ,1)
        a_mat_aϵ  = repeat(a_grid,1,n_ϵ)
        # Labor income matrices 
        y_mat     = p.w*ϵ_mat.*ξ_mat
        y_mat_fine= p.w*ϵ_mat_fine.*ξ_mat_fine
        # Value and policy functions
        V         = Array{Float64}(undef,n_a,n_ϵ,p.Max_Age)       # Value Function
        G_ap      = Array{Float64}(undef,n_a,n_ϵ,p.Max_Age)       # Policy Function for capital
        G_c       = Array{Float64}(undef,n_a,n_ϵ,p.Max_Age)       # Policy Function
        V_fine    = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age)  # Policy Function on fine grid
        # Distribution
        for i_h=1;p.Max_Age
        for i_ϵ=1:n_ϵ
        Γ[:,i_ϵ,i_h] .= 1/(n_a_fine)*MP_ϵ.PDF[i_ϵ]*p.Age_PDF[i_h] # Distribution (initiliazed to uniform)
        end
        end
        H_ind     = Array{Int64}(undef,n_a_fine,n_ϵ,p.Max_Age)              # Index for discretization of savings choice 
        # H_ω       = Array{Int64}(undef,n_a_fine,n_ϵ,p.Max_Age)                    # Probability of transition lo low index in discretization
        H_ω_lo_s   = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age,n_ϵ)
        H_ω_hi_s   = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age,n_ϵ)
        H_ω_lo_d   = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age)
        H_ω_hi_d   = Array{Float64}(undef,n_a_fine,n_ϵ,p.Max_Age)
        # Misc
        read_flag = false # Boolean for reading results from file 
    end

M = Model();

# Load functions in Functions_ModelSolution (solve the model and find stationary distribution)
include("Functions_ModelSolution.jl")


# Execute model solution 
println("\n===============================================\n Solving Aiyagari with EGM-Histogram(loop)")
    
    @time M_Aiyagari = Aiyagari_Equilibrium(Model(read_flag=true));

println("===============================================\n")


# # Get stats and graphs for the solution of the model 
# include("PrintStats_MakeGraphs.jl")


# # Get moments from histogram method
# include("CalculateMoments_Histogram.jl")


# # Get moments from simulation
# include("CalculateMoments_MonteCarlo.jl")

# Add Simulation Functions
include("Functions_MonteCarlo.jl")


# Run Draft Moments for Graphs and Tables 
include("Draft_Results.jl")


println("\n===============================================\n\n    End of Script \n\n===============================================")