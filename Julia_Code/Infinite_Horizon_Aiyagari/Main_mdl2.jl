# Computing Longitudinal Moments for Heterogeneous Agent Models
# Sergio Ocampo, Baxter Robinson, and Emmanuel Murray Leclair
# April 2022
# Aiyagari economy with inelastic labor supply
# 3 partial equilibrium models:
#       1. Aiyagari economy with infinitely lived agents and constant rate of returns
#       2. Aiyagari economy with infinitely lived agents and stochastic rate of returns
#       3. Overlaping generation Aiyagari economy with constant rate of returns
# This scripts computes longitudinal moments for model 2 (Aiyagari economy with infinitely lived agents and stochastic rate of returns)

## Change to your home directory 
# Sergio's Computer 
    cd()
    cd("./Dropbox/Research/Histogram_Iteration/Julia_Code/Infinite_Horizon_Aiyagari/")
# Emmanuel's Computer
    # cd("???")



## Make auxiliary directores
Fig_Folder = "Figures"
mkpath(Fig_Folder) 

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
        # Model prices (partial equilibrium) and aggregates
        r::Float64 = 0.0379 ; # Average real return on net-worth (Fagereng et al. 2020)
        w::Float64 = 53.624 ; # U.S. (2019) - tens of thousands $
        # Borrowing constraint
        a_min::Float64 = 1E-4 ; # Borrowing constraint
        # VFI Parameters
        max_iter::Int64   = 100000; # Maximum number of iterations
        dist_tol::Float64 = 1E-6  ; # Tolerance for distance
        η                 = 0.0   ; # Dampen factor
        # Howard's Policy Iterations
        # Histogram iteration parameters
        Hist_max_iter     = 10000 ; # Maximum number of iterations
        Hist_tol          = 1E-5  ; # Tolerance for distance
        Hist_η            = 0.0   ; # Dampen factor
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
        a_max::Float64  = 1000                       # Max node of a_grid
        θ_a::Float64    = 2.5                        # Curvature of a_grid
        n_a::Int64      = 200                        # Size of a_grid
        n_a_fine::Int64 = 1000                       # Size of fine grid for interpolation and distribution
        a_grid          = Make_Grid(n_a     ,θ_a,p.a_min,a_max,"Poly")  # a_grid for model solution
        a_grid_fine     = Make_Grid(n_a_fine,θ_a,p.a_min,a_max,"Poly")  # Fine grid for interpolation
        # Interest rate process
        n_ζ       = 3                                # Size of ζ_grid
        MP_ζ      = Rouwenhorst95(p.ρ_ζ,p.σ_ζ,n_ζ)   # Markov Process for ζ
        ζ_ref     = n_ζ/sum(exp.(MP_ζ.grid))         # Reference level for interest rate
        ζ_grid    = ζ_ref*exp.(MP_ζ.grid)            # Grid in levels
        # Productivity process
        n_ϵ       = 3                                # Size of ϵ_grid
        MP_ϵ      = Rouwenhorst95(p.ρ_ϵ,p.σ_ϵ,n_ϵ)   # Markov Process for ϵ
        ϵ_ref     = n_ϵ/sum(exp.(MP_ϵ.grid))         # Reference level for labor efficiency 
        ϵ_grid    = ϵ_ref*exp.(MP_ϵ.grid)            # Grid in levels
        # State matrices
        a_mat     = repeat(a_grid,1,n_ϵ,n_ζ)
        a_mat_fine= repeat(a_grid_fine,1,n_ϵ,n_ζ)
        ϵ_mat     = repeat(ϵ_grid',n_a,1,n_ζ)
        ζ_mat     = repeat(reshape(ζ_grid,(1,1,n_ζ)),n_a,n_ϵ,1)
        # Value and policy functions
        V         = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Value Function
        G_ap      = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Policy Function for capital
        G_c       = Array{Float64}(undef,n_a,n_ϵ,n_ζ)       # Policy Function
        V_fine    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Value Function on fine grid
        G_ap_fine = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Policy Function on fine grid
        G_c_fine  = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ)  # Policy Function on fine grid
        # Distribution
        Γ         = 1/(n_ϵ*n_a_fine*n_ζ)*ones(n_a_fine,n_ϵ,n_ζ)     # Distribution (initiliazed to uniform)
        H_ind     = Array{Int64}(undef,n_a_fine,n_ϵ,n_ζ)            # Index for discretization of savings choice 
        H_ω_lo    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)  # Transition probabilities to future states (lower bound)
        H_ω_hi    = Array{Float64}(undef,n_a_fine,n_ϵ,n_ζ,n_ϵ,n_ζ)  # Transition probabilities to future states (lower bound)
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
    @unpack p, n_ϵ, n_ζ, n_a, n_a_fine, θ_a, a_grid, a_grid_fine,method = M
    # PFI paramters
    @unpack max_iter, dist_tol, r, w, PFI_η = p
    # Initialize variables for loop
    if G_ap_old==nothing
        G_ap_old  = (r*M.ζ_mat.+1).*M.a_mat
    end
    G_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("PFI - n_ϵ=$n_ϵ, n_ζ=$n_ζ, n_a=$n_a - θ_a=$θ_a - r=$r")
    for iter=1:max_iter
        # Update value function
        G_ap_new, G_c = T(Model(M,G_ap=copy(G_ap_old)))
        # Update distance and iterations
        G_dist = sqrt(norm(G_ap_new-G_ap_old,2))
        # Update old function
        G_ap_old  = (1-PFI_η)*G_ap_new .+ PFI_η*G_ap_old
        # Report progress
        if mod(iter,250)==0
            println("   PFI Loop: iter=$iter, dist=",G_dist)
        end
        # Check convergence and return results
        if G_dist<=dist_tol
            println("PFI - n_ϵ=$n_ϵ, n_ζ=$n_ζ, n_a=$n_a - θ_a=$θ_a - r=$r")
            println("Iterations = $iter and Distance = ",G_dist)
            println("------------------------")
            println(" ")
            # Check borrowing constraint
            if any( G_ap_new.<p.a_min )
                error("Borrowing Constraint Violated")
            end
            # Interpolate to fine grid
            G_ap_fine = zeros(n_a_fine,n_ϵ,n_ζ)
            G_c_fine  = zeros(n_a_fine,n_ϵ,n_ζ)
            for i_ϵ=1:n_ϵ
            for i_ζ=1:n_ζ
            G_ap_ip = ScaledInterpolations(a_grid,G_ap_new[:,i_ϵ,i_ζ] , BSpline(Cubic(Line(OnGrid()))))
                G_ap_fine[:,i_ϵ,i_ζ].= G_ap_ip.(collect(a_grid_fine))
            G_c_ip  = ScaledInterpolations(a_grid,G_c[:,i_ϵ,i_ζ]  , BSpline(Cubic(Line(OnGrid()))))
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
        C_ip        = Spline1D(A_aux,C_aux)
        G_c[:,i_ϵ,i_ζ] .= C_ip.(M.a_grid)
        Ap_aux      = (r*M.ζ_mat[i_ζ].+1).*collect(M.a_grid) .+ w*M.ϵ_grid[i_ϵ] .- G_c[:,i_ϵ,i_ζ]
    end
    end
    # Update policy function
    G_ap .= (r*M.ζ_mat.+1).*M.a_mat .+ w*M.ϵ_mat .- G_c
        # Adjust for numerical error
        for ind = findall(<=(1e-10),abs.(G_ap.-a_min))
            G_ap[ind] = a_min
            G_c[ind]  = (r*M.ζ_mat[ind].+1)*M.a_mat[ind] + w*M.ϵ_mat[ind] - a_min
        end
        # Check for borrowing constraint
        #if any( G_ap.<a_min )
        #    error("Borrowing Constraint Violated")
        #end
    # Return Results
    return G_ap, G_c
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Histogram method
function Histogram_Method_Loop(M::Model,N_H=nothing,Γ_0=nothing)
    @unpack p, n_ϵ, n_ζ, MP_ϵ, MP_ζ, n_a_fine, a_grid_fine, θ_a, G_ap_fine, H_ind, H_weight = M
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
            H_ind[i_a,i_ϵ,i_ζ] = Grid_Inv(G_ap_fine[i_a,i_ϵ,i_ζ],n_a_fine,θ_a,a_min,a_max)
            ω_lo               = min(1,max(0,(G_ap_fine[i_a,i_ϵ,i_ζ]-a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]])/(a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]+1]-a_grid_fine[H_ind[i_a,i_ϵ,i_ζ]])))
        end

        # Sotre weights for lower and upper bounds on approxiamtion interval, including transition to future states 
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
            println("   Histogram Loop: iter=$i_H, dist=$H_dist")
        end
        # Check convergence
        if H_dist<Hist_tol
            println("Histogram iteartion converged in iteration $i_H. H_dist=$H_dist\n--------------------------------\n")
            M = Model(M; Γ=Γ, H_ind=H_ind, H_ω_lo=H_ω_lo, H_ω_hi=H_ω_hi)
            return M
        end
    end
end

# Graphs and Stats
function  Aiyagari_Graph(M::Model)
    gr()
    # Capital Policy Function
        plt = plot(title="Policy Function - a' - n_a=$(M.n_a) - θ_a=$(M.θ_a)",legend=:bottomright,foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(M.a_grid_fine,M.a_grid_fine,lw=1,linecolor=RGB(0.6,0.6,0.6),label=nothing)
        # Different values of ϵ and ζ
        m_ϵ = [1,convert(Int64,round(M.n_ϵ/2)),M.n_ϵ]
        m_ζ = [1,convert(Int64,round(M.n_ζ/2)),M.n_ζ]
        for i_ϵ=2:2
        plot!(M.a_grid,M.G_ap[:,m_ϵ[i_ϵ],m_ζ[2]],linetype=:scatter,ms=1.5,label="G_ap(ϵ_$m_ϵ[$i_ϵ])")
        plot!(M.a_grid_fine,M.G_ap_fine[:,m_ϵ[i_ϵ],m_ζ[2]],linewidth=1,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        end
        # for i_ζ=1:3
        # plot!(M.a_grid,M.G_ap[:,m_ϵ[2],m_ζ[i_ζ]],linetype=:scatter,ms=1.5,label="G_ap(ϵ_$m_ζ[$i_ζ])")
        # plot!(M.a_grid_fine,M.G_ap_fine[:,m_ϵ[2],m_ζ[i_ζ]],linewidth=1,linestyle=(:dash),linecolor=RGB(0.4,0.4,0.4),label=nothing)
        # end
        xlabel!("Assets")
        ylabel!("Assets")
        savefig("Figures/Aiyagari_G_ap.pdf")
end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Aiyagari Equilibrium
function Aiyagari_Equilibrium(M_in::Model)

    # Compute Policy Functions
    M_in   = PFI_Fixed_Point(T_EGM_G,M_in) ;

    # Compute Distribution
    #M  = Histogram_Method_Loop(M)

    # Graphs
    Aiyagari_Graph(M_in)

    M_in  = Histogram_Method_Loop(M_in) ;

    return M_in
end

println("\n===============================================\n Solving Aiyagari with EGM-Histogram(loop)")

@time M_Aiyagari = Aiyagari_Equilibrium(Model()); 

## Check time with for loop in T for PFI  
# kronecker product does better

## Add stop for PFI when step size is small 
# distance between convergence criterias
println("===============================================\n")


## Define marginal distributions 
    Γ_a = dropdims( sum( M_Aiyagari.Γ , dims=(3,2) ) , dims=(3,2) ) ; # Assets 
    Γ_ϵ = dropdims( sum( M_Aiyagari.Γ , dims=(1,3) ) , dims=(1,3) ) ; # Labor Efficiency 
    Γ_ζ = dropdims( sum( M_Aiyagari.Γ , dims=(1,2) ) , dims=(1,2) ) ; # Returns 

## Plot asset distribution 
    gr()
    scatter(M_Aiyagari.a_grid_fine,Γ_a,marker=(:circle ,3,:cornflowerblue  ),markerstrokewidth=0,label=nothing)   
    xlabel!("Assets (thousands of dollars)",labelsize=18)
    savefig("./"*Fig_Folder*"/Distribution_Wealth.pdf")

## Plot Saving Functions 
    gr() 
    plot(M_Aiyagari.a_grid_fine,M_Aiyagari.a_grid_fine,w=1,linecolor=:gray70,label=nothing,aspect_ratio=1,xlims=(M_Aiyagari.a_grid[1],M_Aiyagari.a_grid[end]))
    plot!(M_Aiyagari.a_grid_fine,M_Aiyagari.G_ap_fine[:,2,2],w=2,linecolor=:cornflowerblue,label=nothing,aspect_ratio=1)
    xlabel!("Assets (thousands of dollars)",labelsize=18)
    ylabel!("Assets (thousands of dollars)",labelsize=18)
    savefig("./"*Fig_Folder*"/Policy_Function_Savings_Level.pdf")

    gr()
    hline( [0] ,c=:gray70  ,w=1,label=nothing) 
    plot!(M_Aiyagari.a_grid_fine[25:end],100*(M_Aiyagari.G_ap_fine[25:end,2,2]./M_Aiyagari.a_grid_fine[25:end].-1),w=2,linecolor=:cornflowerblue,label=nothing,xlims=(M_Aiyagari.a_grid[1],M_Aiyagari.a_grid[end]))
    xlabel!("Assets (thousands of dollars)",labelsize=18)
    ylabel!("Saving Rate (%)",labelsize=18)
    savefig("./"*Fig_Folder*"/Policy_Function_Savings_Rate.pdf")

## Plots of policy function of assets (Choose high-median-low values of ϵ and ζ)


## Check time with for loop in T for PFI

## Add stop for PFI when step size is small 

## Verify moments of returns and income (expected values)



## [Optional Check] Euler Errors






println("\n===============================================\n\n    End of Script \n\n===============================================")