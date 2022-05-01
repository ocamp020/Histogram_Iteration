# Toolbox for VFI in Julia
# Sergio Ocampo
# September 2020

# module VFI_Toolbox

#-----------------------------------------------------------
#-----------------------------------------------------------
    using Parameters # Pkg.add("Parameters") # https://github.com/mauro3/Parameters.jl
    using Distributions #Pkg.add("Distributions")
    using LinearAlgebra
    using Random
    using Interpolations

#-----------------------------------------------------------
#-----------------------------------------------------------


#-----------------------------------------------------------
#-----------------------------------------------------------
# Cubic spline interpolation
#-----------------------------------------------------------
#-----------------------------------------------------------

#-----------------------------------------------------------
# Solve tri-diagonal system
    # Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) in NR
    # using a serial algorithm.
    # Input vectors b (diagonal elements) and r (right-hand sides) have size N,
    # while a and c (off-diagonal elements) are size N − 1.
function tridag(a,b,c,r)
    # Check that sizes agree
    n = length(a)+1 # Lenght of vectors
    if any( [length(b) length(c)+1 length(r)].!=n)
        error("Interpolation requires length(x)==length(y)")
    end
    # Check boundary conditions
    bet = b[1]
    if (bet == 0)
        error("tridag: Error at code stage 1")
        # If this happens then you should rewrite your equations as a set of order N − 1,
        # with u2 trivially eliminated.
    end
    # Solution of system
    u = zeros(n)
    g = zeros(n)
    u[1]=r[1]/bet
    for j=2:n
        g[j] = c[j-1]/bet
        bet  = b[j]-a[j-1]*g[j]
        if (bet == 0)
        error("tridag_ser: Error at code stage 2")
            # Decomposition and forward substitution.
            # Algorithm fails; see below routine in Vol. 1. of NR
        end
        u[j]=(r[j]-a[j-1]*u[j-1])/bet
    end
    for j=n-1:-1:1
        u[j]=u[j]-g[j+1]*u[j+1]
    end
    return u
end


#-----------------------------------------------------------
# Get second derivatives
    function spline_ypp(x,y,yp1=nothing,ypn=nothing)
        # Check that x grid is ordered
        if any(diff(x).<0)
            error("Grid for x must be increasing")
        end
        # Check that sizes agree
        if length(x)!=length(y)
            error("Interpolation requires length(x)==length(y)")
        end
        n = length(x) # Lenght of vectors
        # Set up tri-diagonal system
        a = Array{Float64}(undef,n)
        b = Array{Float64}(undef,n)
        c = Array{Float64}(undef,n)
        r = Array{Float64}(undef,n)
        # Fill in elements for tri-diagonal system
        c[1:n-1].= x[2:n].-x[1:n-1]
        r[1:n-1].= 6*((y[2:n].-y[1:n-1])./c[1:n-1])
        r[2:n-1].= r[2:n-1].-r[1:n-2]
        a[2:n-1].= c[1:n-2]
        b[2:n-1].= 2*(c[2:n-1].+a[2:n-1])
        b[1]     = 1
        b[n]     = 1
        # Lower Boundary
        if yp1==nothing # Use natural spline
            r[1] = 0
            a[1] = 0
        else # User supplied a derivative
            r[1] = (3/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1)
            c[1] = 0.5
        end
        # Upper Boundary
        if ypn==nothing # Use natural spline
            r[n] = 0
            c[n] = 0
        else # User supplied a derivative
            r[n] = (-3/(x[n]-x[n-1]))*((y[n]-y[n-1])/(x[n]-x[n-1])-ypn)
            a[n] = 0.5
        end
        # Compute second derivatives from tri-diagonal system
        ypp = tridag(a[2:n],b[1:n],c[1:n-1],r[1:n])
        return ypp
    end
#-----------------------------------------------------------

#-----------------------------------------------------------
# Get interpolation
    # Given the arrays x and y, which tabulate a function
    # (with the xi’s in increasing order),
    # and given the array ypp, which is the output from spline above,
    # and given a value of x, this routine returns a cubic-spline interpolated value.
    # The arrays x, y and ypp are all of the same size.
function spline_itp(x,y,ypp,z)
    # Check that z∈[x_min,x_max]
    if z<minimum(x) || z>maximum(x)
        error("Interpolation point z must be inside the grid bounds")
    end
    # Defien grid size and bracket z in grid
    n   = length(x) # Lenght of vectors
    khi = findmax(sign.(x .- z))[2] # Index of higher brakcet
    klo = khi-1                     # Index of lower braket
    h   = x[khi]-x[klo]             # Size of bracket
    # Define convexx weights
    a   = (x[khi]-z)/h
    b   = (z-x[klo])/h
    # Evaluate cubic spline
    itp = a*y[klo]+b*y[khi]+((a^3-a)*ypp[klo]+(b^3-b)*ypp[khi])*(h^2)/6
    return itp
end
#-----------------------------------------------------------

#-----------------------------------------------------------
# Get first derivative of interpolation
function spline_ditp(x,y,ypp,z)
    # Check that z∈[x_min,x_max]
    if z<minimum(x) || z>maximum(x)
        error("Interpolation point z must be inside the grid bounds")
    end
    # Defien grid size and bracket z in grid
    n   = length(x) # Lenght of vectors
    khi = findmax(sign.(x .- z))[2] # Index of higher brakcet
    klo = khi-1                     # Index of lower braket
    h   = x[khi]-x[klo]             # Size of bracket
    # Define convexx weights
    a   = (x[khi]-z)/h
    b   = (z-x[klo])/h
    # Evaluate cubic spline
    ditp = (y[khi]-y[klo])/(x[khi]-x[klo]) - ((3*a^2-1)*ypp[klo] + (3*b^2-1)*ypp[khi])*(x[khi]-x[klo])/6
    return ditp
end
#-----------------------------------------------------------

#-----------------------------------------------------------
# Wraper to get function
function spline_NR(x,y,yp1=nothing,ypn=nothing)
    # Get second derivatives
    ypp = spline_ypp(x,y,yp1,ypn)
    # Define interpolation function
    F(z)  = spline_itp(x,y,ypp,z)
    # Define derivative interpolation function
    dF(z) = spline_ditp(x,y,ypp,z)
    return F,dF
end
#-----------------------------------------------------------



#-----------------------------------------------------------
#-----------------------------------------------------------
# Markov Processes Block
#-----------------------------------------------------------
#-----------------------------------------------------------
    # MP: struct to hold Markov Processes
    # Tauchen86: Discretize AR(1) with Tauchen 1986 method
    # Rouwenhorst95: Discretize AR(1) with Rouwenhorst 1995 method
    # SigmaMatchR95: Finds a Sigma that will create a log process with a target std
    # SigmaLogR95: Find Std for an AR(1) process in logs with Rouwenhorst 1995
    # Simulation_MC: Simulate a Markov Chain from a Markov Process MP

#-----------------------------------------------------------
#-----------------------------------------------------------
# Defien a markov process struct
    # Generate structure for markov processes using Parameters module
    @with_kw struct MP
        # Model Parameters
        N::Int64 # Number of states
        grid     # Grid of discrete markov process
        Π        # Transition matrix
        PDF      # Stationary distribution
        CDF      # Stationary distribution
    end

#-----------------------------------------------------------
#-----------------------------------------------------------
# Tauchen (1986)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
        # Ω - Grid expansion in number of standard devaitions (Optional)
    # Outputs:
        # z - Grid of N equally spaced points covering [-Ωσ,Ωσ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Tauchen86(ρ,σ,N,Ω::Any=3)
    # Check if N>1
    if N==1
        z = 0; Π = 1; PDF_z = 1;
        return MP(N=N,grid=0,Π=1,PDF=1,CDF=1)
    end
    # Create z grid
        z = range(-Ω*σ/sqrt(1-ρ^2),Ω*σ/sqrt(1-ρ^2),length=N)
    # Define intermediate step length
        h = (z[2]-z[1])/2
    # Define auxiliary matrices
        z_0 = repeat(z ,1,N) # Matrix of today's z each row is a value, columns are equal
        z_1 = repeat(z',N,1) # Matrix of tomorrow's z each column is a value, rows are equal
    # Define intervals
        z_lim = zeros(N,N,2) # First matrix is lower bounds. Second matrix is uppor bounds.
        z_lim[:,1      ,1] .= -Inf
        z_lim[:,2:end  ,1] .=  ( z_1[:,2:end  ] - ρ*z_0[:,2:end  ] .- h )./σ
        z_lim[:,1:end-1,2] .=  ( z_1[:,1:end-1] - ρ*z_0[:,1:end-1] .+ h )./σ
        z_lim[:,end    ,2] .=  Inf
    # Define reference distribution
        # This line uses "Distributions"
        F(x) = cdf.(Normal(),x)
    # Fill in transition matrix
        Π_z = F.(z_lim[:,:,2]) - F.(z_lim[:,:,1])
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    # Get stationary distribution of markov chain
        PDF_z = real(eigvecs(Π_z')[:,end]); PDF_z = PDF_z/sum(PDF_z) ;
        CDF_z = cumsum(PDF_z)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Rouwenhorst (1995)
    # Objective is to discretize AR(1) process: z'=ρz+η, η~N(0,σ)
    # Code from Kopecky & Suen (2010)
    # Inputs:
        # ρ - Process persisntence
        # σ - Innovation standard deviation
        # N - Size of the grid
    # Outputs:
        # z - Grid of N equally spaced points covering [-ψ,ψ]
        # Π - Transition matrix, a stochastic matrix (sums to 1 across columns)
        # PDF_z, CDF_z - Stationary distribution of z
function Rouwenhorst95(ρ,σ,N)
    # Check if N>1
    if N==1
        z = 0; Π = 1; PDF_z = 1;
        return MP(N=N,grid=0,Π=1,PDF=1,CDF=1)
    end
    # Define paramters for Rouwenhorst's approximation
        p = (1+ρ)/2
        q = p                   # Note: I am leaving q here for comparability with source
        ψ = σ*sqrt((N-1)/(1-ρ^2))
        s = (1-q)/(2-(p+q))     # Note: s=0.5, I leave it for comparability with source
    # Fill in transition matrix
    if N==2
        Π_z = [p 1-p ; 1-q q]
    else
        MP_aux = Rouwenhorst95(ρ,σ,N-1)
        o = zeros(N-1)
        Π_z = p*[MP_aux.Π o ; o' 0] + (1-p)*[o MP_aux.Π ; 0 o'] + (1-q)*[o' 0 ; MP_aux.Π o] + q*[0 o' ; o MP_aux.Π]
        # Adjust scale for double counting
        Π_z = Π_z./repeat(sum(Π_z,dims=2),1,N)
    end
    # Distribution
        PDF_z = pdf.(Binomial(N-1,1-s),(0:N-1))
        CDF_z = cumsum(PDF_z)
    # Create z grid
        z    = range(-ψ,ψ,length=N)
    # Return
        return MP(N=N,grid=z,Π=Π_z,PDF=PDF_z,CDF=CDF_z)
end



#-----------------------------------------------------------
#-----------------------------------------------------------
# SigmaMatchR95
    # Objective is to find sigma_z so that the process:
    #   x=exp(z) where  z'=ρz+η, η~N(0,sigma_z)
    # has a given sigma_x
    # Inputs:
        # mean_x
        # ρ_z - Process persisntence
        # sigma_x - Target standard deviation for x
        # N - Size of the grid
    # Outputs:
        # sigma_z - Standard deviation of z
        function SigmaMatchR95(mean_x,ρ,sigma_x,N)
            # Check if N>1
            if N==1
                sigma_z=0
                return sigma_z
            end
            
            # Find Lower Bound
            sLB=sqrt(log(1+sigma_x^2*(1-ρ^2)^2/mean_x^2))
            fLB=SigmaLogR95(mean_x,ρ,sLB,N)
            if (fLB>sigma_x) 
                println("Error in SigmaMatchR95")
                println("LB=$sLB   gives fLb=$fLB>$sigma_x")
                sigma_z=sLB
                return sigma_z
            end

            
            # Find Upper Bound
            sUB=sLB*2
            fUB=SigmaLogR95(mean_x,ρ,sUB,N)
            if (fUB<sigma_x)
                println("Error in SigmaMatchR95")
                println("UB=$sUB   gives fUb=$fUB<$sigma_x")
                sigma_z=sUB
                return sigma_z
            end 

            # Bracket
            iter=0
            dist=1
            while (iter<50) & (dist>1E-6)
                iter=iter+1

                global sMP=(sLB+sUB)/2
                fMP=SigmaLogR95(mean_x,ρ,sMP,N)

                if (fMP>sigma_x) 
                    sUB=sMP
                    fUB=fMP
                else
                    sLB=sMP
                    fLB=fMP
                end 

                dist=abs(sLB-sUB)
            end
            sigma_z=sMP
            # Return
                return sigma_z
        end
#-----------------------------------------------------------
#-----------------------------------------------------------
# SigmaLogR95
    # Calculates the standard deviation of x where:
    #   x=exp(z) where  z'=ρz+η, η~N(0,sigma_z)
    # Inputs:
        # mean_x
        # ρ_z - Process persisntence
        # sigma_z - Standard deviation for z
        # N - Size of the grid
    # Outputs:
        # sigma_x - Standard deviation of x
        function SigmaLogR95(mean_x,ρ,sigma_z,N)

            # Generate Rouwenhorst Process
            MP_z=Rouwenhorst95(ρ,sigma_z,N)

            #Renomalize Mean
            x_grid=exp.(MP_z.grid)*mean_x/sum(exp.(MP_z.grid).*MP_z.PDF)
            
            #Caluculate Std
            sigma_x=sqrt(sum(MP_z.PDF.*(x_grid.-mean_x).^2))

            # Return
            return sigma_x
        end
        

#-----------------------------------------------------------
#-----------------------------------------------------------
# Simulation of Markov processes

# Simulation function for discrete Markov process
# The result of the simulation is a Markov chain
    # Inputs:
        # Ns - Number of simulated periods
        # z  - Grid with levels of Markov process
        # Π  - Transition matrix of Markov process
        # N_0- Number of periods to throw before reporting simulation (Optional)
    # Output:
        # z_MC - Vector of Ns elemenrts with Markov Chain
function Simulation_MC(Ns,MP::MP,N_0=1000)
    # Compute conditional CDF
    Γ = cumsum(MP.Π,dims=2)
    # Allocate simulated vector
    z_ind    = zeros(Int64,N_0+Ns)
    z_MC     = zeros(N_0+Ns)
    # Starting value for simulation
    z_ind[1] = Int(ceil(length(MP.grid)/2))
    z_MC[1]  = MP.grid[z_ind[1]]
    # Simulate
    for i=2:Ns+N_0
        #= Option 1
        # Draw a uniform random number (r). Compare it with conditional CDF.
        # Conditional CDF given by cumsum of current row of Π, call it Γ
        # Γ[i,j] gives the probability z'<z_j
        # We are looking for j s.t Γ[i,j-1]<r<Γ[i,j]
        # Equivalently, the lowest j s.t. Γ[i,j]-r>0
            z_ind[i] = findmax(sign.(Γ[z_ind[i-1],:] .- rand()))[2]
        =#
        #= Option 2
        # Alternatively we can draw directly from conditional distributional
        # Distributional is categorical with P[z'=z_j] = Π[i,j]
        =#
        z_ind[i] = rand(Categorical(MP.Π[z_ind[i-1],:]))
        z_MC[i]  = MP.grid[z_ind[i]]
    end
    # Throw out first N_0 elements
    z_MC = z_MC[N_0+1:end]
    # Return result
    return z_MC
end
#-----------------------------------------------------------
#-----------------------------------------------------------
# End of Markov Processes Block
#-----------------------------------------------------------
#-----------------------------------------------------------


################################################################################
#################### Scaled Grid Interpolation #################################
################################################################################
    # The following code generates three objects:
        # 1) A "TypeScale" type that transforms ranges into a type of scaled ranges
        # 2) A "TypeRange" type that foroms scaled ranges of a given type.
        #   These ranges behave "just as" normal ranges
        # 3) A ScaledInterpolations function that wraps "interpolations.jl" to be used
        #   with TypeRange ranges
    # Types of Scales: parameters bounds for grid (a,b) and curvature (θ)
        # p = PolyScale(a,b,θ)
            # Maps a number x∈[0,1] to a polynomial grid with curvature θ: p(x)∈[a,b]
            # Usage: p(x), where x is a number, or p.(x) where x is vector/range
            # p(x) = a + (b - a) * x^θ
        # p = InvPolyScale(a,b,θ)
            # Maps a number y∈[a,b] from a polynomial grid with curvature θ to the unit interval: p(x)∈[0,1]
            # Usage: p(y), where x is a number, or p.(y) where y is vector/range
            # p(y) = ((y-a) / (b-a) )^(1/θ)
        # p = ExpScale(a,b,θ)
            # Maps a number x∈[0,1] to a exponential grid with curvature θ: p(x)∈[a,b]
            # Usage: p(x), where x is a number, or p.(x) where x is vector/range
            # p(x) = a + (b - a) * ( (exp(θ*x) - 1)/(exp(θ) - 1) )
        # p = LogScale(a,b,θ)
            # Maps a number y∈[a,b] from an exponential grid with curvature θ to the unit interval: p(x)∈[0,1]
            # Usage: p(y), where x is a number, or p.(y) where y is vector/range
            # p(y) = log( (y-a)*(exp(θ) - 1)/(b - a) + 1 ) / θ
    # Types of Scaled Ranges
        # grid = PolyRange(a,b,θ=t,N=n) = PolyScale(a,b,t)(range(0,1,length=n))
        # grid =  ExpRange(a,b,θ=t,N=n) =  ExpScale(a,b,t)(range(0,1,length=n))
    # ScaledInterpolations
        # itp = ScaledInterpolations(grid,f_grid,Itp_Type)
            # grid is either a PolyRange or ExpRange type
            # f_grid is a vector with the values of the function at the grid nodes
            # Itp_Type is a type from the interpolations.jl package
                # Examples: BSpline(Cubic(Line(OnGrid()))) // shortcut: CubicSpline
                #           BSpline(Cubic(Flat(OnGrid())))
                #           FritschButlandMonotonicInterpolation() // shortcut: MonotoneSpline
                #           BSpline(Linear()) // shortcut: LinearInterp
    # ScaledInterpolations with multiple dimensions
        # itp_md = ScaledInterpolations( (grid_1,...,grid_N) , (f_grid_1,...,f_grid_N) , (Type_1,...,Type_N)  ) 
    # Jacob Adenbaum, 2020


    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Define Scaled types
        #-------------------------------------------------------------------------------
        # Polynomial Scaling
        abstract type Scale{T}       end
        abstract type ScaledRange{T} end

        for Foo in [:PolyScale, :InvPolyScale]

            @eval struct $Foo{T,TF} <: Scale{T}
                a::T
                b::T
                θ::TF
            end

            @eval function $Foo(a, b, θ = 1)
                T = reduce(promote_type, map(typeof, (a,b)))
                $Foo(convert(T, a), convert(T, b), θ)
            end
        end


        (p::PolyScale)(x)           = p.a + (p.b - p.a) * x^p.θ
        (p::InvPolyScale)(y)        = ((y-p.a) / (p.b-p.a) )^(1/p.θ)
        Base.inv(p::PolyScale)      = InvPolyScale(p.a, p.b, p.θ)
        Base.inv(p::InvPolyScale)   = PolyScale(p.a, p.b, p.θ)

        #-------------------------------------------------------------------------------
        # Exponential Scaling
        for Foo in [:ExpScale, :LogScale]

            @eval struct $Foo{T, TF} <: Scale{T}
                a::T
                b::T
                θ::TF
                et::T
                s::T
            end

            @eval function $Foo(a, b, θ = 1)
                s = (exp(θ) - 1)/(b - a)
                et= exp(θ)
                T = reduce(promote_type, map(typeof, (a,b,θ, et,s)))
                $Foo(convert(T, a), convert(T, b), θ, et,s)
            end
        end

        function (p::ExpScale)(x)
            t  = (exp(p.θ * x) - 1)/(p.et - 1)
            return p.a + (p.b - p.a) * t
        end

        function (p::LogScale)(y)
            θ  = p.θ
            log( (y-p.a)*p.s + 1 ) / p.θ
        end

        Base.inv(p::ExpScale)       = LogScale(p.a, p.b, p.θ)
        Base.inv(p::LogScale)       = ExpScale(p.a, p.b, p.θ)

        #-------------------------------------------------------------------------------
        # Finalize Types
        scales = Dict(:PolyRange => :PolyScale, :ExpRange => :ExpScale)


        Base.inv(::Type{ExpScale}) = LogScale
        Base.inv(::Type{LogScale}) = ExpScale
        Base.inv(::Type{PolyScale})= InvPolyScale
        Base.inv(::Type{InvPolyScale}) = PolyScale

    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Define Scaled Ranges types
        #-------------------------------------------------------------------------------
        # Polynomial Range
        for Typ in [:PolyRange] @eval begin
            struct $Typ{T,TF} <: AbstractRange{T}
                a::T
                b::T
                θ::TF
                N::Int
                vals::Vector{Float64}
                scale::$(scales[Typ]){T,TF}
                invscale::inv($(scales[Typ])){T,TF}
            end

            function $Typ(a,b; θ = 1, N = 100)
                T       = reduce(promote_type, map(typeof, (a,b,θ, 1.0)))
                Ta      = convert(T, a)
                Tb      = convert(T, b)

                vals    = domscale($Typ)(Ta,Tb,θ).(LinRange(0, 1, N))
                scale   = $(scales[Typ])(Ta, Tb, θ)
                invscale= inv($(scales[Typ])(Ta, Tb, θ))
                $Typ(convert(T, a), convert(T, b), θ, N, vals, scale, invscale)
            end

            domscale(p::$Typ)      = p.scale
            domscale(::Type{$Typ}) = $(scales[Typ])

            function rangescale(r::$Typ, x)
                r.invscale(x)*(length(r)-1) + 1
            end

            Base.getindex(p::$Typ, i::Int) = p.vals[i]
            Base.length(p::$Typ)  = p.N
            Base.show(io::IO, p::$Typ) = begin
                Typ = $Typ
                print(io, "$Typ($(p.a), $(p.b), $(p.θ), $(p.N))")
            end

            function Base.searchsortedfirst(p::$Typ, x, args...)
                searchsortedfirst(p.vals, x, args...)
            end


        end end

        #-------------------------------------------------------------------------------
        # Exponential Range
        for Typ in [:ExpRange] @eval begin
            struct $Typ{T,TF} <: AbstractRange{T}
                a::T
                b::T
                θ::TF
                N::Int
                vals::Vector{Float64}
                scale::$(scales[Typ]){T,TF}
                invscale::inv($(scales[Typ])){T,TF}
            end

            function $Typ(a,b; θ = 1, N = 100)
                vals    = domscale($Typ)(a,b,θ).(LinRange(0, 1, N))
                scale   = $(scales[Typ])(a, b, θ)
                invscale= inv($(scales[Typ])(a, b, θ))
                T       = reduce(promote_type, map(typeof, (a,b,θ,scale.s)))
                $Typ(convert(T, a), convert(T, b), θ, N, vals, scale, invscale)
            end

            domscale(p::$Typ)      = p.scale
            domscale(::Type{$Typ}) = $(scales[Typ])

            function rangescale(r::$Typ, x)
                r.invscale(x)*(length(r)-1) + 1
            end

            Base.getindex(p::$Typ, i::Int) = p.vals[i]
            Base.length(p::$Typ)  = p.N
            Base.show(io::IO, p::$Typ) = begin
                Typ = $Typ
                print(io, "$Typ($(p.a), $(p.b), $(p.θ), $(p.N))")
            end

            function Base.searchsortedfirst(p::$Typ, x, args...)
                searchsortedfirst(p.vals, x, args...)
            end
        end end

        #-------------------------------------------------------------------------------
        # Finalize Ranges
        function rangescale(r::AbstractRange, x)
            ((x - first(r)))/(last(r) - first(r)) * (length(r) -1) + 1
        end

        rangescale(r::UnitRange, x::Int) = x


    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Scaled Interpolation Wraper
        mutable struct ScaledInterpolations{T, N, IT <: Interpolations.AbstractInterpolation,
                                    RT <: NTuple{N, AbstractRange},
                                    VT <: AbstractArray{T, N}, AT}
            r::RT       # Interpolation grid (ranges)
            v::VT       # Interpolation values
            itp::IT     # Interpolation object
            args::AT    # Arguments to pass to construct interpolant
        end

        function ScaledInterpolations(r, v, args...)
            itp = extrapolate(interpolate(v, args...), Interpolations.Flat())
            return ScaledInterpolations(tuplefy(r), v, itp, args)
        end

        tuplefy(x::Tuple) = x
        tuplefy(x)        = tuple(x)

        """
        ```
        fit!(s::ScaledInterpolations, [v])
        ```
        Re-fit the interpolation after the underlying data has updated.  Works well if
        you use a view into the original array.  If an array of values is passed, copies
        the values to the scaled interpolations value array and then refits the
        interpolation.
        """
        function fit!(s::ScaledInterpolations)
            s.itp = extrapolate(interpolate(s.v, s.args...), Interpolations.Flat())
        end

        function fit!(s::Array{T,N}) where {T <: ScaledInterpolations, N}
            for v in s
                fit!(v)
            end
        end

        dim(::Type{ScaledInterpolations{T, N, IT, RT, VT, AT}}) where {T, N, IT, RT, VT, AT} = N
        Base.eltype(::Type{ScaledInterpolations{T, N, IT, RT, VT, AT}}) where {T, N, IT, RT, VT, AT} = T
        dim(itp::ScaledInterpolations) = dim(typeof(itp))
        Base.eltype(itp::ScaledInterpolations) = eltype(typeof(itp))

        @generated function (s::ScaledInterpolations)(x::Vararg)
            Ns = dim(s)
            Nx = length(x)
            Ns == Nx || begin
                return quote
                    Ns = $Ns
                    Nx = $Nx
                    throw(ArgumentError("Must have $Ns arguments.  You passed $Nx"))
                end
            end

            ex = Expr(:call, :(s.itp))
            for i = 1:Ns
                push!(ex.args, :(rangescale(s.r[$i], x[$i])))
            end
            ex
        end

    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Scaled Interpolation Shortcuts for Itp_Type
        const CubicSpline    = BSpline(Cubic(Natural(OnGrid())))
        const MonotoneSpline = FritschButlandMonotonicInterpolation()
        const LinearInterp   = BSpline(Linear())
        const LinearInterp2  = BSpline(Linear())

        function ScaledInterpolations(grid, values, type::FritschButlandMonotonicInterpolation)
            return interpolate(collect(grid), values, type)
        end
################################################################################
################ End of Scaled Grid Interpolation ##############################
################################################################################


#-----------------------------------------------------------
#-----------------------------------------------------------
# Miscelaneous Functions
#-----------------------------------------------------------
#-----------------------------------------------------------
    # Make_Grid: makes curved grid, uses ScaledInterpolation
    # Grid_Inv: gets index of closest node in grid
    # mnbrak: Bracket a minimum for a one-diimensional function expanding an orignal interval [a,b]

#-----------------------------------------------------------
#-----------------------------------------------------------
# Grid
function Make_Grid(n,θ_x,x_min,x_max,scale_type="Poly")
    # Get x_grid
    if θ_x≠1
        if scale_type=="Poly"
        x_grid = PolyRange(x_min,x_max;θ=θ_x,N=n) ; # Curved grid between x_min and x_max
        elseif scale_type=="Exp"
        x_grid = ExpRange(x_min,x_max;θ=θ_x,N=n) ; # Curved grid between x_min and x_max
        else
        error("scale_type must be either Poly or Exp")
        end
    else
    x_grid = range(x_min,x_max,length=n)
    end
    # Return
    return x_grid
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Grid_Inv
function Grid_Inv(x,n,θ_x,x_min,x_max,scale_type="Poly")

    # Check corner solution
    if x<x_min
        x_ind = 1
    elseif x>x_max
        x_ind = n
    else
    # Get index of x in grid if x∈[x_min,x_max]
        if θ_x≠1
            if scale_type=="Poly"
            x_ind  = Int(floor(((x.-x_min)/(x_max-x_min)).^(1/θ_x) *(n-1)+1))
            elseif scale_type=="Exp"
            x_ind  = Int(floor( (1/θ_x).*log.(((x-x_min)*(exp(θ_x)-1)/(x_max-x_min)).+1) *(n-1)+1))
            else
            error("scale_type must be either Poly or Exp")
            end
        else
        x_ind = Int(floor(((x-x_min)/(x_max-x_min))*(n-1)+1))
        end
    end
    # Return
    return x_ind
end


#-----------------------------------------------------------
#-----------------------------------------------------------
# Minimum Bracketing Function
    # This function comes from Numerical Recipes, Section 10.1
    # Input: numbers a,b that initiliaze bracket and objective function F
        # Optional: Bounds for the brackeet so that a,b,c ∈ [x_min,x_max]
    # Output: numbers (a,b,c) and images (fa,fb,fc) such that a<b<c and fb<fa,fc
    function mnbrak(a,b,F::Function,x_min=-Inf,x_max=Inf)
        # Auxiliariy variables
        Tiny    = 1E-20
        G_limit = 100
        # Define golden ratio
        γ = 1.618034
        # Evaluate function at end points
        fa, fb = F(a), F(b)
        # Swap a and b so that we can go downhill in the direction from a to b.
            # This way we know for certain that fb<fa, we only need to find c such that fb<fc
        if fb>fa
        a , b = b , a
        fa,fb = fb, fa
        end
        # Guess for value c expanding bracket away from b -> (a,b,c) or (c,b,a)
            # Check bounds
        c  = max( min( b + γ*(b-a) , x_max ) , x_min)
        fc = F(c)
        # While fb>fc we need to keep on bracketing
        while fb>fc
            # Compute u (new candidate) by parabolic extrapolation from a, b, c.
                # TINY is used to prevent any possible division by zero.
            r = (b-a)*(fb-fc)
            q = (b-c)*(fb-fa)
            u = min(max(b-((b-c)*q-(b-a)*r)/(2*sign(q-r)*max(abs(q-r),Tiny)),x_min),x_max)
            u_lim = min(max(b+G_limit*(c-b),x_min),x_max)
            # Test various cases for new candidate
            if ((b-u)*(u-c) > 0) # Parabolic u is between b and c
                fu=F(u)
                if (fu < fc) # Got a minimum between b and c so bounds are (b,u,c)
                    a, fa, b, fb = b, fb, u, fu
                    break
                elseif (fu > fb) # Got a minimum between a and u so bounds are (a,b,u)
                    c, fc = u, fu
                    break
                end
                # If you got here is because candidate u failed
                # Get new candidate by expanding interval with golden ratio
                # Check bounds
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            elseif ((c-u)*(u-u_lim) > 0.0) # Parabolic u is between c and its limit (ulim)
                fu=F(u)
                if (fu < fc) # Drop c and replace it with u, get new u with golden expansion
                    b, c, fb, fc = c, u, fc, fu
                    u  = max(min( c+γ*(c-b) , x_max ),x_min)
                    fu = F(u)
                end
            elseif ((u-u_lim)*(u_lim-c) >= 0.0) # U is above its limit, reign it in!
                u  = u_lim
                fu = F(u)
            else # Nothing worked, use golden expansion
                u  = max(min( c+γ*(c-b) , x_max ),x_min)
                fu = F(u)
            end
            # If no break then forget the oldest point and go onto next iteration
                # This means assigning b->a, c->b, u-> and forgetting about a
            a, b, c, fa, fb, fc = b, c, u, fb, fc, fu
        end
        # Return solution once out of the loop
            # Swap a and c so that a<b<c
            if a>c
            a , c  = c, a
            fa, fc = fc, fa
            end
        # Minimum bracketed in [a,c] with intermediate point b so that fb<fa,fc
        # println("Minimum bracketed in [$a,$c] with intermediate point b=$b \n Function values: F(a)=$fa, F(b)=$fb, F(c)=$fc")
        return a,c,b,fa,fc,fb
    end


# end
