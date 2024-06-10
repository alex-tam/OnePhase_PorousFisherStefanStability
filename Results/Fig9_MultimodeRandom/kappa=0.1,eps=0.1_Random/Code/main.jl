# Control function for 2D Porous-Fisher-Stefan level-set solutions
# Nizhum Rahman and Alex Tam, 08/12/2023

# Load packages
using Parameters
using Printf
using Dierckx
using LinearAlgebra
using DifferentialEquations
using Measures
using LaTeXStrings
using DelimitedFiles
using Roots
using Distributions

# Include external files
include("twic.jl")
include("domain.jl")
include("porous-fisher.jl")
include("velocity_extension.jl")
include("interface_density.jl")
include("interface_speed.jl")
include("level-set.jl")
include("reinitialisation.jl")

############################
##### Model parameters #####
############################
"Data structure for parameters"
@with_kw struct Params
    D::Float64 = 1.0 # [-] Diffusion coefficient
    m::Float64 = 1.0 # [-] Nonlinear diffusion exponent, D(u) = u^m
    λ::Float64 = 1.0 # [-] Reaction rate
    κ::Float64 = 0.1 # [-] Inverse Stefan number
    γ::Float64 = 0.0 # [-] Surface tension coefficient
    β::Float64 = 5.0 # [-] Initial interface position
    uf::Float64 = 0.0 # [-] Background density at interface
    θb::Float64 = 0.01 # [-] Threshold for whether a grid point is close to interface (relative to Δx)
    θ::Float64 = 1.99 # [-] Parameter for minmod flux-limiter
    Lx::Float64 = 10.0 # [-] Spatial domain limit (x)
    Ly::Float64 = 10.0 # [-] Spatial domain limit (y)
    Lξ::Float64 = 10.0 # [-] Domain width for travelling wave (ξ)
    T::Float64 = 10.0 # [-] End time
    Nx::Int = 401 # [-] Number of grid points (x)
    Ny::Int = 401 # [-] Number of grid points (y)
    Nt::Int = 10001 # [-] Number of time steps
    plot_interval::Int = 1000 # [-] Time steps between successive output
    Nξ::Int = 101 # [-] Number of grid points for travelling wave (ξ)
    a::Float64 = 1e-2 # [-] Parameter for geometric progression
    V_Iterations::Int = 20 # [-] Number of iterations for velocity extrapolation PDE
    ϕ_Iterations::Int = 20 # [-] Number of iterations for reinitialisation PDE
    ε::Float64 = 0.1 # [-] Small amplitude of perturbations
end

#################################
##### Generate perturbation #####
#################################
"Function to generate a random periodic perturbation"
function perturbation(par::Params)
    d = Uniform(-1.0, 1.0) # Uniform distribution with mean zero
    len::Int = (par.Ny + 1)/2 # Half-length for random perturbation
    pl = rand(d, len) # Random perturbation for half-width
    pls = movingaverage(pl, 5) # Smooth
    plsn = pls./maximum(abs.(pls)) # Normalise perturbation amplitude
    return vcat(plsn, reverse(plsn[1:end-1])) # Reflect random perturbation to ensure periodicity
end

"Function to generate a multi-mode periodic perturbation"
function perturbation(y::StepRangeLen)
    q1 = 2*pi/10 # Mode 1
    q2 = 4*pi/10 # Mode 2
    p = cos.(q1.*y) .+ cos.(q2.*y)
    pn = p./maximum(abs.(p)) # Normalise perturbation amplitude
    return pn
end

"Function to generate a random perturbed correction term"
function correction(ξ)
    d = Uniform(-1.0, 1.0) # Uniform distribution with mean zero
    p = rand(d, length(ξ)) # Random perturbation
    p[1] = 0.0; p[end] = 0.0 # Boundary conditions for u1
    ps = movingaverage(p, 5) # Smooth
    psn = ps./maximum(abs.(ps)) # Normalise perturbation amplitude
    return psn
end

"Calculate a moving average"
function movingaverage(x::Vector, numofele::Int)
    BackDelta = div(numofele,2)
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(x)
    y = similar(x)
    for n = 1:len
        lo = max(1, n - BackDelta)
        hi = min(len, n + ForwardDelta)
        y[n] = mean(x[lo:hi])
    end
    return y
end

####################################
##### Obtain initial condition #####
####################################
"Interpolate to obtain initial condition"
function ic(par, x, y)
    U = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of U
    ϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of ϕ
    u0, ξ = twic(par)
    # Construct linear splines for interpolation
    spl_0 = Spline1D(ξ, u0; k=1) # Generate 1D linear spline
    # Generate perturbation (correction)
    u1 = correction(ξ)
    spl_1 = Spline1D(ξ, u1; k=1) # Linear spline for correction u_1(ξ)
    # Generate perturbation (transverse)
    p = perturbation(par) # Random perturbation p(y)
    # p = perturbation(y) # Multi-mode perturbation p(y)
    ##### Obtain initial front (level-set function) #####
    for i in eachindex(x)
        for j in eachindex(y)
            # ϕ[i,j] = x[i] - par.β # Unperturbed front position
            ϕ[i,j] = x[i] - par.β - par.ε.*p[j] # Perturbed front position
        end
    end
    ##### Obtain initial densities ####
    for i in eachindex(x)
        for j in eachindex(y)
            if ϕ[i,j] < 0 # If grid point is in Ω(0)
                # U[i,j] = spl_0(ϕ[i,j]) # Unperturbed density
                U[i,j] = spl_0(ϕ[i,j]) + par.ε.*spl_1(ϕ[i,j]).*p[j] # Perturbed density
            else # If grid point is not in Ω(0)
                U[i,j] = par.uf
            end
        end
    end
    return U, ϕ
end

###############################
##### Reshaping functions #####
###############################
"Build vector from matrix, ordered by entries in D"
function build_vector(U::Array{Float64}, D)
    u = Vector{Float64}() # Pre-allocate empty vector
    for gp in D
        push!(u, U[gp.xInd, gp.yInd])
    end
    return u
end

"Build matrix from vector ordered by entries in D"
function build_u_matrix(u::Vector, y, par, D)
    U = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on right boundary)
    U[par.Nx,:] .= par.uf # Dirichlet condition on right boundary
    U[1,:] .= 1.0 # Dirichlet condition on left boundary
    for i in eachindex(D)
        U[D[i].xInd, D[i].yInd] = u[i]
    end
    return U
end

"Build matrix from vector ordered by entries in D"
function build_v_matrix(v::Vector, par, D)
    V = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on computational boundary)
    for i in eachindex(D)
        V[D[i].xInd, D[i].yInd] = v[i]
    end
    return V
end

##############################
##### Inteface detection #####
##############################
"Construct y-dependent interface from level-set function"
function interface(ϕ, x, y, dx)
    # Pre-allocate
    L = Vector{Float64}(undef, length(y))
    # Loop over discrete values of y to find x-position of interface
    for j in eachindex(y)
        ϕx = ϕ[:, j] # Obtain 1D vector of ϕ
        for i in eachindex(x)
            if (ϕx[i] < 0) && (ϕx[i+1] >= 0)
                θ = ϕx[i]/(ϕx[i] - ϕx[i+1]) # Linear interpolation
                L[j] = x[i] + θ*dx
            end
        end
    end
    return mean(L), (maximum(L)-minimum(L))/2
end

#########################
##### Main function #####
#########################
"Compute a solution"
function porous_fisher_stefan_2d()
    ##### Parameters and domain #####
    par = Params() # Initialise data structure of model parameters
    nx::Int = (par.Nx-1)/2 # Indices for slice plots (x)
    ny::Int = (par.Ny-1)/2 # Indices for slice plots (y)
    x = range(0, par.Lx, length = par.Nx); dx = x[2] - x[1] # Computational domain (x)
    y = range(0, par.Ly, length = par.Ny); dy = y[2] - y[1] # Computational domain (y)
    t = range(0, par.T, length = par.Nt); dt = t[2] - t[1] # Time domain
    writedlm("x.csv", x); writedlm("y.csv", y); writedlm("t.csv", t) # Write data to files
    ##### Initial condition #####
    U, ϕ = ic(par, x, y) # Obtain initial density and ϕ
    Φ = U.^(par.m + 1) # Φ(x,y,t)
    writedlm("U-0.csv", Φ)
    writedlm("Phi-0.csv", ϕ) # Write data to files
    plot_times = Vector{Int}() # Vector of time-steps at which data is obtained
    writedlm("plot_times.csv", plot_times)
    L = Vector{Float64}() # Preallocate empty vector of interface position
    Amp = Vector{Float64}() # Preallocate empty vector of perturbation amplitude
    Li, amp = interface(ϕ, x, y, dx) # Mean front position and perturbation amplitude of the initial condition
    push!(L, Li); push!(Amp, amp)
    ##### Time stepping #####
    for i = 1:par.Nt-1
        # 1. Find Ω, dΩ, and irregular grid points
        D = find_domain(par, ϕ)
        dΩ = find_interface(par, D, ϕ)
        # 2. Solve Porous-Fisher equation on Ω
        uf = interface_density(dΩ, ϕ, par, dx, dy) # Density on interface for BC
        @time Φ = pf(D, dΩ, Φ, ϕ, uf, y, par, dx, dy, dt, i)
        # 3. Compute extension velocity field
        V = extend_velocity(D, dΩ, Φ, ϕ, par, dx, dy)
        # 4. Solve level-set equation
        ϕ = level_set(V, ϕ, par, dx, dy, dt)
        # 5. Re-initialise level-set function as a signed-distance function
        if mod(i, 1) == 0
            ϕ = reinitialisation(ϕ, par, dx, dy, par.ϕ_Iterations)
        end
        # 6. (Optional) Post-processing
        if mod(i, par.plot_interval) == 0
            writedlm("ux-$i.csv", Φ[:,ny])
            writedlm("uy-$i.csv", Φ[nx,:])
            writedlm("U-$i.csv", Φ)
            writedlm("V-$i.csv", V)
            writedlm("Phi-$i.csv", ϕ)
            push!(plot_times, i)
            writedlm("plot_times.csv", plot_times)
        end
        Li, amp = interface(ϕ, x, y, dx) # Mean front position and perturbation amplitude of the solution
        push!(L, Li)
        push!(Amp, amp)
        writedlm("L.csv", L)
        writedlm("Amp.csv", Amp)
        writedlm("m.csv", par.m)
    end
end

@time porous_fisher_stefan_2d()