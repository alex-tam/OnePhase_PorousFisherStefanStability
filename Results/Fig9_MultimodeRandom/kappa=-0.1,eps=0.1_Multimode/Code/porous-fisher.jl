# Solve 2D Porous-Fisher-Stefan model with  Dirichlet (x) and periodic (y) boundary conditions
# Nizhum Rahman and  Alex Tam, 08/08/2023

"Solve Porous-Fisher equation"
function pf(D, dΩ, U, ϕ, uf, y, par, dx, dy, dt, i)
    # Create vector from matrix data
    u = build_vector(U, D)
    # Construct and solve ODE problem using DifferentialEquations.jl
    prob = ODEProblem((du, u, p, t) -> pf_rhs!(du, u, p, t, D, dΩ, ϕ, uf, y, par, dx, dy), u, ((i-1)*dt, i*dt))
    sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-16, saveat = i*dt)
    # Reshape solution to matrix
    U = build_u_matrix(sol[:,end], y, par, D)
    return U
end

"Construct right-hand vector for use in DifferentialEquations.jl"
function pf_rhs!(du, u, p, t, D, dΩ, ϕ, uf, y, par, dx, dy)
    U = build_u_matrix(u, y, par, D) # Generate matrix
    for i in eachindex(D) # Loop over grid points
        if (D[i].yInd != 1) && (D[i].yInd != par.Ny) && (D[i].Ω == true) && (D[i].dΩ == false) # Interior grid points inside Ω, away from dΩ
            # Obtain density at stencil points
            u_mid = max(U[D[i].xInd, D[i].yInd], 0.0)
            if ϕ[D[i].xInd-1, D[i].yInd] >= 0
                u_left = get_interface_density(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, uf)
            else
                u_left = U[D[i].xInd-1, D[i].yInd]
            end
            if ϕ[D[i].xInd+1, D[i].yInd] >= 0
                u_right = get_interface_density(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, uf)
            else
                u_right = U[D[i].xInd+1, D[i].yInd]
            end
            if ϕ[D[i].xInd, D[i].yInd-1] >= 0
                u_bottom = get_interface_density(D[i].xInd, D[i].yInd-1, D[i].xInd, D[i].yInd, dΩ, uf)
            else
                u_bottom = U[D[i].xInd, D[i].yInd-1]
            end
            if ϕ[D[i].xInd, D[i].yInd+1] >= 0
                u_top = get_interface_density(D[i].xInd, D[i].yInd, D[i].xInd, D[i].yInd+1, dΩ, uf)
            else
                u_top = U[D[i].xInd, D[i].yInd+1]
            end
            # Compute Laplacian and source term
            uxx = 2.0*( u_left/(D[i].θxm*(D[i].θxm + D[i].θxp))   - u_mid/(D[i].θxm*D[i].θxp) + u_right/(D[i].θxp*(D[i].θxm + D[i].θxp)) )/(dx^2)
            uyy = 2.0*( u_bottom/(D[i].θym*(D[i].θym + D[i].θyp)) - u_mid/(D[i].θym*D[i].θyp) + u_top/(D[i].θyp*(D[i].θym + D[i].θyp))   )/(dy^2)
            # End hack
            du[i] = par.D*(uxx + uyy)*u_mid^(par.m/(1 + par.m)) + par.λ*(1 + par.m)*u_mid*(1 - u_mid^(1/(1 + par.m)))
        elseif ((D[i].yInd == 1) || (D[i].yInd == par.Ny)) && (D[i].Ω == true) && (D[i].dΩ == false) # Boundary grid points inside Ω, away from dΩ
            # Obtain density at stencil points
            u_mid = max(U[D[i].xInd, D[i].yInd], 0.0)
            if ϕ[D[i].xInd-1, D[i].yInd] >= 0
                u_left = get_interface_density(D[i].xInd-1, D[i].yInd, D[i].xInd, D[i].yInd, dΩ, uf)
            else
                u_left = U[D[i].xInd-1, D[i].yInd]
            end
            if ϕ[D[i].xInd+1, D[i].yInd] >= 0
                u_right = get_interface_density(D[i].xInd, D[i].yInd, D[i].xInd+1, D[i].yInd, dΩ, uf)
            else
                u_right = U[D[i].xInd+1, D[i].yInd]
            end
            if ϕ[D[i].xInd, par.Ny-1] >= 0
                u_bottom = get_interface_density(D[i].xInd, par.Ny-1, D[i].xInd, par.Ny, dΩ, uf)
            else
                u_bottom = U[D[i].xInd, par.Ny-1]
            end
            if ϕ[D[i].xInd, 2] >= 0
                u_top =  get_interface_density(D[i].xInd, 1, D[i].xInd, 2, dΩ, uf)
            else
                u_top = U[D[i].xInd, 2]
            end
            # Compute Laplacian and source term
            uxx = 2.0*( u_left/(D[i].θxm*(D[i].θxm+D[i].θxp))   - u_mid/(D[i].θxm*D[i].θxp) + u_right/(D[i].θxp*(D[i].θxm+D[i].θxp)) )/(dx^2)
            uyy = 2.0*( u_bottom/(D[i].θym*(D[i].θym+D[i].θyp)) - u_mid/(D[i].θym*D[i].θyp) + u_top/(D[i].θyp*(D[i].θym+D[i].θyp))   )/(dy^2)
            du[i] = par.D*(uxx + uyy)*u_mid^(par.m/(1 + par.m)) + par.λ*(1 + par.m)*u_mid*(1 - u_mid^(1/(1 + par.m)))
        else # Grid points outside Ω or close to dΩ
            du[i] = 0.0
        end
    end
end