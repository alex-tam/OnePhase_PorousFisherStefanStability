# Solve BVP to generate initial condition
# Nizhum Rahman and Alex Tam, 8/12/2023

"Main function to generate initial condition"
function twic(par)
    # Obtain grid using geometric progression
    n = par.Nξ-1
    dz2 = zeros(n)
    xi = zeros(n+1)
    result = fzero((r) -> equation(par, r), 1.01)
    xi[n+1] = 0
    for i=n:-1:1
        xi[i] = xi[i+1] - par.a*result^(n-i)
        dz2[i] = xi[i+1] - xi[i]
    end
    xi[1] = -20
    # Solve leading-order problem
    c, θ0 = tws(par, dz2, 0.68, 1e-6)
    u0 = θ0.^(1/(1 + par.m))
    return vec(u0), vec(xi)
end

"Implement shooting method for travelling wave solutions"
function tws(par, dz2, c, ϵ)
    # Implement Newton's method
    for j = 1:10000
        c_old = c # Store previous guess
        # Solve boundary-value problems
        u = bvp_u(par, dz2, c)
        up = bvp_u(par, dz2, c+ϵ)
        # Implement Newton's method
        f = c + 1/(1+par.m)*par.κ*dudz(u,dz2) # Interface condition
        fp = c + ϵ + 1/(1+par.m)*par.κ*dudz(up,dz2)  # Perturbed interface condition
        d = (fp-f)/ϵ # Approximate df/dϵ
        c = c - f/d # Newton iteration
        if norm(c - c_old) < 1e-6
            return c, bvp_u(par, dz2, c)
        end
    end
    @printf("Shooting method for TWS did not converge.\n")
end

"Finite difference approximation for du/dz at z = 0"
function dudz(u, dz2)
    hp = dz2[end]
    hm = dz2[end-1]    
    return hp/(hm*hp + hm^2)*(u[end-2] - u[end-1]*(1 + 2*hm/hp + hm^2/hp^2) - u[end]*(2*hm/hp + hm^2/hp^2))
end

"Nonliner boundary-value problem for u"
function bvp_u(par, dz2, c)
    u = ones(par.Nξ) # Initial guess
    # Newton's method
    for i = 1:10000
        u_old = u # Store previous iteration u
        F = Fu(par, dz2, u, c) # Construct F
        J = Ju(par, dz2, u, c) # Construct Jacobian
        u = u - J\F # Newton iteration
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    @printf("BVP for u0 did not converge.\n")
end

"Vector function for u"
function Fu(par, dz2, u, c)
    F = Vector{Float64}(undef, par.Nξ) # Pre-allocate F
    F[1] = u[1] - 1.0 # Dirichlet condition on left boundary
    for i = 2:par.Nξ-1        
        hp = dz2[i]
        hm = dz2[i-1]
        C = 2/(hm*(hm+hp))
        B = -2/(hm*hp)
        A = 2/(hp*(hm+hp))
        F[i] = (C*(u[i-1]) + B*u[i] + A*u[i+1])*u[i]^(par.m/(1 + par.m)) + c*(u[i+1] - u[i-1])/(hm + hp) + (par.m + 1)*u[i]*(1 - u[i]^(1/(1+par.m)))
    end
    F[par.Nξ] = u[par.Nξ]- par.uf # Dirichlet condition at interface
    return F
end

"Jacobian matrix for u"
function Ju(par, dz2, u, c)
    J = zeros(par.Nξ, par.Nξ) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nξ, par.Nξ] = 1.0 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nξ-1
        hp = dz2[i]
        hm = dz2[i-1]
        C = 2/(hm*(hm+hp))
        B = -2/(hm*hp)
        A = 2/(hp*(hm+hp))
        J[i,i-1] = u[i]^(par.m/(1+par.m))*C - c/(hm+hp)
        J[i,i] = B*((2*par.m+1)/(1+par.m))*u[i]^(par.m/(1+par.m))  + (par.m+1) - (par.m+1)*((2*par.m+1)/(1+par.m))*u[i]^(par.m/(1+par.m))
        J[i,i+1] = A*u[i]^(par.m/(1+par.m)) + c/(hm+hp)
    end
    return J
end

"Function for finding grid for geometric progression"
function equation(par, r)
    n = par.Nξ-1
    target = -par.Lξ/par.a  # The right-hand side of the equation
    return (1-r^n)/(1-r) + target
end