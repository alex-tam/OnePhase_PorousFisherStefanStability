# Solve BVP to generate initial condition
# Nizhum Rahman and Alex Tam, 8/12/2023

"Main function to generate initial condition"

"Data structure for model and numerical parameters"


function twicc(par)
    # 1. Parameters and domain
    par = Params(D=1.0) # Initialise structure for parameters
    zl = range(-par.L, 0.0, length = par.Nz)
    #zr = range(0.0, par.L, length = par.Nz)
    
   

   
    n=par.Nz-1
    

# Define the equation as a function
    function equation(r)
    n = par.Nz-1
       # Replace with your specific value of n
    target = -par.L / par.a  # The right-hand side of the equation
    #return (1+r)^(n-1) + target
    return (1-r^n)/(1-r)+target
    end
    dz2=zeros(n)
    xi=zeros(n+1)

# Use fzero to find the root (value of r)
    result = fzero(equation,1.01)  # You may adjust the initial interval [a, b] as needed

    #println("The value of r is approximately: ", result)
#dz[1]=a
    xi[n+1]=0
    
    for i=n:-1:1
    #term=term+1.0
    xi[i]=xi[i+1]-par.a*result^(n-i)
    dz2[i]=xi[i+1]-xi[i]
    #dz[i]=  -( (result)^(i-1))*a
    end
    xi[1]=-20
    #println(dz2)
    x=reverse(xi)
    #u=bvp_u(par, dz2, 0)
    #println(U)
    c,theta_0 = tws(par, dz2, 0.68, 1e-6)
    u0=(theta_0).^(1/(1+par.m))
   q=par.q#2*pi/5
    
ω, u2= lsa(par, c, u0, q, dz2, 1.0, 1e-6)
u1=u2 ./u0.^par.m 


 k=par.κ
 m=par.m

return vec(u0), vec(u1), vec(xi)
end

"Implement shooting method for travelling wave solutions"
function tws(par, dz2, c, ϵ)
    # Implement Newton's method
    for j = 1:20
        c_old = c # Store previous guess
        # Solve boundary-value problems
        u = bvp_u(par, dz2, c)
        #v = bvp_v(par, dz, c)
        up = bvp_u(par, dz2, c+ϵ)
        #vp = bvp_v(par, dz, c+ϵ)
        # Implement Newton's method
        f = c + 1/(1+par.m)*par.κ*dudz(u,dz2) # Interface condition
        fp = c + ϵ + 1/(1+par.m)*par.κ*dudz(up,dz2)  # Perturbed interface condition
        d = (fp-f)/ϵ # Approximate df/dϵ
        c = c - f/d # Newton iteration
        if norm(c - c_old) < 1e-6
            return c, bvp_u(par, dz2, c)
        end
    end
    @printf("Shooting method did not converge.\n")
end

"Finite difference approximation for du/dz at z = 0"
 function dudz(u, dz2)
          
         hp=dz2[end]
         hm=dz2[end-1]
        
    return hp/(hm*hp+hm^2)*(u[end-2]-u[end-1]*(1+2*hm/hp+hm^2/hp^2)-u[end]*(2*hm/hp+hm^2/hp^2))
    
 end

"Finite difference approximation for du/dz at z = 0"


"Nonliner boundary-value problem for u"
function bvp_u(par, dz2, c)
    u = ones(par.Nz) # Initial guess
    # Newton's method
    for i = 1:20000000
        u_old = u # Store previous iteration u
        F = Fu(par, dz2, u, c) # Construct F
        J = Ju(par, dz2, u, c) # Construct Jacobian
        u = u - J\F # Newton iteration
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    
    @printf("BVP for u did not converge.\n")
end

"Nonlinear boundary value problem for v"


"Vector function for u"
function Fu(par, dz2, u, c)
    
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = u[1] - 1.0 # Dirichlet condition on left boundary
    for i = 2:par.Nz-1        
            hp=dz2[i]
            hm=dz2[i-1]
            C=2/(hm*(hm+hp))
            B=-2/(hm*hp)
            A=2/(hp*(hm+hp))
       
        F[i] = (C*(u[i-1])+B*u[i]+A*u[i+1]) * (abs(u[i]))^(par.m/(1+par.m))+  c*(u[i+1]-u[i-1])/(hm+hp) + (par.m+1)* u[i]*(1-(abs(u[i]))^(1/(1+par.m)))
    end
    F[par.Nz] = u[par.Nz]- par.uf # Dirichlet condition at interface
    #println(F)
    return F
    
end

"Jacobian matrix for u"
function Ju(par, dz2, u, c)
    #c=0
    J = zeros(par.Nz, par.Nz) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nz, par.Nz] =1# 0.00001 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nz-1
        hp=dz2[i]
            hm=dz2[i-1]
            C=2/(hm*(hm+hp))
            B=-2/(hm*hp)
            A=2/(hp*(hm+hp))
        J[i,i-1] =(abs(u[i]))^(par.m/(1+par.m)) *C- c/(hm+hp) #- c/(2*dz)
        J[i,i] = B*((2*par.m+1)/(1+par.m))*(abs(u[i]))^(par.m/(1+par.m))  + (par.m+1) - (par.m+1)* ((2*par.m+1)/(1+par.m))*(abs(u[i]))^(par.m/(1+par.m))
        J[i,i+1] = A*(abs(u[i]))^(par.m/(1+par.m))  +c/(hm+hp)#+ c/(2*dz)
    end
    #println(J)
    return J
    
end



"Main function to compute travelling wave"


"Implement shooting method for linear stability analysis"
function lsa(par, c, u0, q, dz2, ω, ϵ)
    # Implement Newton's method
    for j = 1:300
        ω_old = ω # Store previous guess
        #
        u=u_direct(par, dz2, c, q, ω, u0)
        
        up = u_direct(par, dz2, c, q, ω+ϵ, u0)
        
        f = ω + par.κ*dudz(u,dz2)  # Interface condition
        fp = ω + ϵ +par.κ *dudz(up,dz2)  # Perturbed interface condition
        d = (fp-f)/ϵ # Approximate df/dϵ
        ω = ω - f/d # Newton iteration
        if abs(ω - ω_old) < 1e-6
            return ω, u_direct(par, dz2, c, q, ω, u0)
        end
    end
    @printf("LSA Shooting method did not converge.\n")
end

"Nonlinear boundary-value problem for u1"
function bvp_u1(par, dz2, c, q, ω, u0)
    u =ones(par.Nz) # Initial guess
   
    # Newton's method
    for i = 1:20000000
        u_old = u # Store previous iteration u
        F = Fu1(par, dz2, u, c, q, ω, u0) # Construct F
        J = Ju1(par, dz2, c, q, ω, u0) # Construct Jacobian
        u = u - J\F # Newton iteration
        if norm(u - u_old) < 1e-12
            return u
        end
    end
    @printf("TWS BVP for u did not converge.\n")
end





function u_direct(par, dz2, c, q, ω, u0)

    function fillmat!(M,par, dz2, c, q, ω, u0)

        M[1,1]=1.0
        M[par.Nz,par.Nz]=1.0
        
        for i=2:par.Nz-1
            hp=dz2[i]
            hm=dz2[i-1]
            C=2/(hm*(hm+hp))
            B=-2/(hm*hp)
            A=2/(hp*(hm+hp))
          M[i,i] = B*u0[i].^(par.m+1)-c*par.m*(u0[i+1]-u0[i-1])/(hm+hp)-ω*u0[i]-u0[i]*(((q^2*u0[i]^(par.m-1)+2)*u0[i]-1))               # Diagonal
          if i<par.Nz
            M[i,i+1] = A*u0[i].^(par.m+1)+c*u0[i]/(hm+hp)           # Sup-diagonal
          end
          if i>1
            M[i,i-1] = C*u0[i].^(par.m+1)-c*u0[i]/(hm+hp)           # Sub-diagonal
          end
        end
      end
      
    M = zeros(par.Nz,par.Nz)
  fillmat!(M,par, dz2, c, q, ω, u0)
    
  # Create Right-hand side  
  F = Array{Float64}(undef,par.Nz)  
  for i=2:par.Nz-1
    hp=dz2[i]
    hm=dz2[i-1]
    F[i] = -u0[i].^(par.m+1)* (u0[i+1]-u0[i-1])/(hm+hp)* (ω+q^2*u0[i]^par.m)
  end
  F[1] =  0
  F[par.Nz] =  0

  u = M\F;
  return u


end



