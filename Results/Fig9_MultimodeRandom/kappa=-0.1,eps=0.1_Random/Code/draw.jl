# Plot solutions to 2D Porous-Fisher-Stefan
# Alex Tam, 19/05/2024

# Load packages
using Plots
using Measures
using LaTeXStrings
using DelimitedFiles
using FFTW

#########################################
##### Solution processing functions #####
#########################################
"Construct y-dependent interface from level-set function"
function interface(ϕ, x, y, ny, dx)
    # Pre-allocate
    L = Vector{Float64}(undef, ny)
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
    return L
end

##############################
##### Plotting functions #####
##############################
"Plot initial condition as 2D heat maps"
function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    # Draw heat map (u)
    p1 = heatmap(x, y, transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    contour!(p1, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.000001], clim=(0.0, 1.0))
    savefig("u-$i.pdf")
    # Draw heat map (ϕ)
    p2 = heatmap(x, y, transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    contour!(p2, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.0])
    savefig("phi-$i.pdf")
end

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    # Draw heat map (u)
    p1 = heatmap(x, y, transpose(U), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma, clims=(0.0, 1.0))
    contour!(p1, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.000001], clim=(0.0, 1.0))
    savefig("u-$i.pdf")
    # Draw heat map (v)
    heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    savefig("V-$i.pdf")
    # Draw heat map (ϕ)
    p2 = heatmap(x, y, transpose(ϕ), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    contour!(p2, x, y, transpose(ϕ), c=:red, linewidth = 2, levels=[0.0])
    savefig("phi-$i.pdf")
end

"Density slice plots"
function draw_slices(x, y, nx, ny, Lx, Ly, plot_times, m)
    # Import initial condition
    Φ = readdlm("U-0.csv")
    U = Φ.^(1/(1+m))
    # Import and plot x-slice
    plot(x, U[:,ny], xlabel = L"$x$", ylabel = L"$u(x,5,t)$", linecolor = :black, linewidth = 2, aspect_ratio = 20.0, grid = false, margin=3mm, legend = false, xlims=(0,Lx), ylims=(0,1))
    for i in plot_times
        Φx = (readdlm("ux-$i.csv"))
        ux = Φx.^(1/(1+m))
        plot!(x, ux, linecolor = :black, linestyle = :solid, linewidth = 1, legend=false)
    end
    savefig("u_slice_x.pdf")
    # Import and plot y-slice
    plot(y, U[nx,:], xlabel = L"$y$", ylabel = L"$u(5,y,t)$", linecolor = :black, linewidth = 2, aspect_ratio = 20.0, grid = false, margin=3mm, legend = false, xlims=(0,Ly), ylims=(0,1))
    for i in plot_times
        Φy = (readdlm("uy-$i.csv"))
        uy = Φy.^(1/(1+m))
        plot!(y, uy, linecolor = :black, linestyle = :dash, linewidth = 2)
    end
    savefig("u_slice_y.pdf")
end

"Interface and power spectrum plots"
function draw_interface(ϕ, x, y, dx, Lx, Ly, mean_front_pos, i)
    # Obtain interface
    L = interface(ϕ, x, y, length(y), dx)
    # Plot interface
    plot(y, L, xlabel = L"$y$", ylabel = L"$L(y, 0)$", linecolor = :black, linewidth = 2, margin=3mm, legend = false, xlims=(0, Ly), ylims=(0, Lx))
    savefig("interface-$i.pdf")
    writedlm("interface-$i.csv", L)
    # Compute discrete cosine transform
    d = dct(L .- mean_front_pos[i+1], 1) # Compute discrete cosine transform of front position (centred at zero)
    # Plot spectrum of front position
    plot(0:length(d)-1, abs.(d).*sqrt(2/length(y)), xlabel = L"$k$", ylabel = L"$|f_k|$", linecolor = :black, linewidth = 2, margin=3mm, legend = false, xlims=(0, 50), ylims = (0,0.0301))
    savefig("dct-$i.pdf")
    writedlm("k.csv", 0:length(d)-1)
    writedlm("dct-$i.csv", d)
end

"Perturbation amplitude plot"
function draw_amplitude(t, A)
    plot(t, A, xlabel = L"$t$", ylabel = L"$A(t)$", linecolor = :black, linewidth = 2, grid = false, margin=5mm, legend = false, xlims=(0,10), ylims=(0,0.101))
    savefig("Amplitude.pdf")
end

#########################
##### Main function #####
#########################
"Main function for plotting"
function draw()
    ##### Configure plots #####
    gr()
    plot()
    default(fontfamily = "Computer Modern", titlefontsize = 24, guidefontsize = 24, tickfontsize = 24, legendfontsize = 20)
    ##### Import data #####
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv"))
    dx = x[2] - x[1] 
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)+1)/2
    ny::Int = (length(y)+1)/2
    Lx = maximum(x)
    Ly = maximum(y)
    m::Float64 = readdlm("m.csv")[1, 1]
    Φ = (readdlm("U-0.csv"))
    U = Φ.^(1/(m+1))
    ϕ = readdlm("Phi-0.csv")
    mean_front_pos = readdlm("L.csv") # Mean front position
    ##### Plot initial condition and interface ####
    draw_heat(x, y, U, ϕ, 0, Lx, Ly) # Heat maps
    draw_interface(ϕ, x, y, dx, Lx, Ly, mean_front_pos, 0) # Interface plot and spectrum
    ##### Plot solution and interface ####
    for i in plot_times
        Φ = (readdlm("U-$i.csv"))
        U = Φ.^(1/(m+1))
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
        draw_interface(ϕ, x, y, dx, Lx, Ly, mean_front_pos, i) # Interface plot and spectrum
    end
    ##### Slice plots #####
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times, m) # Slice plots
    ##### Amplitude plot #####
    t = vec(readdlm("t.csv"))
    A = vec(readdlm("Amp.csv"))
    draw_amplitude(t, A)
end

@time draw()