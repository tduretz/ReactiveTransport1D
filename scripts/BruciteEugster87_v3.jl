using LinearAlgebra, CairoMakie, Printf, ForwardDiff, CSV, DataFrames, MathTeXEngine
Makie.inline!(true)
Makie.update_theme!(fonts=(regular=texfont(), bold=texfont(:bold), italic=texfont(:italic)))

# Jacobian constructed with automatic differentiation

""" Evaluation of residual """
function Residual!(f, m, D, b)
    for i in eachindex(f)
        if i < 3
            f[i] = D[i, :]' * m - b[i]
        else
            f[i] = D[i, :]' * log10.(m) - b[i]
        end
    end

end

""" Optimize globalization parameter of Newton-Raphson solver """
function LineSearch(f, m, Δm, D, b)

    α = [0.1 0.25 0.5 0.8 1.0]
    F = zero(α)
    m1 = copy(m)

    for i in eachindex(α)
        m1 .= m + α[i] .* Δm
        if minimum(m1) < 0.0
            f .= 1000.0
        else
            Residual!(f, m1, D, b)
        end
        F[i] = norm(f)
    end
    _, ind = findmin(F)

    return α[ind]
end

# Solve for log(activities)
species = ["MgCl₂", "HCl", "MgCl⁺", "Mg²⁺", "Cl⁻", "H⁺", "OH⁻"]

""" Main function for speciation """
function Speciation()

    # Read data into a dataframe
    df_Eugster87 = (CSV.read("data/MatrixBruciteEugster87.csv", DataFrame))

    niter = 100
    iter = 0

    # Convert to Matrix
    D = Float64.(Matrix(df_Eugster87[:, 2:end]))
    b = [0.0; 1.0; 5.04; -1.68; -3.62; -2.79; -10.2]
    n = length(species)
    m = 0.01 * ones(length(b)) # Initial condition
    f = zero(m)
    J = zeros(n, n)
    e = zeros(niter)

    for _ = 1:niter

        iter += 1

        # Residual evaluation
        Residual!(f, m, D, b)
        e[iter] = norm(f)
        if norm(f) < 1e-12
            break
        end

        # Automatic Jacobian generation
        r = (f, m) -> Residual!(f, m, D, b) # closure
        J = ForwardDiff.jacobian(r, f, m)

        # Update molalities
        Δm = -J \ f
        α = LineSearch(f, m, Δm, D, b) # Optminise with line search
        m .+= α * Δm
        @printf("it. %03d --- f = %1.4e --- α = %1.2f\n", iter, norm(f), α)
    end

    # Figure 
    f = Figure(resolution=(1200, 600), fontsize=25, aspect=2.0)
    ax1 = Axis(f[1, 1], title=L"$$Molalities", xlabel=L"$$Species", ylabel=L"$$Molality [-]", xgridvisible=false, ygridvisible=false)
    ax1.xticks = (collect(1:length(species)), species)
    scatter!(ax1, 1:length(m), m, label="a")
    ax2 = Axis(f[1, 2], title=L"$$Convergence", xlabel=L"$$Iterations", ylabel=L"$$Error", xgridvisible=false, ygridvisible=false)
    lines!(ax2, 1:iter, log10.(e[1:iter]))
    DataInspector(f)
    display(f)

end

Speciation()

