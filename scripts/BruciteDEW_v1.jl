using LinearAlgebra, CairoMakie, Printf, ForwardDiff, CSV, DataFrames, MathTeXEngine
Makie.inline!(true)
Makie.update_theme!(fonts=(regular=texfont(), bold=texfont(:bold), italic=texfont(:italic)))

# Comparison with fluid speciation from EQ3 using the DEW thermodynamic database
# Jacobian constructed with automatic differentiation

"""Compute activity coefficient for fluid species """
function ActivityCoeff(z, Aγ, Bγ, coeff, å, I, T)
    logγ = ones(length(z))
    for i in eachindex(z)
        if (z[i] == 0)        # Neutral fluid species
            logγ[i] = ((coeff[1] + coeff[2]*T + coeff[3]/T)*I - (coeff[4] + coeff[5]*T)*(I/(I+1)))/2.303
        else                  # Charged fluid species
            logγ[i] = -((Aγ * z[i]^2 * I^0.5)/(1 + å[i] * Bγ * I^0.5))
        end
    end
    return logγ
end

""" Evaluation of residual """
function Residual!(f, m, D, b, Ḃ, T)

    # Compute activity coefficients
    I            = (1/2) * m'*(Ḃ.z.^2)
    logγ_charged = -((Ḃ.Aγ * Ḃ.z.^2 * I) ./ (1 .+ Ḃ.å * Ḃ.Bγ * I))
    logγ_neutral = ((Ḃ.coeff[1] + Ḃ.coeff[2]*T + Ḃ.coeff[3]/T)*I - (Ḃ.coeff[4] + Ḃ.coeff[5]*T)*(I/(I+1)))/2.303
    logγ         = logγ_charged.*(Ḃ.z .!= 0) .+ logγ_neutral.*(Ḃ.z .== 0)

    for i in eachindex(f)
        if i < 3
            f[i] = D[i, :]' * m - b[i]
        else
            f[i] = D[i, :]' * log10.(m) + D[i, :]' * logγ - b[i]
        end
    end
end

""" Optimize globalization parameter of Newton-Raphson solver """
function LineSearch(f, m, Δm, D, b, Ḃ, T)

    α  = [0.1 0.25 0.5 0.8 1.0]
    F  = zero(α)
    m1 = copy(m)

    for i in eachindex(α)
        m1 .= m + α[i] .* Δm
        if minimum(m1) < 0.0 # Do not compute if m<0 otherwise 'Domain Error'
            f .= 1000.0      # In that case set residual to artificially high value (makes sure that it's never identifed as a minimum)
        else
            Residual!(f, m1, D, b, Ḃ, T)
        end
        F[i] = norm(f)
    end
    
    _, ind = findmin(F)
    return α[ind]
end

# Solve for log(activities)

""" Main function for speciation """
function Speciation()

    # Read data into a dataframe
    df = (CSV.read("data/MatrixBruciteDEW.csv", DataFrame))
    species      = names(df)[2:end] # Read column labels

    # Parameters
    Clᵗᵒᵗ = 0.01     # Total chlorinity
    T     = 673      # Temperature in Kelvin

    # Newton-Raphson solver
    niter = 100     # Max number of iterations
    ϵ     = 1e-12   # Non-linear tolerance
    iter  = 0       # Iteration count

    # Thermodynamic parameters for activity coefficients
    Ḃ = (
        Aγ    = 0.7879,
        Bγ    = 0.3711,
        coeff = [-1.0312, 0.0012806, 255.9, 0.445, -0.001606],   # Coefficients for neutral species C, F, G, E and H
        å     = [3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7, 3.7],        # Size of fluid species (including hydration shell)
        z     = collect(df[1,2:end]),                 # Charges for the fluid species
    )

    # Arrays
    D    = Float64.(Matrix(df[:, 2:end]))       # Read coefficients and convert to Matrix
    b    = [0.0; Clᵗᵒᵗ; 6.8466; -1.0841; -0.6078; -8.1764; 6.6296; 4.9398]  
    n    = length(species)
    m    = 0.01 * ones(length(b))               # Initial condition
    f    = zero(m)
    Δm   = zero(m)
    J    = zeros(n, n)
    e    = zeros(niter)
    
    for _ = 1:niter

        iter += 1

        # Residual evaluation
        Residual!(f, m, D, b, Ḃ, T)   
        e[iter] = norm(f)
        if norm(f) < ϵ
            break
        end

        # Automatic Jacobian generation
        r = (f, m) -> Residual!(f, m, D, b, Ḃ, T) # closure
        J = ForwardDiff.jacobian(r, f, m)

        # Update molalities
        Δm .= -J \ f
        α   = LineSearch(f, m, Δm, D, b, Ḃ, T) # Optminise with line search
        m .+= α * Δm
        @printf("it. %03d --- f = %1.4e --- α = %1.2f\n", iter, norm(f), α)
    end

    for i in 1:n
        @printf("Molality of %s is %1.4e\n", species[i], m[i])
    end

    # Figure 
    f = Figure(size=(1200, 600), fontsize=25, aspect=2.0)
    ax1 = Axis(f[1, 1], title=L"$$Molalities", xlabel=L"$$Species", ylabel=L"$$Molality [-]", xgridvisible=false, ygridvisible=false)
    ax1.xticks = (collect(1:length(species)), species)
    scatter!(ax1, 1:length(m), m, label="a")
    ax2 = Axis(f[1, 2], title=L"$$Convergence", xlabel=L"$$Iterations", ylabel=L"$$Error", xgridvisible=false, ygridvisible=false)
    lines!(ax2, 1:iter, log10.(e[1:iter]))
    DataInspector(f)
    display(f)

end

Speciation()