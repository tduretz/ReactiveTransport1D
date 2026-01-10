using LinearAlgebra, CairoMakie, Printf, ForwardDiff, CSV, DataFrames, MathTeXEngine
using MAGEMin_C
Makie.inline!(true)
Makie.update_theme!(fonts=(regular=texfont(), bold=texfont(:bold), italic=texfont(:italic)))

# Comparison with fluid speciation from EQ3 using the DEW thermodynamic database
# Jacobian constructed with automatic differentiation

"""Compute activity coefficient for fluid species """
function ActivityCoeff(z, bdot, coeff, å, I, T)
    logγ = ones(length(z))
    for i in eachindex(z)
        if (z[i] == 0)        # Neutral fluid species
            logγ[i] = ((coeff[1] + coeff[2] * T + coeff[3] / T) * I - (coeff[4] + coeff[5] * T) * (I / (I + 1))) / 2.303
            # logγ[i] = 0
        else                  # Charged fluid species
            logγ[i] = -((bdot[1] * z[i]^2 * I^0.5) / (1 + å[i] * bdot[2] * I^0.5))
        end
    end

    return logγ
end

"""Compute water activity """
function ActivityWater(m, bdot, å, I, Ω, T)
    logaw = (-(sum(m) / 2.303) + ((2 * bdot[1] * I^(3 / 2) * Sigma(å * bdot[2] * I^(0.5))) / 3)) / Ω

    return logaw
end

"""Compute sigma equation for water activity """
function Sigma(x)
    σ = (3 / (x^3)) * (1 + x - (1 / (1 + x)) - 2 * log(1 + x))

    return σ
end

""" Evaluation of residual """
function Residual!(f, m, logγ, D, cH2O, logaw, b)
    for i in eachindex(f)
        if i < 3
            f[i] = D[i, :]' * m - b[i]
        else
            f[i] = D[i, :]' * log10.(m) + D[i, :]' * logγ + cH2O[i] * logaw - b[i]
        end
    end
end

""" Optimize globalization parameter of Newton-Raphson solver """
function LineSearch(f, m, logγ, Δm, D, cH2O, logaw, b)

    α = [0.1 0.25 0.5 0.8 1.0]
    F = zero(α)
    m1 = copy(m)

    for i in eachindex(α)
        m1 .= m + α[i] .* Δm
        if minimum(m1) < 0.0 # Do not compute if m<0 otherwise 'Domain Error'
            f .= 1000.0      # In that case set residual to artificially high value (makes sure that it's never identifed as a minimum)
        else
            Residual!(f, m1, logγ, D, cH2O, logaw, b)
        end
        F[i] = norm(f)
    end

    _, ind = findmin(F)
    return α[ind]
end

# Solve for log(activities)

""" Main function for speciation """
function Speciation(logaoxides, T_calc)

    # Read data into a dataframe
    df = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MatrixBruciteAntigorite_400C_BackCalc_SiO2(aq)_wo_Mg(OH)2_Fe.csv", DataFrame))
    species = names(df)[2:end-2] # Read column labels

    # Thermodynamic parameters for activity coefficients
    Aᵧ = 0.7879
    Bᵧ = 0.3711
    bdot = [Aᵧ, Bᵧ]                                         # B-dot coefficients Aᵧ and Bᵧ at 400 C and 5 kbar
    coeff = [-1.0312, 0.0012806, 255.9, 0.445, -0.001606]   # Coefficients for neutral species C, F, G, E and H

    # Parameters
    Clᵗᵒᵗ = 0.1                                    # Total chlorinity
    T = T_calc + 273                               # Temperature in Kelvin
    å_w = 4                                        # Hard core diameter for water 4.0 Å
    Ω = 55.55                                      # Water constant (1000 divide by the moelcular weight of water)

    # Newton-Raphson solver
    niter = 100                                    # Max number of iterations
    ϵ = 1e-10                                      # Non-linear tolerance
    iter = 0                                       # Iteration count

    # Arrays
    D = Float64.(Matrix(df[:, 2:end-2]))           # Read coefficients and convert to Matrix
    z = collect(df[1, 2:end-2])                    # Charges for the fluid species
    b = collect(df[1:end, end])                    # Get the Cltot and log K for each reaction/law of mass action
    b[3] = logaoxides[1]
    b[4] = logaoxides[2]
    b[5] = logaoxides[3]
    coeffH2O = collect(df[1:end, end-1])
    # print(b[8])
    å = 3.7 * ones(length(b))                      # Size of fluid species (including hydration shell)
    # b    = [0.0; Clᵗᵒᵗ; 6.8466; -1.0841; -0.6078; -8.1764; 6.6296; 4.9398]  
    n = length(species)
    m = 0.01 * ones(length(b))                     # Initial condition
    m[6] = 10^(b[8]/2)
    m[7] = 10^(b[8]/2)
    # m = [0.02; 0.01; 0.01; 0.1; 0.0005; 0.0008; 0.3; 0.01; 0.01]
    logγ = ones(length(b))                         # Initial activity coefficients equal to 1
    f = zero(m)
    Δm = zero(m)
    J = zeros(n, n)
    e = zeros(niter)
    logaw = 0

    for _ = 1:niter

        iter += 1

        # Compute activity coefficients
        I = (1 / 2) * sum(m .* z .^ 2)
        logγ = ActivityCoeff(z, bdot, coeff, å, I, T)
        @printf("Ionic strength = %1.6e\n", I)

        # Residual evaluation
        Residual!(f, m, logγ, D, coeffH2O, logaw, b)
        e[iter] = norm(f)
        if norm(f) < ϵ
            break
        end

        # Automatic Jacobian generation
        r = (f, m) -> Residual!(f, m, logγ, D, coeffH2O, logaw, b) # closure
        J = ForwardDiff.jacobian(r, f, m)

        # Update molalities
        Δm .= -J \ f
        α = LineSearch(f, m, logγ, Δm, D, coeffH2O, logaw, b) # Optminise with line search
        m .+= α * Δm
        @printf("it. %03d --- f = %1.4e --- α = %1.2f\n", iter, norm(f), α)
        for i in 1:n
            @printf("Molality of %s is %1.4e\n", species[i], m[i])
        end

        # Compute water activity
        logaw = ActivityWater(m, bdot, å_w, I, Ω, T)
        aw = 10^logaw
        # aw = 1
        # logaw = 0
        @printf("Water activity = %1.4e\n", aw)
    end

    for i in 1:n
        @printf("Molality of %s is %1.4e\n", species[i], m[i])
    end
    @printf("pH is %1.4e\n", -log10(m[6]))
    @printf("Mg dissolved is %1.4e\n", m[2]+m[3]+m[7])

    # Figure 
    # f = Figure(size=(1200, 600), fontsize=25, aspect=2.0)
    # ax1 = Axis(f[1, 1], title=L"$$Molalities", xlabel=L"$$Species", ylabel=L"$$Molality [-]", xgridvisible=false, ygridvisible=false)
    # ax1.xticks = (collect(1:length(species)), species)
    # scatter!(ax1, 1:length(m), m, label="a")
    # ax2 = Axis(f[1, 2], title=L"$$Convergence", xlabel=L"$$Iterations", ylabel=L"$$Error", xgridvisible=false, ygridvisible=false)
    # lines!(ax2, 1:iter, log10.(e[1:iter]))
    # DataInspector(f)
    # display(f)

end

function GetChemicalPotentials(X, Xoxides, data, T_calc, P, sys_in)
    out = single_point_minimization(P, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    for n in 1:length(out.sol_name)
        @printf("%s is stable with %2.10e vol\n", out.sol_name[n],out.ph_frac_vol[n])
    end
    for m in 1:length(Xoxides)
        @printf("The chemical potential of %s = %2.10e (kJ/mol)\n", Xoxides[m],out.Gamma[m])
    end
    # print(out.SS_vec[6].siteFractionsNames)
    # print(out.SS_vec[6].siteFractions)

    return out.Gamma[1], out.Gamma[2], out.Gamma[3], out.Gamma[4]

end

# Initialization for the MAGEMin conditions
P = 5.0
T_calc = 400.0
data = Initialize_MAGEMin("ume", verbose=false);
Xoxides = ["SiO2"; "FeO"; "MgO"; "H2O"; "Al2O3"; "O"; "S"];  # System of reduced serpentinite of Evans & frost (2021)
X_comp = [34.146613; 6.415533; 33.41302; 23.883372; 1.808672; 0.060068; 0.272721];   # Composition of reduced serpentinite of Evans & frost (2021) in wt.%
sys_in = "wt"
# Xoxides = ["SiO2"; "FeO"; "MgO"; "H2O"];  # System of component for simple Olivine + H2O
# X_comp = [22.74; 4.15; 38.33; 34.78];   # Composition of simple Olivine + H2O in mol%
# sys_in = "mol"

# Get the chemical potentials values for each component
µ_SiO₂, µ_FeO, µ_MgO, µ_H₂O = GetChemicalPotentials(X_comp, Xoxides, data, T_calc, P, sys_in)
@printf("Chemical potentials (kJ/mol) of SiO2 = %2.10e, FeO = %2.10e, MgO = %2.10e and of H2O = %2.10e\n", µ_SiO₂, µ_FeO, µ_MgO, µ_H₂O)

# Compute the log(aM/aH+)
S0_SiO₂ = 223.96 
S0_MgO = 135.255
S0_FeO = 129.855
S0_H₂O = 233.255
# µ_SiO₂ = -968.919674    # SiO2 chemical potential from PerpleX with HSC convention
# µ_MgO = -636.3732783    # MgO chemical potential from PerpleX with HSC convention
# µ_H₂O = -335.0539783    # H2O chemical potential from PerpleX with HSC convention
# µ_FeO = -310.579        # FeO chemical potential from PerpleX with HSC convention
µ_FeO = -310.9142683    # FeO chemical potential from PerpleX with HSC convention at 400 °C
# µ_FeO = -331.615        # FeO chemical potential from PerpleX with HSC convention at 500 °C
R = 8.314               # Gas constant (J/mol/K) 
G_Mg⁰ = -417093.5       # Gibbs free energy (J/mol) of Mg2+ at P and T from PerpleX
G_Fe⁰ = -63615.6        # Gibbs free energy (J/mol) of Fe2+ at P and T from PerpleX
G_SiO₂⁰ = -851437.5     # Gibbs free energy (J/mol) of SiO2(aq) at P and T from PerpleX
logMgH = (1000*µ_MgO + S0_MgO*298.15 - 1000*µ_H₂O - G_Mg⁰) / (2.303 * R * (T_calc+273.15))
logFeH = (1000*µ_FeO + S0_FeO*298.15 - 1000*µ_H₂O - G_Fe⁰) / (2.303 * R * (T_calc+273.15))
logSiO₂H = (1000*µ_SiO₂ + S0_SiO₂*298.15 - G_SiO₂⁰) / (2.303 * R * (T_calc+273.15))
@printf("Log(aSiO2) = %2.10e, log(aFe2+/aH+) = %2.10e and log(aMg2+/aH+) = %2.10e\n", logSiO₂H, logFeH, logMgH)
logaoxides = [logMgH; logSiO₂H; logFeH]

# Compute the speciation using the log(aMg2+/aH+) from MAGEMin
Speciation(logaoxides, T_calc)

