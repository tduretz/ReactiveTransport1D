using LinearAlgebra, CairoMakie, Printf, ForwardDiff, CSV, DataFrames, MathTeXEngine
using MAGEMin_C
Makie.inline!(true)
Makie.update_theme!(fonts=(regular=texfont(), bold=texfont(:bold), italic=texfont(:italic)))

# Comparison with fluid speciation from EQ3 using the DEW thermodynamic database
# Jacobian constructed with automatic differentiation

function Itp1D_scalar1(xlt, varlt, xdata, dx, xmin)
    """Interpolation function through one variable:
    Inputs:
    - xlt: x points [1D array]
    - varlt: y variable at the x points [1D array]
    - xdata: x value at which we want to know the value of y [scalar]
    - dx: spacing between two values [scalar]
    - xmin: minimum value of x [scalar]

    Return the y value at the point of interest x [scalar]
    """
    iW = Int(floor((xdata - xmin) / dx) + 1)
    wW = 1.0 - (xdata - xlt[iW]) / dx
    return wW * varlt[iW] + (1.0 - wW) * varlt[iW+1]
end

@views function Itp1D_rev_scalar1(xlt, varlt, xdata)
    xinf_id = sum(xlt .- xdata .< 0)
    if xinf_id < 1
        xinf_id = 1
    end
    if xinf_id > length(xlt) - 1
        xinf_id = length(xlt) - 1
    end
    xinf_dist = (xlt[xinf_id+1] - xdata) / (xlt[xinf_id+1] - xlt[xinf_id])
    return xinf_dist * varlt[xinf_id] + (1.0 - xinf_dist) * varlt[xinf_id+1]
end

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
function Speciation(logaoxides, T_calc, P)

    # Read data into a dataframe
    df = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MatrixBruciteAntigorite_400C_BackCalc_SiO2(aq)_wo_Mg(OH)2.csv", DataFrame))
    species = names(df)[2:end-1] # Read column labels

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
    ϵ = 1e-12                                      # Non-linear tolerance
    iter = 0                                       # Iteration count

    # Log Ks lookup tables at variable pressure and 400 °C
    MgCl = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MgCl+_dissociation_500C.csv", DataFrame))
    P_logK    = collect(MgCl[:,1])                 # Pressure vector for lookup tables of log K values
    LogK_MgCl = collect(MgCl[:,2])                 # Log K values for MgCl+ dissociation
    HCl = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/HCl_dissociation_500C.csv", DataFrame))
    LogK_HCl = collect(HCl[:,2])                   # Log K values for HCl dissociation
    H2O = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/H2O_dissociation_500C.csv", DataFrame))
    LogK_H2O = collect(H2O[:,2])                   # Log K values for H2O dissociation
    MgOH = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MgOH+_dissociation_500C.csv", DataFrame))
    LogK_MgOH = collect(MgOH[:,2])                 # Log K values for MgOH dissociation


    # Arrays
    D = Float64.(Matrix(df[:, 2:end-1]))           # Read coefficients and convert to Matrix
    z = collect(df[1, 2:end-1])                    # Charges for the fluid species
    # b = collect(df[1:end, end])                    # Get the Cltot and log K for each reaction/law of mass action (previous version)
    coeffH2O = collect(df[1:end, end])
    b = 0.01 * ones(length(coeffH2O))
    b[1] = 0.0
    b[2] = Clᵗᵒᵗ
    b[3] = logaoxides[1]
    b[4] = logaoxides[2]
    b[5] = Itp1D_rev_scalar1(P_logK, LogK_MgCl, 1000*P)  # Interpolate log Ks for MgCl+ dissociation constant
    b[6] = Itp1D_rev_scalar1(P_logK, LogK_HCl, 1000*P)   # Interpolate log Ks for HCl dissociation constant
    b[7] = Itp1D_rev_scalar1(P_logK, LogK_H2O, 1000*P)   # Interpolate log Ks for H2O dissociation constant
    b[8] = Itp1D_rev_scalar1(P_logK, LogK_MgOH, 1000*P)  # Interpolate log Ks for MgOH+ dissociation constant
    å = 3.7 * ones(length(b))                      # Size of fluid species (including hydration shell)
    # b    = [0.0; Clᵗᵒᵗ; 6.8466; -1.0841; -0.6078; -8.1764; 6.6296; 4.9398]  
    n = length(species)
    m = 0.01 * ones(length(b))                     # Initial condition
    m[5] = 10^(b[7]/2)
    m[6] = 10^(b[7]/2)
    print(m)
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
        @printf("pH is %1.4e\n", -log10(m[5]))
        @printf("Dissolved Mg is %1.4e\n", m[2]+m[3]+m[7])
        @printf("Dissolved Si is %1.4e\n", m[8])
    end

    for i in 1:n
        @printf("Molality of %s is %1.4e\n", species[i], m[i])
    end

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
        @printf("%s is stable with %2.2e vol\n", out.sol_name[n],out.ph_frac_vol[n])
    end
    for m in 1:length(Xoxides)
        @printf("The chemical potential of %s = %2.5e (kJ/mol)\n", Xoxides[m],out.Gamma[m])
    end

    return out.Gamma[1], out.Gamma[3], out.Gamma[5]

end

# Initialization for the MAGEMin conditions
P = 5.0
T_calc = 500.0
data = Initialize_MAGEMin("ume", verbose=false);
Xoxides = ["SiO2"; "FeO"; "MgO"; "Al2O3"; "H2O"; "O"];  # System of component for reduced serpentinites from Evans & Frost (2021)
X_comp = [34.146613; 6.415533; 33.41302; 1.808672; 23.883372; 0.060068];   # Composition of reduced serpentinites from Evans & Frost (2021)
sys_in = "wt"
# Xoxides = ["SiO2"; "FeO"; "MgO"; "H2O"];  # System of component for simple Ol (Fo90) + H2O
# X_comp = [22.74; 4.15; 38.33; 34.78];     # Composition of simple olivine (Fo90) + H2O to saturate at Br + Atg
# sys_in = "mol"

# Get the chemical potentials values for each component
µ_SiO₂, µ_MgO, µ_H₂O = GetChemicalPotentials(X_comp, Xoxides, data, T_calc, P, sys_in)
@printf("Chemical potentials (kJ) of SiO2 = %2.5e, MgO = %2.5e and of H2O = %2.5e\n", µ_SiO₂, µ_MgO, µ_H₂O)

# Thermodynamic properties
# Entropies needed to convert HSC convention Gibbs free energies/chemical potentials to SUPCRT convention
S0_SiO₂    = 223.96              # Entropy of SiO2 at 298 K and 1 bar (J/mol/K)
S0_MgO     = 135.255             # Entropy of MgO at 298 K and 1 bar (J/mol/K)
S0_H₂O     = 233.255             # Entropy of H2O at 298 K and 1 bar (J/mol/K)

R          = 8.314               # Gas constant (J/mol/K) 
# G_Mg⁰    = -417093.5           # Gibbs free energy (J/mol) of Mg2+ at P and T from PerpleX
# G_SiO₂⁰  = -851437.5           # Gibbs free energy (J/mol) of SiO2(aq) at P and T from PerpleX

# Interpolate Gibbs free energy (J/mol) of Mg2+ and SiO2(aq) at P and T from PerpleX lookup table
G0_Mg_lt   = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_Mg2+_500C.csv", DataFrame))
G0_SiO2_lt = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_SiO2(aq)_500C.csv", DataFrame))
P_var      = collect(G0_Mg_lt[:,1])                                    # P values for lookup tables of Gibbs free energies for Mg2+ and SiO2(aq)
G0_Mg      = collect(G0_Mg_lt[:,2])                                    # Gibbs free energies for Mg2+ at different pressures (1 to 25 kbar)
G0_SiO2    = collect(G0_SiO2_lt[:,2])                                  # Gibbs free energies for SiO2(aq) at different pressures (1 to 25 kbar)
G_Mg⁰      = Itp1D_rev_scalar1(P_var, G0_Mg, 1000*P)      # G0 for Mg2+ at the pressure of interest
G_SiO₂⁰    = Itp1D_rev_scalar1(P_var, G0_SiO2, 1000*P)    # G0 for SiO2(aq) at the pressure of interest
@printf("Gibbs free energy of Mg2+ = %2.5e and SiO2(aq) = %2.5e\n", G_Mg⁰, G_SiO₂⁰)

# Compute the log(aM/aH+)
logMgH     = (1000*µ_MgO + S0_MgO*298.15 - 1000*µ_H₂O - G_Mg⁰) / (2.303 * R * (T_calc+273.15))
logSiO₂H   = (1000*µ_SiO₂ + S0_SiO₂*298.15 - G_SiO₂⁰) / (2.303 * R * (T_calc+273.15))
@printf("Log(aSiO2) = %2.5e and log(aMg2+/aH+) = %2.5e\n", logSiO₂H, logMgH)
logaoxides = [logMgH;logSiO₂H]

# Compute the speciation using the log(aMg2+/aH+) from MAGEMin
Speciation(logaoxides, T_calc, P)

