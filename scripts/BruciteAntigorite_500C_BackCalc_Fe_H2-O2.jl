using LinearAlgebra, CairoMakie, Printf, ForwardDiff, CSV, DataFrames, MathTeXEngine
using MAGEMin_C
Makie.inline!(true)
Makie.update_theme!(fonts=(regular=texfont(), bold=texfont(:bold), italic=texfont(:italic)))

# Comparison with fluid speciation from EQ3 using the DEW thermodynamic database
# Jacobian constructed with automatic differentiation

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
    df = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/Matrix_500C_BackCalc_Fe_H2-O2.csv", DataFrame))
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
    ϵ = 1e-10                                      # Non-linear tolerance
    iter = 0                                       # Iteration count

    # Log Ks lookup tables at variable pressure and 400 °C
    MgCl          = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MgCl+_dissociation_500C.csv", DataFrame))
    P_logK        = collect(MgCl[:,1])                 # Pressure vector for lookup tables of log K values
    LogK_MgCl     = collect(MgCl[:,2])                 # Log K values for MgCl+ dissociation
    HCl           = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/HCl_dissociation_500C.csv", DataFrame))
    LogK_HCl      = collect(HCl[:,2])                  # Log K values for HCl dissociation
    H2O           = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/H2O_dissociation_500C.csv", DataFrame))
    LogK_H2O      = collect(H2O[:,2])                  # Log K values for H2O dissociation
    MgOH          = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/MgOH+_dissociation_500C.csv", DataFrame))
    LogK_MgOH     = collect(MgOH[:,2])                 # Log K values for MgOH dissociation
    FeCl          = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/FeCl+_dissociation_500C.csv", DataFrame))
    LogK_FeCl     = collect(FeCl[:,2])                 # Log K values for FeCl+ dissociation
    FeCl2         = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/FeCl2_dissociation_500C.csv", DataFrame))
    LogK_FeCl2    = collect(FeCl2[:,2])                # Log K values for FeCl+ dissociation
    FeOH          = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/FeOH+_dissociation_500C.csv", DataFrame))
    LogK_FeOH     = collect(FeOH[:,2])                 # Log K values for MgOH dissociation
    H2_H2O        = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/H2(aq)-H2O_dissociation_500C.csv", DataFrame))
    LogK_H2_H2O   = collect(H2_H2O[:,2])               # Log K values for MgOH dissociation

    # Arrays
    D = Float64.(Matrix(df[:, 2:end-1]))                       # Read coefficients and convert to Matrix
    z = collect(df[1, 2:end-1])                                # Charges for the fluid species
    coeffH2O = collect(df[1:end, end])
    b = 0.01 * ones(length(coeffH2O))
    b[3] = logaoxides[1]
    b[4] = logaoxides[2]
    b[5] = logaoxides[3]
    b[6] = logaoxides[4]
    b[7] = Itp1D_rev_scalar1(P_logK, LogK_MgCl, 1000*P)        # Interpolate log Ks for MgCl+ dissociation constant
    b[8] = Itp1D_rev_scalar1(P_logK, LogK_HCl, 1000*P)         # Interpolate log Ks for HCl dissociation constant
    b[9] = Itp1D_rev_scalar1(P_logK, LogK_H2O, 1000*P)         # Interpolate log Ks for H2O dissociation constant
    b[10] = Itp1D_rev_scalar1(P_logK, LogK_MgOH, 1000*P)       # Interpolate log Ks for MgOH+ dissociation constant
    b[11] = Itp1D_rev_scalar1(P_logK, LogK_FeCl2, 1000*P)      # Interpolate log Ks for MgOH+ dissociation constant
    b[12] = Itp1D_rev_scalar1(P_logK, LogK_FeCl, 1000*P)       # Interpolate log Ks for MgOH+ dissociation constant
    b[13] = Itp1D_rev_scalar1(P_logK, LogK_FeOH, 1000*P)       # Interpolate log Ks for MgOH+ dissociation constant
    b[14] = Itp1D_rev_scalar1(P_logK, LogK_H2_H2O, 1000*P)     # Interpolate log Ks for MgOH+ dissociation constant
    å = 3.7 * ones(length(b))                      # Size of fluid species (including hydration shell)
    # b    = [0.0; Clᵗᵒᵗ; 6.8466; -1.0841; -0.6078; -8.1764; 6.6296; 4.9398]  
    n = length(species)
    m = 0.01 * ones(length(b))                               # Initial condition
    m[6] = 10^(b[9]/2)                                       # Initial condition for H+ to be half of Ke
    m[7] = 10^(b[9]/2)                                       # Initial condition for OH- to be half of Ke
    m[14] = 10^(b[6])                                       # Initial condition for O2 to be 10^(log fO2)
    # m = [0.02; 0.01; 0.01; 0.1; 0.0005; 0.0008; 0.3; 0.01; 0.01]
    logγ = ones(length(b))                                   # Initial activity coefficients equal to 1
    f = zero(m)
    Δm = zero(m)
    J = zeros(n, n)
    e = zeros(niter)
    logaw = 0                                                # Initial log of water activity equal to 0 (activity of water equal to 1)
    for i in 1:n
        @printf("Initial conditions for molality of %s is %1.4e\n", species[i], m[i])
    end

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
    @printf("Mg dissolved is %1.2e\n", m[2]+m[3]+m[8])
    @printf("Fe dissolved is %1.2e\n", m[10]+m[11]+m[12])
    @printf("Si dissolved is %1.2e\n", m[9])
    @printf("Log fO2 is %1.2f\n", log10(m[14]))

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
    out = single_point_minimization(P, T_calc, data, X=X, Xoxides=Xoxides, sys_in=sys_in; name_solvus=true)
    @printf("Deviation residuals is %1.10f\n", out.bulk_res_norm)
    for n in 1:length(out.ph)
        @printf("%s is stable with %1.3f vol\n", out.ph[n],out.ph_frac_vol[n])
    end
    for j in 1:length(out.sol_name)-1
        @printf("%s is stable with %1.3f vol\n", out.sol_name[j],out.ph_frac_vol[j])
        print(out.SS_vec[j].emNames)
        @printf(" \n")
        print(out.SS_vec[j].emFrac)
        @printf(" \n")
        # @print(out.SS_vec[j].siteFractionsNames)
        # @print("The site fractions proportions are %s is stable with %2.2e\n", out.SS_vec[j].siteFractions)
    end
    for m in 1:length(Xoxides)
        @printf("The chemical potential of %s = %2.5e (kJ/mol)\n", out.oxides[m],out.Gamma[m])
    end
    # print(out.SS_vec[6].siteFractionsNames)
    # print(out.SS_vec[6].siteFractions)

    return out.Gamma[1], out.Gamma[4], out.Gamma[3], out.Gamma[6], out.Gamma[5]

end

# Initialization for the MAGEMin conditions
P = 5.0
T_calc = 500.0
data = Initialize_MAGEMin("ume", verbose=false; solver=2);
Xoxides = ["SiO2"; "FeO"; "MgO"; "H2O"; "Al2O3"; "O"; "S"];  # System of reduced serpentinite of Evans & frost (2021)
X_comp = [34.146613; 6.415533; 33.41302; 23.883372; 1.808672; 0.060068; 0.272721];   # Composition of reduced serpentinite of Evans & frost (2021) in wt.%
sys_in = "wt"
# Xoxides = ["SiO2"; "FeO"; "MgO"; "H2O"];  # System of component for simple Olivine + H2O
# X_comp = [22.74; 4.15; 38.33; 34.78];   # Composition of simple Olivine + H2O in mol%
# sys_in = "mol"

# Get the chemical potentials values for each component
µ_SiO₂, µ_FeO, µ_MgO, µ_H₂O, µ_O = GetChemicalPotentials(X_comp, Xoxides, data, T_calc, P, sys_in)
@printf("Chemical potentials (kJ/mol) of SiO2 = %2.5e, FeO = %2.5e, MgO = %2.5e, of H2O = %2.5e and of O = %2.5e\n", µ_SiO₂, µ_FeO, µ_MgO, µ_H₂O, µ_O)

# Thermodynamic properties
# Entropies needed to convert HSC convention Gibbs free energies/chemical potentials to SUPCRT convention
S0_SiO₂    = 223.96              # Entropy of SiO2 at 298 K and 1 bar (J/mol/K)
S0_MgO     = 135.255             # Entropy of MgO at 298 K and 1 bar (J/mol/K)
S0_H₂O     = 233.255             # Entropy of H2O at 298 K and 1 bar (J/mol/K)
S0_FeO     = 129.855             # Entropy of FeO at 298 K and 1 bar (J/mol/K)
S0_O       = 129.855             # Entropy of O at 298 K and 1 bar (J/mol/K)

R          = 8.314               # Gas constant (J/mol/K) 
# G_Mg⁰ = -417093.5       # Gibbs free energy (J/mol) of Mg2+ at P and T from PerpleX
# G_Fe⁰ = -63615.6        # Gibbs free energy (J/mol) of Fe2+ at P and T from PerpleX
# G_SiO₂⁰ = -851437.5     # Gibbs free energy (J/mol) of SiO2(aq) at P and T from PerpleX

# Interpolate Gibbs free energy (J/mol) of Mg2+ and SiO2(aq) at P and T from PerpleX lookup table
G0_Mg_lt   = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_Mg2+_500C.csv", DataFrame))
G0_SiO2_lt = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_SiO2(aq)_500C.csv", DataFrame))
G0_Fe_lt   = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_Fe2+_500C.csv", DataFrame))
G0_O2_lt   = (CSV.read("/Users/guillaumesiron/Documents/Julia_scripts/ReactiveTransport1D/data/G0_O2_500C.csv", DataFrame))
P_var      = collect(G0_Mg_lt[:,1])                          # P values for lookup tables of Gibbs free energies for Mg2+ and SiO2(aq)
G0_Mg      = collect(G0_Mg_lt[:,2])                          # Gibbs free energies for Mg2+ at different pressures (1 to 25 kbar)
G0_SiO2    = collect(G0_SiO2_lt[:,2])                        # Gibbs free energies for SiO2(aq) at different pressures (1 to 25 kbar)
G0_Fe      = collect(G0_Fe_lt[:,2])                          # Gibbs free energies for Mg2+ at different pressures (1 to 25 kbar)
G0_O2      = collect(G0_O2_lt[:,2])                          # Gibbs free energies for O2(g) at different pressures (1 to 25 kbar)
G_Mg⁰      = Itp1D_rev_scalar1(P_var, G0_Mg, 1000*P)         # G0 for Mg2+ at the pressure of interest
G_SiO₂⁰    = Itp1D_rev_scalar1(P_var, G0_SiO2, 1000*P)       # G0 for SiO2(aq) at the pressure of interest
G_Fe⁰      = Itp1D_rev_scalar1(P_var, G0_Fe, 1000*P)         # G0 for Mg2+ at the pressure of interest
G_O₂⁰      = Itp1D_rev_scalar1(P_var, G0_O2, 1000*P)         # G0 for Mg2+ at the pressure of interest
@printf("Gibbs free energy of Mg2+ = %2.5e, SiO2(aq) = %2.5e, Fe2+ = %2.5e and O2(g) = %2.5e\n", G_Mg⁰, G_SiO₂⁰, G_Mg⁰, G_O₂⁰)

# µ_SiO₂ = -968.919674    # SiO2 chemical potential (kJ/mol) from PerpleX with HSC convention
# µ_MgO = -636.3732783    # MgO chemical potential (kJ/mol) from PerpleX with HSC convention
# µ_H₂O = -335.0539783    # H2O chemical potential (kJ/mol) from PerpleX with HSC convention
# µ_FeO = -310.579        # FeO chemical potential (kJ/mol) from PerpleX with HSC convention
# µ_FeO = -310.9142683    # FeO chemical potential (kJ/mol) from PerpleX with HSC convention at 400 °C
# µ_FeO = -292.899 - S0_FeO*298.15        # FeO chemical potential (kJ/mol) from PerpleX with HSC convention at 500 °C
# µ_SiO₂  = -896345        # SiO2 chemical potential (kJ/mol) from PerpleX with SUPCRT convention at 500 °C
# µ_MgO   = -609284        # MgO chemical potential (kJ/mol) from PerpleX with SUPCRT convention at 500 °C
# µ_H₂O   = -262845        # H2O chemical potential (kJ/mol) from PerpleX with SUPCRT convention at 500 °C
# µ_FeO   = -292899        # FeO chemical potential (kJ/mol) from PerpleX with SUPCRT convention at 500 °C
# µ_O     = -257654        # O chemical potential (kJ/mol) from PerpleX with SUPCRT convention at 500 °C

# With HSC convention and chemical potentials from MAGEMin
logMgH = (1000*µ_MgO + S0_MgO*298.15 - 1000*µ_H₂O - S0_H₂O*298.15 - G_Mg⁰) / (2.303 * R * (T_calc+273.15))
logFeH = (1000*µ_FeO + S0_FeO*298.15 - 1000*µ_H₂O - S0_H₂O*298.15- G_Fe⁰) / (2.303 * R * (T_calc+273.15))
logSiO₂H = (1000*µ_SiO₂ + S0_SiO₂*298.15 - G_SiO₂⁰) / (2.303 * R * (T_calc+273.15))
logO₂ = (2*(1000*µ_O + S0_O*298.15) - G_O₂⁰) / (2.303 * R * (T_calc+273.15))

# # With SUPCRT convention (from PerpleX)
# logMgH = (µ_MgO - µ_H₂O - G_Mg⁰) / (2.303 * R * (T_calc+273.15))
# logFeH = (µ_FeO - µ_H₂O - G_Fe⁰) / (2.303 * R * (T_calc+273.15))
# logSiO₂H = (µ_SiO₂ - G_SiO₂⁰) / (2.303 * R * (T_calc+273.15))
# logO₂ = (2*µ_O - G_O₂⁰) / (2.303 * R * (T_calc+273.15))

@printf("Log(aSiO2) = %2.5e, Log(fO2) = %2.5e, log(aFe2+/aH+) = %2.5e and log(aMg2+/aH+) = %2.5e\n", logSiO₂H, logO₂, logFeH, logMgH)
logaoxides = [logMgH; logSiO₂H; logFeH; logO₂]

# Compute the speciation using the log(aMg2+/aH+) from MAGEMin
Speciation(logaoxides, T_calc, P)

