using LinearAlgebra, GLMakie, Printf
Makie.inline!(false)

# Jacobian constructed manually

function Residual!(f, m, lK)
    #       MgCl2   HCl     MgCl+   Mg2+    Cl-     H+      OH-
    f[1] = [ 0.;     0.;     1.;     2.;    -1.;     1.;    -1.]' *         m   - lK[1]
    f[2] = [ 2.;     1.;     1.;     0.;     1.;     0.;     0.]' *         m   - lK[2]
    f[3] = [ 1.;    -2.;     0.;     0.;     0.;     0.;     0.]' * log10.((m)) - lK[3]
    f[4] = [-1.;     0.;     1.;     0.;     1.;     0.;     0.]' * log10.((m)) - lK[4]
    f[5] = [ 0.;     0.;    -1.;     1.;     1.;     0.;     0.]' * log10.((m)) - lK[5]
    f[6] = [ 0.;    -1.;     0.;     0.;     1.;     1.;     0.]' * log10.((m)) - lK[6]
    f[7] = [ 0.;     0.;     0.;     0.;     0.;     1.;     1.]' * log10.((m)) - lK[7]
end

function LineSearch(f, m, Δm, lK)

    α  = [0.1 0.25 0.5 .8 1.0]
    F  = zero(α)
    m1 = copy(m)

    for i in eachindex(α)
        m1 .= m + α[i].*Δm
        if minimum(m1) < 0. 
            f .= 1000.
        else
            Residual!(f, m1, lK )
        end
        F[i] = norm(f)
    end

    _, ind = findmin(F)

    return α[ind]

end

# Solve for log(activities)
species = ["Mg2+",  "MgCl+",   "HCl",    "Cl-",    "H+",     "OH-",  "Mg(OH)₂"]
# unknowns: log10 molalities of each specie

# 1: Chlorinity #Cl- per specie: 
# log10(m_MgCl⁺) + log10(m_HCl) + log10(m_Cl⁻) = log10(Clᵗᵒᵗ) 

# 2: Charges: neutrality  
# 2*log10(m_Mg2⁺) + log10(m_MgCl⁺) - log10(m_Cl⁻) + log10(m_H⁺) - log10(m_OH⁻) = 0. 

function Speciation()

    lK   = [0.; 1.0; 5.04; -1.68; -3.62; -2.79; -10.2]
    n    = length(species)
    m    = 0.01*ones(length(lK)) 
    f    = zero(m)
    dfdm = zero(m)
    J    = zeros(n,n)
    
    for iter=1:10
        
        # Residual evaluation
        Residual!(f, m, lK )

        if norm(f)<1e-12 break end

        dfdm = 1.0 ./ (m*log(10))
        J = [ 0.     0.     1.     2.    -1.     1.    -1.    ; # 1
              2.     1.     1.     0.     1.     0.     0.    ; # 2
              1.0   -2.     0.     0.     0.     0.     0.    ;
             -1.0    0.     1.     0.     1.     0.     0.    ;
              0.     0.    -1.     1.     1.     0.     0.    ;
              0.    -1.     0.     0.     1.     1.     0.    ; 
              0.     0.     0.     0.     0.     1.     1.    ;
        ]
        J[3:end,:] .*= repeat(dfdm, 1, n-2)'

        Δm  = -J\f
        α   = LineSearch(f, m, Δm, lK)
        m .+= α*Δm
        @printf("it. %03d --- f = %1.4e --- α = %1.2f\n", iter, norm(f), α)
    end


    return m

end

f   = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0 )
ax1 = Axis(f[1, 1], title = L"$$Molalities", xlabel = L"$$Species", ylabel = L"$$Molality [-]",  xreversed = false, xgridvisible = false, ygridvisible = false)
ax1.xticks = ( collect(1:length(species)), species)
m   = Speciation()
scatter!(ax1, 1:length(m), m, label="a")
DataInspector(f)
display(f)