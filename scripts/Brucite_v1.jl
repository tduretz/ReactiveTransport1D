using LinearAlgebra, GLMakie

# Solve for log(activities)
species = ["Mg2+"  "MgCl+"   "HCl"    "Cl-"    "H+"     "OH-"  "Mg(OH)₂" "H2O"]
# unknowns: log10 molalities of each specie

# Results from EQ3: 
#    Mg2+     MgCl+    HCl     Cl-     H+      OH-      Mg(OH)₂  H₂O
#    -2.0371  -2.3170  -5.967  -1.1975 -4.9862 -4.3028  -.0994   -0.00142

# 1: Chlorinity #Cl- per specie: 
# m_MgCl⁺ + m_HCl + m_Cl⁻ = m_Clᵗᵒᵗ 

# 2: Charges: neutrality  
# 2*(m_Mg2⁺) + (m_MgCl⁺) - (m_Cl⁻) + (m_H⁺) - (m_OH⁻) = 0. 

function f1(Sc)

    #    Mg2+  MgCl+   HCl    Cl-    H+     OH-    Mg(OH)₂  H₂O
    M = [ 0.     1.     1.     1.     0.     0.     0.       0.; # 1
          2.     1.     0.    -1.     1.    -1.     0.       0.; # 2
          0.     1.    -1.     0.     0.     1.     0.       1.;
          0.     0.     0.     0.     1.     1.     0.      -1.;
          1.    -1.     0.     1.     0.     0.     0.       0.;
          0.     0.    -1.     1.     1.     0.     0.       0.; 
          1.     0.     0.     0.    -2.     0.    -1.       0.;
          0.     0.     0.     0.     0.     0.     0.       1.; # a_H20 = 1.0 or  log10(a_H20) = 0.0
        ]

    # M[3:end,:] .*= Sc

    # # log(K) of reactions P = 5 kbar, T = 300 °C 
    # Clᵗᵒᵗ = 0.1
    # #         Cl      Z     Br      H20      MgCl+   HCl     Mg(OH)₂ aH₂O
    # lK    = [Clᵗᵒᵗ; 0.0;  7.9325; -9.2876; -.9176; -.2167; 8.0319; 0.0;]

    # # Solve for activities log10(a_i)
    # @show lai   = M\lK
    # mi    = 10.0.^(lai)

    return mi

end

Sc = [1e0 ]

f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
ax1 = Axis(f[1, 1], title = L"$$z", xlabel = L"$t$ [Ga]", ylabel = L"$$Continental crust volume [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)

for i in eachindex(Sc)
    a = f1(Sc[i])
    lines!(ax1, 1:length(a), a, label="a")


end

@show norm(M*lai .- lK)