using LinearAlgebra, GLMakie

# Solve for log(activities)
species = ["Mg2+"  "MgCl+"   "HCl"    "Cl-"    "H+"     "OH-"  "Mg(OH)₂" "H2O"]
# unknowns: log10 molalities of each specie

# 1: Chlorinity #Cl- per specie: 
# log10(m_MgCl⁺) + log10(m_HCl) + log10(m_Cl⁻) = log10(Clᵗᵒᵗ) 

# 2: Charges: neutrality  
# 2*log10(m_Mg2⁺) + log10(m_MgCl⁺) - log10(m_Cl⁻) + log10(m_H⁺) - log10(m_OH⁻) = 0. 

function f1(Sc)

    #    MgCl2   HCl   MgCl+  Mg2+    Cl-    H+     OH-
    M = [ 0.     0.     1.     2.    -1.     1.    -1.    ; # 1
          2.     1.     1.     0.     1.     0.     0.    ; # 2
          1.    -2.     0.     0.     0.     0.     0.    ;
         -1.     0.     1.     0.     1.     0.     0.    ;
          0.     0.    -1.     1.     1.     0.     0.    ;
          0.    -1.     0.     0.     1.     1.     0.    ; 
          0.     0.     0.     0.     0.     1.     1.    ;
        ]

    M[3:end,:] .*= Sc
    @show M

    # log(K) of reactions P = 5 kbar, T = 300 °C 
    # lK    = [-0.02; 0.95;  3.043; 0.4823; -1.3025; -.6277; -6.2149;]

    lK    = [0; 1; 5.04; -1.68; -3.62; -2.79; -10.2]

    # Solve for activities log10(a_i)
    lai   = M\lK

    # # mi 
    # mi = 0.01*ones(length(lK))
    # for iter=1:1
    #     @show f  = M*mi - lK
    #     @show norm(f)
    #     mi .+= -M\f
    # end
    # @show mi

    mi = 0.01*ones(length(lK)) 
    
    for iter=1:10000
        f  = zero(mi)
           #    MgCl2   HCl     MgCl+   Mg2+    Cl-     H+      OH-
        f[1] = [ 0.;     0.;     1.;     2.;    -1.;     1.;    -1.]' * mi - lK[1]
        f[2] = [ 2.;     1.;     1.;     0.;     1.;     0.;     0.]' * mi - lK[2]
        f[3] = [ 1.;    -2.;     0.;     0.;     0.;     0.;     0.;]' * log10.((mi)) - lK[3]
        f[4] = [-1.;     0.;     1.;     0.;     1.;     0.;     0.;]' * log10.((mi)) - lK[4]
        f[5] = [ 0.;     0.;    -1.;     1.;     1.;     0.;     0.;]' * log10.((mi)) - lK[5]
        f[6] = [ 0.;    -1.;     0.;     0.;     1.;     1.;     0.;]' * log10.((mi)) - lK[6]
        f[7] = [ 0.;     0.;     0.;     0.;     0.;     1.;     1.;]' * log10.((mi)) - lK[7]
        
        println("it.", iter, " --- f = ", norm(f))
        if norm(f)<1e-12 break end

        # J = [ 0.     0.     1.     2.    -1.     1.    -1.    ; # 1
        #       2.     1.     1.     0.     1.     0.     0.    ; # 2
        #       1.0/(mi[1]*1)    -2/(mi[2]*1)     0.     0.     0.     0.     0.    ;
        #       -1.0/(mi[1]*1)     0.     1/(mi[3]*1)     0.     1/(mi[5]*1)     0.     0.    ;
        #       0.     0.    -1/(mi[3]*1)     1/(mi[4]*1)     1/(mi[5]*1)     0.     0.    ;
        #       0.    -1/(mi[2]*1)     0.     0.     1/(mi[5]*1)     1/(mi[6]*1)     0.    ; 
        #       0.     0.     0.     0.     0.     1/(mi[6]*1)     1/(mi[7]*1)    ;
        # ]

        J = [ 0.     0.     1.     2.    -1.     1.    -1.    ; # 1
              2.     1.     1.     0.     1.     0.     0.    ; # 2
              1.0/(mi[1]*log(10))    -2/(mi[2]*log(10))     0.     0.     0.     0.     0.    ;
              -1.0/(mi[1]*log(10))     0.     1/(mi[3]*log(10))     0.     1/(mi[5]*log(10))     0.     0.    ;
              0.     0.    -1/(mi[3]*log(10))     1/(mi[4]*log(10))     1/(mi[5]*log(10))     0.     0.    ;
              0.    -1/(mi[2]*log(10))     0.     0.     1/(mi[5]*log(10))     1/(mi[6]*log(10))     0.    ; 
              0.     0.     0.     0.     0.     1/(mi[6]*log(10))     1/(mi[7]*log(10))    ;
        ]

        Δm = J\f
        if iter<20
            mi .-= 0.1*Δm
        else
            mi .-= Δm
        end

    end


    return mi

end

Sc = [1e0]

f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
ax1 = Axis(f[1, 1], title = L"$$z", xlabel = L"$$Species", ylabel = L"$$Activity [-]",  xreversed = false, xgridvisible = false, ygridvisible = false)

for i in eachindex(Sc)
    m = f1(Sc[i])
    @show m
    lines!(ax1, 1:length(m), m, label="a")
end

display(f)