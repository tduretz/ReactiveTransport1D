using CairoMakie
function MainPorosityWaves1D()
    # Physics
    Ly   = 600.
    y0   = -Ly/4.0
    ϕmax = 1.0
    ϕbg  = 1.0
    η0   = 1.0
    λ    = 5.0
    k_μf = 1.0
    β    = 1e-2
    Δρg  = -1.0
    R    = 1.0
    t    = 0.0
    # Numerics
    ncy  = 802
    nt   = 1000000
    nout = 1000
    # Preprocessing
    Δy   = Ly/ncy
    yv   = LinRange(-Ly/2, Ly/2, ncy+1)
    P    = zeros(ncy+1)
    P0   = zeros(ncy+1)
    ϕ    = ϕbg .+ ϕmax.*exp.(-((yv.-y0/2)/λ).^2) .+ ϕmax/4.0.*exp.(-((yv.-y0/4)/λ).^2)
    ϕ0   = copy(ϕ)
    ηϕ   = η0.*ones(ncy+1)
    K    = k_μf .* ϕ.^3
    Kc   = 0.5*(K[1:end-1] .+ K[2:end-0])
    qD   = zeros(ncy)
    dϕdt = zeros(ncy+1)
    dPdt = zeros(ncy+1)
    # Time loop
    for it=1:nt
        # Time
        Δt = β * Δy^2 / maximum(Kc) / 4.5
        t += Δt
        # BC: periodic
        P[1] = P[end-1]
        P[end] = P[2]
        # Business
        @. ηϕ            = η0
        @. ηϕ[P>0.0]     = η0/R
        @. K             = k_μf * ϕ^3
        @. Kc            = 0.5*(K[1:end-1] .+ K[2:end-0])
        @. qD            = -Kc * ((P[2:end] - P[1:end-1])/Δy + Δρg)
        @. dϕdt          = P / ηϕ
        @. dPdt[2:end-1] = -((qD[2:end] - qD[1:end-1])/Δy +  dϕdt[2:end-1])/β 
        @. P            += Δt*dPdt
        @. ϕ            += Δt*dϕdt
        @. ϕ[ϕ<0.0] = 0.0
        if it%nout==0 || it==1
            f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
            ax1 = Axis(f[1, 1], title = L"$$Pressure", xlabel = L"$P$ [-]", ylabel = L"$y$ [-]",  xreversed = false, xgridvisible = false, ygridvisible = false)
            lines!(ax1, P, yv, label="Pressure")
            ax2 = Axis(f[1, 2], title = L"$$Porosity", xlabel = L"$ϕ$ [-]", ylabel = L"$y$ [-]",  xreversed = false, xgridvisible = false, ygridvisible = false)
            lines!(ax2, ϕ, yv, label="Porosity")
            ax3 = Axis(f[1, 3], title = L"$$Portrait", xlabel = L"$ϕ$ [-]", ylabel = L"$P$ [-]",  xreversed = false, xgridvisible = false, ygridvisible = false)
            lines!(ax3, ϕ, P, label="Portrait")
            display(f)
        end
    end
end

MainPorosityWaves1D()