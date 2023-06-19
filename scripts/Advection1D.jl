using GLMakie
Makie.inline!(true)

function Itp1D_Centers2Markers( Tm, xm, Tc, xc, Δx, xmin )
    for i = 1:length(xm)
        iW    = trunc(Int, (xm[i] - xmin) / Δx) + 1 # index of west node
        Δxm   = xm[i] - xc[iW] # distance to marker-node
        wW    = 1.0 - Δxm/Δx
        wE    = Δxm/Δx
        # Your task: Update marker field as function of grid values and weights
        Tm[i] = wW*Tc[iW] + wE*Tc[iW+1]
    end
end

function Itp1D_Markers2Centers( Tc, xc, Tm, xm, Δx, xmin )
    Tc .= 0.0 
    Wc = zero(Tc)
    for i = 1:length(xm)
        iC      = trunc(Int, (xm[i] - xmin) / Δx) + 1
        w       = 1.0 - (xm[i] - xc[iC]) / Δx 
        # Your task: Update grid field values and weight array as function of marker values
        Tc[iC] += w*Tm[i]
        Wc[iC] += w
    end
    Tc ./= Wc
end

function main()

# Physics
xmin = -1/2
xmax = 1/2
Tmax = 100  # max. amplitude
σ    = 0.1  # bandwidth
ρ    = 1000 # density
c    = 1000 # heat capacity
Twest = 0.0
Teast = 0.0 
V     = -1.0

# Numerics
nx   = 100 # vertices
ncx  = nx-1 # centroids
nt   = 100  # time steps
nout = 10   # plot each nout

# Pre-compute
Δx   = (xmax-xmin)/ncx
xv   = LinRange(xmin, xmax, nx)
xc   = 0.5*(xv[1:end-1] .+ xv[2:end])
xce  = LinRange(xmin-Δx/2, xmax+Δx/2, ncx+2)  # New for marker interpolation
Vx   = V*ones(nx) 
Δt   = Δx/abs(V)/1.1

# Allocate arrays
Tup  = zeros(ncx)
Tmc  = zeros(ncx)
TWE  = zeros(ncx+2)
qx   = zeros(nx)
Tana = zeros(ncx)

# Initial condition
t    = 0.0
Tup .= Tmax.*exp.(-(xc.-V*t).^2/2/σ^2)
Tmc .= Tmax.*exp.(-(xc.-V*t).^2/2/σ^2)
Tini = copy(Tup)

# Marker data
nxm  = 4*ncx 
Δxm  = (xmax-xmin)/nxm
xm   = collect(LinRange(xmin+Δxm/2.0, xmax-Δxm/2.0, nxm))
Tm   = zeros(nxm)
Vxm  = V*ones(nxm) 

# Populate TWE for marker scheme
TWE[2:end-1] .= Tmc
TWE[1]        = Tmc[end]
TWE[end]      = Tmc[1]
# YOUR TASK: Finish this interpolation function (TOP OF FILE)
Itp1D_Centers2Markers( Tm, xm, TWE, xce, Δx, xmin )

# Time loop
for it=1:nt

    t += Δt

    # UPWIND: Populate TWE for upwind scheme
    TWE[2:end-1] .= Tup
    TWE[1]        = Tup[end]
    TWE[end]      = Tup[1]

    # UPWIND: Flux
    for i in eachindex(qx)  # i=1:nx, i=1:size(qx,1)
        if Vx[i]>0.0
            qx[i] = Vx[i] * TWE[i]
        end
        if Vx[i]<0.0
            qx[i] = Vx[i] * TWE[i+1]
        end
    end

    # UPWIND: Flux concise
    qx .= (Vx.>0.0) .* (Vx.*TWE[1:end-1]) .+ (Vx.<0.0) .* (Vx.*TWE[2:end])

    # UPWIND: update using conservation law in Eulerian form
    Tup .= Tup .- Δt .* diff(qx, dims=1) ./ Δx

    # MARKERS: update using conservation law in Lagrangian form
    xm .+= Vxm.*Δt   # similar to  xm .= xm .+ Vxm.*Δt
    
    # Set periodicity
    xm[xm.<xmin] .= xmax .- abs.(xm[xm.<xmin] .- xmin) # Periodic
    xm[xm.>xmax] .= xmin .+ abs.(xm[xm.>xmax] .- xmax) # Periodic
   
    # MARKERS: interpolate T from markers that were advected
    Itp1D_Markers2Centers( Tmc, xc, Tm, xm, Δx, xmin )

    # Visualise
    if mod(it,nout)==0 || it==1

        f = Figure(resolution = (1200,600), fontsize=25, aspect = 2.0)
        ax1 = Axis(f[1, 1], title = L"$$Advection/Transport", xlabel = L"$t$ [Ga]", ylabel = L"$$Continental crust volume [%]",  xreversed = false, xgridvisible = false, ygridvisible = false)
        # Interpolated data
        lines!(ax1, xc, Tini, label="T initial")
        lines!(ax1, xc, Tup,  label="T upwind")
        x0      = V*t
        fact = floor(abs.(x0-xmin)/(xmax-xmin))
        if x0<xmin x0 = xmax -  ( abs(x0 - xmin) - (fact)*(xmax-xmin))  end # Periodic
        fact = floor(abs.(x0-xmax)/(xmax-xmin))
        if x0>xmax x0 = xmin +  ( abs(x0 - xmax) - (fact)*(xmax-xmin))  end # Periodic
        Tana   .= Tmax.*exp.(-(xc.-x0).^2/2/σ^2)
        lines!(ax1, xc, Tana, label="T analytic")
        scatter!(ax1, xm[1:5:end], Tm[1:5:end], label="T on markers")
        lines!(ax1, xc, Tmc, label="T from markers")        
        f[1, 2] = Legend(f, ax1, "Legend", framevisible = false)
        # # Display figure
        DataInspector(f)
        display(f)
    end

end

end

main()

