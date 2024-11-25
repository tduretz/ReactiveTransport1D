function Res!(f, m, D)
    # I = 1
    I   = m[1]  # in practice, should be prop. to sum(m) 
    n   = ones(size(m))
    n  .= I
    for i in eachindex(f)
        f[i] = D[i,:]' * m +  D[i,:]' * n
    end
end

let 
    # Initialise
    D = ones(8,8)
    f = zeros(8)
    m = ones(8)

    # Residual
    Res!(f, m, D)
    @show f

    # Jacobian
    r = (f, m) -> Res!(f, m, D) 
    J = ForwardDiff.jacobian(r, f, m)
end