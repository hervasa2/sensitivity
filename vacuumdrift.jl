function driftcharge!(path, charge)
    q = 0
    e = 1.602176634e-19
    m = 9.1093837015e-31
    δt = 1e-10 #s
    if charge == "h+"
        q = e
    end
    if charge == "e-"
        q = -1*e
    end
    vel = zeros(size(path))
    for i in 2:size(path)[1]
        a = q*getEfield(path[i-1,:])/m
        vel[i,:] = vel[i-1,:] + a*δt
        δx = 0.5*a*δt^2 #+ vel[i-1,:]*δt
        path[i,:] = path[i-1,:] + δx
        if path[i,1] >= 0.0348 || path[i,2] >= 0.0398
            path[i,:] = path[i-1,:]
        end
    end
    return vel
end

function signal(path, charge)
    if charge == "h+"
        q = 1
    end
    if charge == "e-"
        q = -1
    end
    Δϕ_h = zeros(size(path)[1])
    t = δt*collect(0:length(Δϕ_h))
    Q = zeros(length(Δϕ_h))
    R = zeros(length(Δϕ_h))
    m = 0
    rep = 0
    for i in 2:length(Δϕ_h)
        Δϕ_h[i] = getWpot(path[i,:])-getWpot(path[i-1,:])
        Q[i] = Q[i-1] + q*Δϕ_h[i]
        if Q[i] == Q[i-1]
            rep += 1
        else
            rep = 0
        end
        R[i] = rep
        if rep == 0
            m = ( Q[i] - Q[i-convert(Int64,R[i-1]+1)] ) / ( R[i-1]+1 )  #interpolation
            for j in 1:convert(Int64,R[i-1])
                Q[i-j] = Q[i] - m*j
            end
        end
    end
    return t,Q
end
