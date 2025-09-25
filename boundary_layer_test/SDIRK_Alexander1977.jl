
function SDIRK(unm1::Vector{Float64}, tnm1::Float64, param::Union{parametros, parametros_5})

    # Tabla de Butcher
    gamma = 1 - 1/sqrt(2)
    A = [gamma 0; 1 - 2*gamma gamma]
    b = [1/2, 1/2]
    c = [gamma, 1 - gamma]

    s = 2
    Nx = param.Nx
    Ny = param.Ny
    ht = param.ht

    ui = zeros(Float64, s, Nx*Ny)
    myI = sparse(1:Nx*Ny, 1:Nx*Ny, ones(Nx*Ny))
    tol = 1e-6
    iter = 50

    # Etapa 1
    function R1(u)
        dudt1 = dudt_blt(tnm1 + c[1]*ht, u, param)
        return u - unm1 - ht*A[1,1]*dudt1
    end

    function JR1(u)
        Jdudt1 = J_blt(u, param)
        return myI - ht*A[1,1]*Jdudt1
    end

    ui[1,:] = newton(R1, JR1, unm1, tol, iter)

    if Nx > 1 && Ny > 1
        uexact1 = uexact_blt(param.epsi, tnm1 + c[1]*ht, param.x, Nx, param.y, Ny)
        for j = 2:Ny-1
            ui[1,(j-1)*Nx + 1] = uexact1[(j-1)*Nx + 1]
        end
        for i = 1:Nx
            ui[1,i] = uexact1[i]
        end
        for j = 2:Ny
            ui[1,j*Nx] = uexact1[j*Nx]
        end
        for i = 1:Nx-1
            ui[1,Nx*(Ny-1) + i] = uexact1[Nx*(Ny-1) + i]
        end
    end

    # Etapa 2
    function R2(u)
        dudt1 = dudt_blt(tnm1 + c[1]*ht, ui[1,:], param)
        dudt2 = dudt_blt(tnm1 + c[2]*ht, u, param)
        return u - unm1 - ht*(A[2,1]*dudt1 + A[2,2]*dudt2)
    end

    function JR2(u)
        Jdudt2 = J_blt(u, param)
        return myI - ht*A[2,2]*Jdudt2
    end

    ui[2,:] = newton(R2, JR2, unm1, tol, iter)

    if Nx > 1 && Ny > 1
        uexact2 = uexact_blt(param.epsi, tnm1 + c[2]*ht, param.x, Nx, param.y, Ny)
        for j = 2:Ny-1
            ui[2,(j-1)*Nx + 1] = uexact2[(j-1)*Nx + 1]
        end
        for i = 1:Nx
            ui[2,i] = uexact2[i]
        end
        for j = 2:Ny
            ui[2,j*Nx] = uexact2[j*Nx]
        end
        for i = 1:Nx-1
            ui[2,Nx*(Ny-1) + i] = uexact2[Nx*(Ny-1) + i]
        end
    end

    dudt1 = dudt_blt(tnm1 + c[1]*ht, ui[1,:], param)
    dudt2 = dudt_blt(tnm1 + c[2]*ht, ui[2,:], param)
    un = unm1 + ht*(b[1]*dudt1 + b[2]*dudt2)

    if Nx > 1 && Ny > 1
        uexactn = uexact_blt(param.epsi, tnm1 + ht, param.x, Nx, param.y, Ny)
        for j = 2:Ny-1
            un[(j-1)*Nx + 1] = uexactn[(j-1)*Nx + 1]
        end
        for i = 1:Nx
            un[i] = uexactn[i]
        end
        for j = 2:Ny
            un[j*Nx] = uexactn[j*Nx]
        end
        for i = 1:Nx-1
            un[Nx*(Ny-1) + i] = uexactn[Nx*(Ny-1) + i]
        end
    end

    return un
end
