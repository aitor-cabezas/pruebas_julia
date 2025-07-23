
#   SDIRK

function SDIRK(unm1::Vector{Float64}, tnm1::Float64, param::parametros)
    
    #   Tabla de Butcher
    
    gamma   =   1 - 1/sqrt(2)
    
    A       =   [gamma              0;
                1-2*gamma       gamma]
    b       =   [1/2    1/2]
    c       =   [gamma  1-gamma]
    
    
    
    s       =   2                                           #   Número de etapas
    
    Nx      =   param.Nx
    Ny      =   param.Ny
    ht      =   param.ht
    
    ui      =   zeros(Float64,s,Nx*Ny)                      #   Solución etapa i
    
    myI = sparse(1:Nx*Ny, 1:Nx*Ny, ones(Nx*Ny))

    
    tol     =   0.001                                       #  Tolerancia Newton
    iter    =   50                                          #  Número de iteraciones máximas Newton
    
    #   Etapa 1 Runge-Kutta

    function R1(u)   
        dudt1        =   dudt_blt(tnm1 + c[1]*ht,u,param)
        return  u - unm1 - ht*A[1,1]*dudt1
    end

    function JR1(u) 
        Jdudt1       =   J_blt(u,param)
        return  myI - ht*A[1,1]*Jdudt1
    end

    ui[1,:]          =   newton(R1,JR1,unm1,tol,iter)

    #   Etapa 2 Runge-Kutta

    function R2(u)
        dudt1            =   dudt_blt(tnm1 + c[1]*ht,ui[1,:],param)
        dudt2            =   dudt_blt(tnm1 + c[2]*ht,u,param)
        return  u - unm1 - ht*(A[2,1]*dudt1 + A[2,2]*dudt2)
    end

    function JR2(u) 
        Jdudt2           =   J_blt(u,param)
        return  myI - ht*A[2,2]*Jdudt2
    end

    ui[2,:]          =   newton(R2,JR2,unm1,tol,iter)
    dudt1            =   dudt_blt(tnm1 + c[1]*ht,ui[1,:],param)
    dudt2            =   dudt_blt(tnm1 + c[2]*ht,ui[2,:],param)
    un               =   unm1 + ht*(b[1]*dudt1 + b[2]*dudt2)    #Solución en el instante n+1

    return un
       
    
end
