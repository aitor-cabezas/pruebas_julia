using LinearAlgebra
using SparseArrays

# Función que calcula du/dt y J del boundary layer test



function dudt_J_blt(ti,ui)
    
    #   Bucle nodos interiores
    
    for i=2:Nx-1, j=2:Ny-1
        
        N           =   (i-1)*Ny + j
        filas_int   =   vcat(filas_int,N,N,N,N,N)
        cols_int    =   vcat(cols_int,N,N+1,N-1,N+Ny,N-Ny)
        vals        =   vcat(vals,
                             -ui[N]*(1/hx+1/hy) - 2*eps*(1/hx^2 + 1/hy^2),
                             eps/hy^2,
                             ui[N]/hy + eps/hy^2,
                             eps/hx^2,
                             ui[N]/hx + eps/hx^2
                             )
        
    end
    
    #   Condiciones de Contorno
    
    #   Bucle Frontera Izquierda
        
    for j=1:Ny
        
        N           =   j
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
        ui[N]       =   1.0
        
    end
    
    #   Bucle Frontera Inferior
        
    for i=2:Nx
        
        N           =   (i-1)*Ny + 1
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
        ui[N]       =   1.0
        
    end
    
    
    #   Bucle Frontera Derecha
        
    for j=2:Ny
        
        N           =   (Nx-1)*Ny + j
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
        ui[N]       =   1.0 - y[j]
        
    end
    
    
    #   Bucle Frontera Superior
        
    for i=2:Nx-1
        
        N           =   i*Ny
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
        ui[N]       =   1.0 - x[i]
        
    end
    
    J               =   sparse(filas_int[:],cols_int[:],vals[:],Nx*Ny,Nx*Ny)
    
    
    
    phi             =   1 + 1/2*(1+tanh((tn-1/2)/delta_t))
    dphidt          =   1/2*((1-(tanh((tn-1/2)/delta_t))^2)/delta_t)
    
    for i=2:Nx-1,j=2:Ny-1
        
        N             =   (i-1)*Ny + j
        fx            =   x[i] + c1*exp((x[i]-1)/eps) + c2
        gy            =   y[j] + c1*exp((y[j]-1)/eps) + c2
        U             =   fx*gy
        dUdx          =   (1+c1/eps*exp((x[i]-1)/eps))*gy
        dUdy          =   (1+c1/eps*exp((y[j]-1)/eps))*fx
        d2Udx2        =   (1+c1/eps^2*exp((x[i]-1)/eps))*gy
        d2Udy2        =   (1+c1/eps^2*exp((y[j]-1)/eps))*fx
        Qij           =   U*dphidt +
                          U*phi^2*(dUdx + dUdy) -
                          eps*phi*(d2Udx2 + d2Udy2)
        dudt[N]       =   -ui[N]*((ui[N]-ui[N-Ny])/hx + 
                          (ui[N] - ui[N-1])/hy) +
                          eps*((ui[N+Ny] - 2*ui[N] + ui[N-Ny])/hx^2 + 
                          (ui[N+1] - 2*ui[N] + ui[N-1])/hy^2) +
                          Qij
                     
    end
    
       
    return dudt, Jdudt
 
 
end
    
