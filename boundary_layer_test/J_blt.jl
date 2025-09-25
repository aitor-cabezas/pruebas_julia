
# Expresión analítica Jacobiana boundary layer test

function J_blt(ui::Vector{Float64},param::parametros)
    
    epsi        =   param.epsi
    x           =   param.x
    Nx          =   param.Nx
    hx          =   param.hx
    y           =   param.y
    Ny          =   param.Ny
    hy          =   param.hy
    
    
    filas_int   =   zeros(Int64,0,1)
    cols_int    =   zeros(Int64,0,1)
    vals        =   zeros(Float64,0,1)
    
    #   Bucle nodos interiores
    
    for i=2:Nx-1, j=2:Ny-1
        
        N           =   (j-1)*Nx + i
        filas_int   =   vcat(filas_int,N,N,N,N,N)
        cols_int    =   vcat(cols_int,N,N+Nx,N-Nx,N+1,N-1)
        vals        =   vcat(vals,
                             -2*ui[N]*(1/hx+1/hy) +   ui[N-Nx]/hy + ui[N-1]/hx -
                             2*epsi*(1/hx^2 + 1/hy^2),
                             epsi/hy^2,
                             ui[N]/hy + epsi/hy^2,
                             epsi/hx^2,
                             ui[N]/hx + epsi/hx^2
                             )
        
    end
    
    #   Condiciones de Contorno
    
    #   Bucle Frontera Izquierda
        
    for j=2:Ny-1
        
        N           =   (j-1)*Nx + 1
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
        
    end
    
    #   Bucle Frontera Inferior
        
    for i=1:Nx
        
        N           =   i
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
    end
    
    
    #   Bucle Frontera Derecha
        
    for j=2:Ny
        
        N           =   j*Nx
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
    end
    
    
    #   Bucle Frontera Superior
        
    for i=1:Nx-1
        
        N           =   Nx*(Ny-1) + i
        filas_int   =   vcat(filas_int,N)
        cols_int    =   vcat(cols_int,N)
        vals        =   vcat(vals,1.0)
    end
    
    J_blt   =   sparse(filas_int[:],cols_int[:],vals[:],Nx*Ny,Nx*Ny)
    
    return J_blt
    
end
