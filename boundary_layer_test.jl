

#   Boundary layer test main

#   Parámetros del problema
        
    eps         =   0.001        #   Coeficiente de Difusión
    # Coeficientes separación de variables U(x)*phi(t) 
    c1          =   1/(exp(-1/eps)-1) 
    c2          =   -c1/2*(1+exp(-1/eps))-1/2
    delta_t     =   0.2         
    
    #   Discretización Espacial (x1x2 x y1y2)

    x1          =   0   # (0,0)
    x2          =   1   # (1,0)
    x           =   Vector(LinRange(x1,x2,Nx))
    hx          =   (x2-x1)/(Nx-1)  #Paso eje x
    y1          =   0   #(0,1)
    y2          =   1   #(1,1)
    y           =   Vector(LinRange(y1,y2,Ny))
    hy          =   (y2-y1)/(Ny-1)  #Paso eje y
    
    #   Discretización Temporal
    
    t   =   LinRange(t0,tf,Nt)
    ht  =   (tf-t0)/(Nt-1)
    
    #   Inicialización vector dudt (la f) y vectores para J(dudt) 
    
    dudt        =   zeros(Nx*Ny)
    filas_int   =   zeros(Int64,0,1)
    cols_int    =   zeros(Int64,0,1)
    vals        =   zeros(Float64,0,1)
    
    #   Condición Inicial
    
    un          =   ones(Nx*Ny)
