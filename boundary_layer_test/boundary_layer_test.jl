using LinearAlgebra
using SparseArrays
using PyPlot

#   Boundary layer test main

#   Parámetros del problema

struct parametros  
        
        epsi::Float64       #   Coeficiente de Difusión
        x::Vector{Float64}
        Nx::Int64           #   Nodos eje x
        hx::Float64         #   Paso eje x
        y::Vector{Float64}
        Ny::Int64           #   Nodos eje y 
        hy::Float64         #   Paso eje y
        t0::Int64           #   Tiempo Inicial
        tf::Int64           #   Tiempo final
        Nt::Int64           #   Nodos temporales
        ht::Float64         #   Paso temporal


end

include("dudt_blt.jl")
include("J_blt.jl")
include("newton.jl")
include("SDIRK_Alexander1977.jl")
include("uexact_blt.jl")


function boundary_layer_test(; Nx::Int64 = 80, Ny::Int64 = 80, t0::Int64 = 0,tf::Int64 = 1,Nt::Int64 = 20, epsi::Float64 = 0.1, SC::Int=0)
    
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
    
    t           =   LinRange(t0,tf,Nt)
    ht          =   (tf-t0)/(Nt-1)
    
    
    param       =   parametros(epsi,x,Nx,hx,y,Ny,hy,t0,tf,Nt,ht)
    
    u           =   zeros(Float64,Nt,Nx*Ny) #Matriz soluciones aproximadas
    uexact      =   zeros(Float64,Nt,Nx*Ny) #Matriz de soluciones exactas
    errv        =   zeros(Float64,Nt)
    
    
    #   Condición Inicial
    
    u0          =   uexact_blt(epsi, t[1], x, Nx, y, Ny)
    u[1,:]      =   u0
    uexact[1,:] =   uexact_blt(epsi, t[1], x, Nx, y, Ny)
    errv[1]     =   maximum(abs.(uexact[1,:] .- u[1,:]))

    unm1        =   copy(u0)        #   Solución en el instante n
    
    for k=2:Nt
        un                  =   SDIRK(unm1,t[k-1],param)
        uex                 =   uexact_blt(epsi, t[k], x, Nx, y, Ny)

        # Frontera izquierda

        for j = 2:Ny-1
            N = (j-1)*Nx + 1
            un[N] = uex[N]
        end

        # Frontera inferior

        for i = 1:Nx
            N = i
            un[N] = uex[N]
        end

        # Frontera derecha

        for j = 2:Ny
            N = j*Nx
            un[N] = uex[N]
        end

        # Frontera superior

        for i = 1:Nx-1
            N = Nx*(Ny-1) + i
            un[N] = uex[N]
        end

        u[k,:]              =   copy(un)
        unm1                =   copy(un)
        uexact[k,:]         =   uexact_blt(epsi,t[k],x,Nx,y,Ny)
        errmax = 0.0
        for j = 2:Ny-1,i = 2:Nx-1
                N = (j-1)*Nx + i
                err = abs(uexact[k,N] - u[k,N])
                if err > errmax
                   errmax = err
                end
        end
        errv[k] = errmax
    end
    
    #   Representación de los resultados
        
    xmat        =   zeros(Nx,Ny)
    ymat        =   zeros(Nx,Ny)
    umat        =   zeros(Nx,Ny)
    umatexact   =   zeros(Nx,Ny)
    
    
    T       =   Nt
    
    for i=1:Nx,j=1:Ny
        N                =   (j-1)*Nx + i
        xmat[i,j]        =   x[i]
        ymat[i,j]        =   y[j]
        umat[i,j]        =   u[T,N]
        umatexact[i,j]   =   uexact[T,N]
    end
    
    figure()
    contourf(xmat, ymat, umat ,100 ,cmap="jet")
    colorbar()
    xlabel("x")
    ylabel("y")
    title("Solución aprox")
    
    figure()
    contourf(xmat, ymat, umatexact ,100 ,cmap="jet")
    colorbar()
    xlabel("x")
    ylabel("y")
    title("Solución analítica")

    figure()
    plot(t,errv)

    show(errv)
    
end
