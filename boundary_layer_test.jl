using LinearAlgebra
using SparseArrays

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


function boundary_layer_test(; Nx::Int64 = 20, Ny::Int64 = 20, t0::Int64 = 0,tf::Int64 = 1,Nt::Int64 = 10, epsi::Float64 = 0.001)
    
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
    
    u           =   zeros(Float64,Nt,Nx*Ny) #   Matriz soluciones   
    
    #   Condición Inicial
    
    u0          =   ones(Nx*Ny)
    unm1        =   copy(u0)        #   Solución en el instante n
    
    for k=1:Nt
        
        un                  =   SDIRK(unm1,t[k],param)
        u[k,:]              =   copy(un)
        unm1                =   copy(un)
        
    end
    
    
    
    
    return u
    
    
    
end
