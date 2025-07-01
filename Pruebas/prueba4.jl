using LinearAlgebra
using SparseArrays
using PyPlot
using LaTeXStrings

function prueba4(Nx::Int64,Ny::Int64)
    #Problema convección-difusión diferencias finitas 2-D

    t1  =   time()

    #Dominio espacial (x1x2 x y1y2)

    x1  =   0 # (0,0)
    x2  =   1 # (Nx,0)
    x   =   Vector(LinRange(x1,x2,Nx))
    hx  =   (x2-x1)/(Nx-1) #Paso eje x
    y1  =   0 # (0,Ny)
    y2  =   1 # (Nx,Ny)
    y   =   Vector(LinRange(y1,y2,Ny))
    hy  =   (y2-y1)/(Ny-1) #Paso eje y



    #Parámetros del problema
    c1  =   1.0 #Coeficiente término convectivo en el eje x
    c2  =   1.0 #Coeficiente término convectivo en el eje y
    e   =   0.01 #Coeficiente término difusión

    #Discretización espacial

    filas_int   =   zeros(Int64,0,1)    #Vector de filas nodos internos
    cols_int    =   zeros(Int64,0,1)    #Vector de columnas nodos internos
    vals_int    =   zeros(Float64,0,1)    #Vector de valores nodos internos

    #Bucle para los nodos interiores
    for i=2:Nx-1, j=2:Ny-1
        N   =   (i-1)*Ny+j
        filas_int   =   vcat(filas_int,N,N,N,N,N)
        cols_int    =   vcat(cols_int,N,N+1,N-1,N+Ny,N-Ny)
        vals_int    =   vcat(vals_int,
                             (c1/hx+c2/hy)+2*e*(1/hx^2+1/hy^2),
                             -e/hy^2,
                             -c2/hy-e/hy^2,
                             -e/hx^2,
                             -c1/hx-e/hx^2)
    end

    bv  =   zeros(Nx*Ny)

    #Bucle frontera inferior
    for i=1:Nx-1
        N   =   (i-1)*Ny+1
        filas_int  =   vcat(filas_int,N)
        cols_int   =   vcat(cols_int,N)
        vals_int   =   vcat(vals_int,1.0)
        bv[N]   =   1.0
    end

    #Bucle frontera derecha
    for j=1:Ny-1
        N   =   Ny*(Nx-1)+j
        filas_int  =   vcat(filas_int,N)
        cols_int   =   vcat(cols_int,N)
        vals_int   =   vcat(vals_int,1)
        bv[N]   =   1-y[j]
    end

    #Bucle frontera superior
    for i=1:Nx
        N   =   i*Ny
        filas_int  =   vcat(filas_int,N)
        cols_int   =   vcat(cols_int,N)
        vals_int   =   vcat(vals_int,1)
        bv[N]   =   1-x[i]
    end

    #Bucle frontera izquierda

    for j=2:Ny-1
        N   =   j
        filas_int  =   vcat(filas_int,N)
        cols_int   =   vcat(cols_int,N)
        vals_int   =   vcat(vals_int,1)
        bv[N]   =   1
    end

    #Construcción matriz de coeficientes global
    A   =   sparse(filas_int[:],cols_int[:],vals_int[:],Nx*Ny,Nx*Ny)

#     figure()
#     spy(A)

    #Resolución del sistema.
    u_i =   A\bv
    t2  =   time()-t1

    #Representación gráfica de la solución.
    xmat    = zeros(Nx,Ny)
    ymat    = zeros(Nx,Ny)
    umat    = zeros(Nx,Ny)

    for i=1:Nx, j=1:Ny
        N           = (i-1)*Ny+j
        xmat[i,j]   = x[i]
        ymat[i,j]   = y[j]
        umat[i,j]   = u_i[N]
    end

    figure()
    plot3D(xmat, ymat, umat, ".b");
    xlabel("x")
    ylabel("y")

    return A
end

