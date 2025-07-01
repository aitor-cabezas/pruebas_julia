using LinearAlgebra
using SparseArrays
using PyPlot
using LaTeXStrings

function prueba2(Nx::Int64)
    #Problema convección-difusión diferencias finitas 1D

    t1 = time()

    #Dominio espacial (a,b)
    a=0 #Límite inferior del dominio.
    b=1 #Límite superior del dominio.
    x=LinRange(a,b,Nx)
    Dx=(b-a)/(Nx-1)

    #Parámetros del problema
    c=1.0; #Coeficiente término convectivo
    e=0.01; #Coeficiente término difusión

    #Discretización espacial
    A=zeros(Nx,Nx); #Inicialización de la matriz de coeficientes.
    u_i=zeros(Nx); #Inicialización del vector de soluciones.
    bv=zeros(Nx); #Inicialización del término independinete.

    A[1,1]=1; #Condición de contorno x=0;
    A[Nx,Nx]=1 #Condición de contorno x=L;
    bv[1]=1; #Condición de contorno x=0;
    bv[Nx]=0; #Condición de contorno x=L;
    @inbounds for i=2:(Nx-1)
        A[i,i-1]=-c/Dx-e/Dx^2;
        A[i,i]= c/Dx + 2*e/Dx^2;
        A[i,i+1]= -e/Dx^2;
    end

    #Resolución del sistema.
    u_i=A\bv;
    t2=time()-t1

    #Representación gráfica de la solución.
    if true
        figure()
        plot(x,u_i, "b");
        u_theor = @. (exp(c*b/e)-exp.(c*x/e))/(exp(c*b/e)-1)
        plot(x, u_theor, "k")
    end

    #Cálculo del error
    u_theor = @. (exp(c*b/e)-exp.(c*x/e))/(exp(c*b/e)-1) #Solución analítica
    err=norm(u_i-u_theor,Inf) #Error L2inf

    return Dx, err, t2

end

function prueba2_sparse(Nx::Int64)
    #Problema convección-difusión diferencias finitas 1D

    t1 = time()

    #Dominio espacial (a,b)
    a=0 #Límite inferior del dominio.
    b=1 #Límite superior del dominio.
    x=LinRange(a,b,Nx)
    Dx=(b-a)/(Nx-1)

    #Parámetros del problema
    c=1.0; #Coeficiente término convectivo
    e=0.01; #Coeficiente término difusión

    #Discretización espacial
    A=zeros(Nx,Nx); #Inicialización de la matriz de coeficientes.
    u_i=zeros(Nx); #Inicialización del vector de soluciones.
    bv=zeros(Nx); #Inicialización del término independinete.

    filas   =   vcat([1],
                     2:Nx-1,
                     2:Nx-1,
                     2:Nx-1,
                     [Nx])
    cols    =   vcat([1],
                     1:Nx-2,
                     2:Nx-1,
                     3:Nx,
                     [Nx])
    vals    =   vcat([1],
                     fill(-c/Dx-e/Dx^2,Nx-2),
                     fill(c/Dx+2*e/Dx^2,Nx-2),
                     fill(-e/Dx^2,Nx-2),
                     [1])
    A       =   sparse(filas,cols,vals,Nx,Nx)

    bv[1]   = 1; #Condición de contorno x=0;
    bv[Nx]  = 0; #Condición de contorno x=L;

    #Resolución del sistema.
    u_i=A\bv;
    t2=time()-t1

    #Representación gráfica de la solución.
    if true
        figure()
        u_theor = @. (exp(c*b/e)-exp.(c*x/e))/(exp(c*b/e)-1)
        plot(x, u_theor, "xb")
        plot(x,u_i, "r");
    end

    #Cálculo del error
    u_theor = @. (exp(c*b/e)-exp.(c*x/e))/(exp(c*b/e)-1) #Solución analítica
    err=norm(u_i-u_theor,Inf) #Error L2inf

    return Dx, err, t2

end

function prueba2_convergencia()

    Nxv     = [1000, 2000, 4000, 8000, 16000]
    NNx     = length(Nxv)
    Dxv     = zeros(NNx)
    errv    = zeros(NNx)
    t2v     = zeros(NNx)

    for i=1:NNx
        Dx, err, t2 = prueba2(Nxv[i])
        Dxv[i]      = Dx
        errv[i]     = err
        t2v[i]      = t2
    end

    figure()
    loglog(Dxv,errv)
    xlabel("h")
    ylabel(latexstring("||u_h-u||_{\\infty}"))
    grid("on")

    figure()
    loglog(Dxv,t2v)
    xlabel("Dx")
    ylabel("t")
    grid("on")

end




