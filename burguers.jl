using LinearAlgebra
using SparseArrays
using PyPlot
using LaTeXStrings



#   Ecuación de Burguers Diferencias Finitas 2D Método de Newton


function burguers(Nx::Int64,Ny::Int64,Nt::Int64,t0::Int64,tf::Int64)


    t1  =   time()

    #   Dominio Espacial (x1x2 x y1y2)

    x1  =   0   # (0,0)
    x2  =   1   # (1,0)
    x   =   Vector(LinRange(x1,x2,Nx))
    hx  =   (x2-x1)/(Nx-1)  #Paso eje x
    y1  =   0   #(0,1)
    y2  =   1   #(1,1)
    y   =   Vector(LinRange(y1,y2,Ny))
    hy  =   (y2-y1)/(Ny-1)  #Paso eje y

    #   Parámetros del problema

    e   =   0.01        #Coeficiente de Difusión
    
    #   Discretización temporal
    
    u0          =   ones(Nx*Ny)         #Condición Inicial
    ht          =   (tf-t0)/(Nt-1)      #Paso temporal
    u_k         =   zeros(Nt,Nx*Ny)     #Solución en el instante k
    u_k[1,:]    =   u0                  
    
    #   Euler implícito
        

    #   Método de Newton

    tol         =   0.001   #Tolerancia
    iter        =   50
    cont        =   0
    
    for k=2:Nt

        while true

        cont        +=  1
        filas_int   =   zeros(Int64,0,1)
        cols_int    =   zeros(Int64,0,1)
        vals        =   zeros(Float64,0,1)
        f_N         =   zeros(Float64,Nx*Ny)

        #   Bucle nodos interiores

            for i=2:Nx-1,   j=2:Ny-1
            N           =   (i-1)*Ny+j
            filas_int   =   vcat(filas_int,N,N,N,N,N)
            cols_int    =   vcat(cols_int,N,N+1,N-1,N+Ny,N-Ny)
            vals        =   vcat(vals,
                                 1+(2*u_k[k,N]-u_k[k,N-Ny])*ht/hx + (2*u_k[k,N]-u_k[k,N-1])*ht/hy +2*e*ht*(1/hx^2+1/hy^2),
                                 -e*ht/hy^2,
                                 -u_k[k,N]*ht/hy-e*ht/hy^2,
                                 -e*ht/hx^2,
                                 -u_k[k,N]*ht/hx-e*ht/hx^2
                                 )

            f_N[N]      =   u_k[k,N]+
                            u_k[k,N]*ht*(u_k[k,N]-u_k[k,N-Ny])/hx+  u_k[k,N]*ht*(u_k[k,N]-u_k[k,N-1])/hy-
                            e*ht*((u_k[k,N+Ny]-2*u_k[k,N]+u_k[k,N-Ny])/hx^2 + (u_k[k,N+1]-2*u_k[k,N]+u_k[k,N-1])/hy^2)
                            -u_k[k-1,N]
                            
            end


        #   Bucle Frontera Inferior

            for i=0:Nx-1
            N           =   i*Ny+1
            filas_int   =   vcat(filas_int,N)
            cols_int    =   vcat(cols_int,N)
            vals        =   vcat(vals,1.0)
            f_N[N]      =   u_k[k,N] - 1.0
            end

        #   Bucle Frontera Derecha
            for j=2:Ny
            N           =   (Nx-1)*Ny+j
            filas_int   =   vcat(filas_int,N)
            cols_int    =   vcat(cols_int,N)
            vals        =   vcat(vals,1.0)
            f_N[N]      =   u_k[k,N]-(1.0 - y[j])
            end

        #   Bucle Frontera Superior
            for i=1:Nx-1
            N           =   i*Ny
            filas_int   =   vcat(filas_int,N)
            cols_int    =   vcat(cols_int,N)
            vals        =   vcat(vals,1.0)
            f_N[N]      =   u_k[k,N]-(1.0 - x[i])
            end

        #   Bucle Frontera Izquierda
            for j=2:Ny-1
            N           =   j
            filas_int   =   vcat(filas_int,N)
            cols_int    =   vcat(cols_int,N)
            vals        =   vcat(vals,1.0)
            f_N[N]      =   u_k[k,N] - 1.0
            end

        J   =  sparse(filas_int[:],cols_int[:],vals[:],Nx*Ny,Nx*Ny)

        u_knw =   u_k[k,:] .-J\f_N

       if maximum(abs.(u_knw.-u_k[k,:])) < tol || cont >= iter
       println("Convergió en $cont iteraciones.")

       break

       end

       u_k[k,:]    =   u_knw


        end
    
    end

    t2  =   time()-t1

#     #   Representación Gráfica de la Solución
# 
#     xmat    =   zeros(Nx,Ny)
#     ymat    =   zeros(Nx,Ny)
#     umat    =   zeros(Nx,Ny)
#     u_ex    =   zeros(Nx,Ny)
# 
#     for i=1:Nx, j=1:Ny
#         N           =   (i-1)*Ny+j
#         xmat[i,j]   =   x[i]
#         ymat[i,j]   =   y[j]
#         umat[i,j]   =   u_i[N]
#     end
# 
#     j_corte     =   round(Int64,Ny/2)
#     u_num       =   umat[:,j_corte]               #Solución Numérica (x,y(j_corte))
#     L           =   1.0
#     u_ex        =   tanh.((L .- x) ./ (2 * e))    #Solución Exacta 1D
#     



#     figure()
#     plot3D(xmat, ymat, umat, ".b");
#     xlabel("x")
#     ylabel("y")
#     title("Solución Numérica Ecuación de Burguers 2D")
#     figure()
#     plot(x, u_num, label="S.Numérica", linewidth=2)
#     plot(x, u_ex, "--", label="Analí­tica", linewidth=2)
#     xlabel("x")
#     ylabel("u(x)")
#     title("Comparación solución numérica vs analí­tica (en y(j_corte))")
#     legend()
#     grid(true)
    
#     #   Cálculo del error
#     
#     err     =   norm(u_num-u_ex,Inf)
#     
#     return err, hx, hy, t2
#     
# end
# 
# function convergencia()
#     
#     Nxv     =   [120, 140, 180, 200, 220]
#     Nyv     =   Nxv
#     hxv     =   zeros(length(Nxv))
#     hyv     =   zeros(length(Nyv))
#     errv    =   zeros(length(Nxv))
#     t2v     =   zeros(length(Nxv))
#     
#     for i=1:length(Nxv)
#         err, hx, hy, t2     =   burguers(Nxv[i],Nyv[i])
#         errv[i]             =   err
#         hxv[i]              =   hx
#         hyv[i]              =   hy
#         t2v[i]              =   t2
#     end
#     
#     figure()
#     loglog(hxv,errv)
#     xlabel("hx")
#     ylabel(latexstring("|| u_{num}-u_{ex}||_{\\infty}"))
#     grid("on")
#     
#     figure()
#     loglog(hxv,t2v)
#     xlabel("hx")
#     ylabel("t")
#     grid("on")
#         
#     
end

