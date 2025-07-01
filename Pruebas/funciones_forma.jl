using LinearAlgebra
using PyPlot
using Polynomials


function funciones_forma_lagrange_1D(N::Int64)
    #Funciones de forma Lagrange Pk 1D Dominio:[a,b]#
    #N:Orden del polinomio#

    a   =   0
    b   =   1
    xi  =   zeros(Float64,N+1) # Vector de nodos de Chebyshev
    for i=0:N
        xi[i+1]   =   (a+b)/2   -   (b-a)/2*cos(i*pi/N)
    end

    A       =   zeros(Float64,length(xi),length(xi))
    A[:,1]  .=   1
    for i=1:length(xi), j=2:length(xi)
        A[i,j]  =   xi[i]^(j-1)
    end

    Delta   =   I(length(xi))

    C       =   A\Delta

    #Representación

    x = LinRange(a,b,300)

    for j=1:length(xi)
        coeff_ij    =   C[:,j]
        psi_i       =   Polynomial(coeff_ij)
        y           =   psi_i.(x)
        plot(x,y,label="ψ$j(x)")
    end

        # Eje x en y=0
        axhline(0, color="black", linewidth=1)

        # Marcar los nodos en el eje x con puntos rojos
        scatter(xi, zeros(length(xi)), color="red", marker="o", zorder=3, label="Chebyshev nodes")

        xlabel("x")
        ylabel("ψᵢ(x)")
        legend(loc="upper left", bbox_to_anchor=(1, 1), borderaxespad=0.)
        grid(true)

        # Opcional: mover ejes para que X se cruce en y=0
        gca().spines["left"].set_position("zero")
        gca().spines["bottom"].set_position("zero")
        gca().spines["top"].set_visible(false)
        gca().spines["right"].set_visible(false)

        tight_layout()

        figure(figsize=(8, 5))

        for i in 1:length(xi)
            y = zeros(length(x))

            if i == 1
                # Primer nodo: activo entre xi[1] y xi[2]
                idx = findall((x .>= xi[1]) .& (x .<= xi[2]))
                y[idx] .= (xi[2] .- x[idx]) ./ (xi[2] - xi[1])
                elseif i == length(xi)
                # Último nodo: activo entre xi[end-1] y xi[end]
                idx = findall((x .>= xi[end-1]) .& (x .<= xi[end]))
                y[idx] .= (x[idx] .- xi[end-1]) ./ (xi[end] - xi[end-1])
            else
                # Nodo interior: activo entre xi[i-1] y xi[i+1]
                idx1 = findall((x .>= xi[i-1]) .& (x .<= xi[i]))
                y[idx1] .= (x[idx1] .- xi[i-1]) ./ (xi[i] - xi[i-1])

                idx2 = findall((x .>= xi[i]) .& (x .<= xi[i+1]))
                y[idx2] .= (xi[i+1] .- x[idx2]) ./ (xi[i+1] - xi[i])
            end

            plot(x, y, label="ψ$i(x)")
        end

        scatter(xi, zeros(length(xi)), color="red", marker="o", label="Chebyshev nodes")
        axhline(0, color="black", linewidth=1)

        xlabel("x")
        ylabel("ψᵢ(x)")

        legend(loc="upper left", bbox_to_anchor=(1.05, 1))
        grid(true)

        gca().spines["left"].set_position("zero")
        gca().spines["bottom"].set_position("zero")
        gca().spines["top"].set_visible(false)
        gca().spines["right"].set_visible(false)

        tight_layout()


end





