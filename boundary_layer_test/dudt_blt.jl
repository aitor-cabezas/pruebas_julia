

# du/dt boundary layer test



function dudt_blt(tnm1::Float64,ui::Vector{Float64},param::parametros)
    
    epsi        =   param.epsi
    x           =   param.x
    y           =   param.y
    Nx          =   param.Nx
    hx          =   param.hx
    Ny          =   param.Ny
    hy          =   param.hy
    dudt        =   zeros(Nx*Ny)
    
    # Coeficientes separación de variables U(x)*phi(t)
    
    c1          =   1/(exp(-1/epsi)-1) 
    c2          =   -c1/2*(1+exp(-1/epsi))-1/2
    delta_t     =   0.2

    phi         =   1/2*(1+tanh((tnm1-1/2)/delta_t))
    dphidt      =   1/2*((1-(tanh((tnm1-1/2)/delta_t))^2)/delta_t)
    
    for i=2:Nx-1,j=2:Ny-1
        
        N             =   (j-1)*Nx + i
        fx            =   x[i] + c1*exp((x[i]-1)/epsi) + c2
        gy            =   y[j] + c1*exp((y[j]-1)/epsi) + c2
        U             =   fx*gy
        dUdx          =   (1+c1/epsi*exp((x[i]-1)/epsi))*gy
        dUdy          =   (1+c1/epsi*exp((y[j]-1)/epsi))*fx
        d2Udx2        =   (c1/epsi^2*exp((x[i]-1)/epsi))*gy
        d2Udy2        =   (c1/epsi^2*exp((y[j]-1)/epsi))*fx
        Qij           =   U*dphidt +
                          (1+U*phi)*phi*(dUdx + dUdy) -
                          epsi*phi*(d2Udx2 + d2Udy2)
        dudt[N]       =   -ui[N]*((ui[N]-ui[N-1])/hx +
                          (ui[N] - ui[N-Nx])/hy) +
                          epsi*((ui[N+1] - 2*ui[N] + ui[N-1])/hx^2 +
                          (ui[N+Nx] - 2*ui[N] + ui[N-Nx])/hy^2) +
                          Qij
                     
    end
    
    return dudt
 
 
end
    
