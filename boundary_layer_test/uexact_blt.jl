function uexact_blt(epsi::Float64,t::Float64,x::Vector{Float64},Nx::Int64,y::Vector{Float64},Ny::Int64)  
    
    uexact          =   zeros(Float64,Nx*Ny)
    c1              =   1/(exp(-1/epsi)-1) 
    c2              =   -c1/2*(1+exp(-1/epsi))-1/2
    delta_t         =   0.2   
    phi             =   1/2*(1+tanh((t-1/2)/delta_t))
    

    
    for i=1:Nx,j=1:Ny
        
        N             =   (j-1)*Nx + i
        fx            =   x[i] + c1*exp((x[i]-1)/epsi) + c2
        gy            =   y[j] + c1*exp((y[j]-1)/epsi) + c2
        uexact[N]     =   1 + phi*fx*gy
    
    end
    
    
    

    return uexact
    
    
end
