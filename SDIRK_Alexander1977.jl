#   SDIRK

#   Tabla de Butcher
    
    gamma   =   1 - 1/sqrt(2)
    
    A       =   [gamma              0;
                1-2*gamma       gamma]
    b       =   [1/2    1/2]
    c       =   [gamma  1-gamma]
    
    s       =   2   #   Número de etapas
    
    ui      =   zeros(Float64,s)
    
    I       =   Matrix{Float64}(I, Nx*Ny, Nx*Ny)
    

function IRK(tn)
    
    #   Residuo Runge-Kutta
        
        dudt, Jdudt     =   dudt_J_blt(tn + c[1]*ht,ui[1])   
        R1              =   ui[1] - un - ht*A[1,1]*dudt
        JR1             =   I - ht*A[1,1]*Jdudt
        
        dudt, Jdudt     =   dudt_J_blt(tn + c[2]*ht,ui[2])   
        R2              =   ui[2] - un - ht*A[2,2]*dudt
        JR2             =   I - ht*A[2,2]*Jdudt
    
  
    
    
    
    
    
    
    
end
