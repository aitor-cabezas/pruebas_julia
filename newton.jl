

using LinearAlgebra
using ForwardDiff


#   Método de Newton 

function Newton(fun::Function,x0::Vector,tol::Float64,iter::Int64)
    
xkm1    =   x0
cont    =   0
   
   while true
       
        cont    += 1
        f,J     =  fun(xkm1)
        xk      =   xkm1 -J\f
        
        if maximum(abs.(xk .- xkm1)) < tol || cont >=iter
            
                   return xk
                   
        end
        
        xkm1    =   xk
    
    end
   
end





