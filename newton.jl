
#   Método de Newton 

function newton(fun::Function, J::Function, x0::Vector{Float64}, tol::Float64, iter::Int64)
    
    xkm1    =   x0      #   Semilla
    cont    =   0
   
    while true
       
        cont    += 1
        f        =  fun(xkm1)
        J_eval   =  J(xkm1)  
        xk       =  xkm1 -J_eval\f
        
        if maximum(abs.(xk .- xkm1)) < tol || cont >=iter
            
                   return xk
        end
        
        xkm1    =   copy(xk)
    
    end
   
end





