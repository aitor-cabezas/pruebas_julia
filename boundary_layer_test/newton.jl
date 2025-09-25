
#   Newton's Method

function newton(fun::Function, J::Function, x0::Vector{Float64}, tol::Float64, iter::Int64)
    
    xkm1    =   x0      #   Seed
    cont    =   0
   
    while true
       
        cont    += 1
        f        =  fun(xkm1)
        J_eval   =  J(xkm1)  
        xk       =  Array(xkm1 -J_eval\f)   #Array() Force the vector to be dense
        
        if maximum(abs.(xk .- xkm1)) < tol || cont >=iter
            
                   return xk
        end
        
        xkm1    =   copy(xk)
    
    end
   
end




