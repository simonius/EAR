# inclusion of needed libraries
include("burgparam.jl")



# File contains HOES for cell mean values, 
# used for stencil selection and error estimation

# return Hoes functional for regular grid, centerd
# order 2*p+1 using m = p-1, i.e. 2p+2 points
# S- = {1, ..., p+1}, S+ = {p+2, ..., 2*p+2}
function findHOES(hw)
        dx = (x[2] - x[1]) / Nstd
        m = 2*hw-1
        						# define grid
        lg = collect(1:m+3)*dx
        c = zeros(m+3)					# annihilation coefficients representing the functional
        p = zeros(m+3)					# the right hand side of the system of annihilation equations
        E = zeros(m+3, m+3)				# the system matrix

        p[m+1] = factorial(m)				# set the nonzero components of the right hand side
        p[m+2] = 0.0					
        p[m+3] = dx^(1-m)
        for k=1:m+1
                for l=1:m+3
                        E[k, l] = lg[l]^(k-1)		# build up the system for the annihilation coefficients
                end
        end
        for l = hw+2:2*hw+2
                E[m+2, l] = 1.0
                E[m+3, l] = 2.0*lg[l]
        end
        c = inv(E)*p					# solve for the coefficients
							#  calculation of the scaling coefficient.
        q = 0
        xcent = (lg[hw+1] + lg[hw+2])/2
        for k = hw+2:2*hw+2
             q = q + c[k]*(lg[k] - xcent)
        end

        return c/q					# return the scaled functional
end

# Applies the hoes functional c to the vector of mean values u and returns the result.
function applyhoes(u, c)
        dx = (x[2] - x[1]) / Nstd
        N = length(c)
        w = zeros(N)
        for k = 2:N
                w[k] = w[k-1] + dx*u[k-1]		# calculate the primitve function of u
        end
        return dot(c, w)
end
                
        
