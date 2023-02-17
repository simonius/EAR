# File contains the continuous and numerical flux functions for burgers equation. 
# the numerical fluxes are not in all cases classical two-point fluxes. 
# most of them take 4 arguments and are some of the following
# uml, umr: left and right average values in the neighboring cells
# ulb, uub: a lower and an upper bound on the point value of u at the cell interface
# up: point value
# eps: uncertainty of this point value

# analytic flux and flux jacobian
f(u) = u^2/2
df(u) = u

# Used entropy
U(u) = u^2/2
# Entropy variables/entropy gradient. 
dU(u) = u


# Dafermos flux for Burgers equation. The evaluation is based on the connection to Godunovs flux.
function hBurgDaf(ulb, uub, uml, umr)
        c = dU(umr) - dU(uml)
        if c >= 0
                if ulb*uub <= 0
                        v = 0
                else
                        v = min(abs(ulb), abs(uub))
                end
        else
                v = max(abs(ulb), abs(uub))
        end
        
        return f(v)
end

# the Lax-Friedrichs modified flux from the publication.
function hBurgKLF(ulb, uub, uml, umr)
        um = 0.5*(ulb+uub)
        return f(um) + max(abs(df(ulb)), abs(df(uub)))*abs(ulb-uub)*sign(uml - umr)/2
end 


# Numerical solution to the one-dimensional Riemann problem
# for a scalar conservation law with conves flux
# xdbt is the position in the fan, i.e. x divided by t
# the a and a inv shall be the derivative of the flux and
# its inverse function
function Riemann1D(xdbt, ul, ur, a, ainv, f)
        s = (f(ul) - f(ur))/(ul-ur)
        if ul > ur
                if xdbt < s
                        return ul
                else
                        return ur
                end
        else
                if xdbt <= a(ul)
                        return ul
                elseif xdbt >= a(ur)
                        return ur
                else
                        return ainv(xdbt)
                end
        end
end
# Solution to the Riemann problem for Burgers' equation
function RiemannBurg(xdbt, ul, ur)
        f(u) = u^2/2
        a(u) = u
        ainv(y) = y
        return Riemann1D(xdbt, ul, ur, a, ainv, f)
end

# Godunov flux for Burgers' equation
function hBurgersGod(ul, ur)
        f(u) = u^2/2
        return f(RiemannBurg(0, ul, ur))
end

