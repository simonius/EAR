# This file contains the recovery operators

using LinearAlgebra

include("burgparam.jl")

# Performs classic ENO reconstruction of order 2, u is the vector [u_{i-1}, u_i, u_{i+1}]
# returns mu, du/dx
function reconENO2(u)
        a1 = u[2] - u[1]
        a2 = u[3] - u[2]
        if a1^2 < a2^2
                return u[2], a1
        else
                return u[2], a2
        end
end


# Returns the reconstructed values u_l and u_r for cell around x_k, when
# (u_{k-1}, u_k, u_{k+1}) is handed to the function.
function classicENO2(u)
        mean, a = reconENO2(u)
        return mean - a/2, mean + a/2
end

# As classicENO, but takes an entire Matrix of values u, where the first index should be space like
function ENO2Mat(u)
	N, K = size(u)
	ul = zeros(N, K)
	ur = zeros(N, K)
	for i=2:N-1
		for k = 1:K
			ul[i, k], ur[i, k] = classicENO2(u[i-1:i+1, k])
		end
	end
	ul[1, :] = u[1, :]
	ur[1, :] = 0.5*(u[1, :] + u[2, :])
	ul[N, :] = 0.5*(u[N-1, :] + u[N, :])
	ur[N, :] = u[N, :]
	return ul, ur
end


# centered reconstruction using 2*hw+1 points
# yields monomial coefficients
# is only used to calculate the error for smooth solutions
function MonCoeffMat(hw)
        A = zeros(Rational, 2*hw+1, 2*hw+1) # Matrix that yields the mean values of the basis
        
        for k=1:2*hw+1
                for l=1:2*hw+1
                        p = l-1
                        xo = (k-hw - 1//2)
                        xu = (k-hw - 3 //2)
                        A[k, l] = ((xo)^(p+1) - (xu)^(p+1))/(p+1)
                end
        end
        return Float64.(inv(A))
end

CRM = MonCoeffMat(min(w, 5))

# Function to construct linear base reconstructions
# the function uses the monomial basis. While this will be problematic for high orders
# at least up to usual orders the condition numbers of the matrix A is acceptable
# as mean values are not as bad as Vandermonde matrices from point evaluations.
# Note: The entire function assumes dx = 1 as the resulting recovery is scale invariant
# w: the width of the stencil used for reconstruction
# returns a w+1:2*w matrix returning all possible w + 1 constructions possible by 
# selecting a coherent w cell stencil out of 2w cells.
function ReconMats(w)
        A = zeros(Rational, w, w) 	# Matrix that yields the mean values of the basis
        B = zeros(Rational, w+1, w) 	# Matrix that yields the point values at the boundaries
        C = zeros(Rational, w+1, w) 	# Matrix that yields the point values 
                                        # of the first derivative at the boundaries        
        
        for k=1:w
                for l=1:w
                        p = l-1		# the degree of the monomial
                        xo = k - 1//2	# the left cell boundary
                        xu = k - 3 //2 	# the right cell boundary
                        A[k, l] = ((xo)^(p+1) - (xu)^(p+1))/(p+1) # the corresponding mean value
                end
        end
 
        for k=0:w
                for l=1:w
                        p = l-1			# degree
                        x = k-0.5		# left edge of the cell
                        B[k+1, l] = x^p		# evaluation 
                        C[k+1, l] = p*x^(p-1)	# evaluation of the derivative
                end
        end
	# M has all reconstructions from a 2*w+1 stencil
	# Taking the inverse of A results in an operator returning the monomial basis
	# representation of the given mean values while B and C evaluate this
	# representation at the respective left cell edges
        M = Float64.(B*inv(A))
        N = Float64.(C*inv(A))
	# We convert this to all reconstructions around the middle cell interface 
	# with with stencil -w:w by shifting the respective functionals
        R = zeros(w+1, 2*w)
        Q = zeros(w+1, 2*w)
        for k=1:w+1
                R[k, w+2-k:2*w+1-k] = M[k, 1:w]
                Q[k, w+2-k:2*w+1-k] = N[k, 1:w]
        end

        return R, Q
end


const Rpoly, Qpoly = ReconMats(w)

# Function performs bounding box/bounding sphere error estimation
# u: 	2w stencil
# w: 	width
# ds:	hoes output (unused)
# returns recovery, radius estimate, derivative recovery, and derivative recovery radius estimate
function bbrecon(u, w, ds)
        vals  = Rpoly*u			# reconstruct the values
        dvals = Qpoly*u
        maxv = maximum(vals)		# determine their maximum and minimum
        minv = minimum(vals)
        maxdv = maximum(dvals)
        mindv = minimum(dvals)
	# return the center of the maximum and minimum and the radius
        return (maxv + minv)/2, (maxv-minv)/2, (maxdv + mindv)/2, (maxdv - mindv)/2
end

# Function performs bounding box error/bounding sphere estimation
# The value lying on the surface of the bounding sphere is discarded 
# two times to allow 
# a some stencil selection.
# u:	2*w stencil
# w:	width
# ds:	hoes output(unused)
function odbbrecon(u, w, ds)
        vals  = Rpoly*u
        dvals = Qpoly*u
        dists = zeros(length(vals))
        for k=1:2			# loop discards two recoveries that lie on the surface
                maxv = maximum(vals)
                minv = minimum(vals)
                maxdv = maximum(dvals)
                mindv = minimum(dvals)
                cent, err, dcent, derr = (maxv + minv)/2, (maxv-minv)/2, (maxdv + mindv)/2, (maxdv - mindv)/2
                if k == 1
                        dists[:] = abs.(cent .- vals)
                end
                outlier = argmax(dists)	# a recovery with maximal distance has to lie on the surface
                vals = vals[1:end .!= outlier]
                dists = dists[1:end .!= outlier]
        end
        maxv = maximum(vals)
        minv = minimum(vals)
        return (maxv + minv)/2, (maxv - minv) /2, 0.0, 0.0
end



# Function performs recovery with a variance based radius estimate.
function srecon(u, w, ds)
        vals  = Rpoly*u
        dvals = Qpoly*u
        mean = sum(vals) / (w+1)
        dmean = sum(dvals) / (w+1)

        var = 0.0*u[1]
        dvar = 0.0*u[1]
        for k=1:w+1			# variance calculation
                var = var + (vals[k] - mean)^2/(w)
                dvar = dvar + (dvals[k] - dmean)^2/(w)
        end
        tvar = var/(w)
        tdvar = dvar/(w)
        terr = sqrt(tvar)		# radius estimate is the square root of the total variance
        tderr = sqrt(tdvar)
        return mean, terr, dmean, tderr
end



#### HOES based algorithms

# Params: u mean values, w width, ds discontinuity sensor output on the stencil. 
# basic hoes algorithm. Not reported in the publication
function hoesrecon(u, w, ds)
        vals  = Rpoly*u
        N = length(vals)
        badness = abs.(Rpoly)*ds
        best = argmin(badness)
        return vals[best], badness[best], 0.0, 0.0 
end

# the best of both worlds - combines the stochastic error estimate with 
# the hoes weighting of the mean
# Not reported in the publication
function hoesstochrecon(u, w, ds)
        vals  = Rpoly*u
        N = length(vals)
        badness = abs.(Rpoly)*ds
        ow = ones(N)
        eps  = 10E-30
      
        wg = ow ./ (eps .+ badness.^2)
        wg = wg ./ sum(wg)
        mean = sum(vals.*wg)
        var = 0.0
        for k=1:N
                var = var + (vals[k] - mean)^2*wg[k]^2
        end
        tvar = var/sum(wg.^2)
        terr = sqrt(tvar)
        return mean, terr, 0.0, 0.0

end

# the best of both worlds - combines the bounding box error estimate with the 
# hoes outlier detection and discard
# not reported.
function hoesbbrecon(u, w, ds)
        vals  = Rpoly*u
        N = length(vals)
        badness = abs.(Rpoly)*ds
        for k=1:2
                outlier = argmax(badness)
                vals = vals[1:end .!= outlier]
                badness = badness[1:end .!= outlier]
        end
        maxv = maximum(vals)
        minv = minimum(vals)
        return (maxv + minv)/2, (maxv - minv) /2, 0.0, 0.0 
end
