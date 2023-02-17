# When included the statements and routines in this file perform a convergence analysis

include("solvers.jl")		# inclusion of the solvers


# File performs a convergence analysis for u(x, 0) = sin(\pi*x)/50 + 1
# The reference solution is calculated by
# solving for an appropiate wave speed u that
# leads back along a characteristic to a position
# on the initial condition with exactly this value
uinit(x) = sin(pi*x) / 50 .+ 1.0
function refsol(x, t, N = 100)		# the sought after value u has to be a fixed point of 
	ures = 1.0			# the map that traces back from (x, t) along the A(u)
	for k = 1:N			# charactersistic to (x_0, 0) and evaluates the initial
                xpred = x - t*ures	# condition at that point.
                ures = uinit(xpred)
        end
        return ures
end

# Evaluates a polynomial given in monomial coefficients at x
# Obviously this is not a good idea for large orders.
function EvalPoly(mcoeff,x)
        val = 0
        for k=1:length(mcoeff)
                val = val + x^(k-1)*mcoeff[k]
        end
        return val
end

# Evaluates the solution at x using a recover with a stencil centerd
# around x
function EvalSol(sol, x)
        N = length(sol)					# the first action of the routine is
        dx = 2/N					# a bisection search for the cell in which x
        starts = (collect(1:N) .- 1.0)*dx		# is located. While overpowered in the case
        lb = 1						# of a equidistant grid, non-equidistant 
        ub = N+1					# grids could be also handled.
        im = 1
        while ub - lb > 1
                im = floor(Int, 0.5*(ub + lb))
                if starts[im] <= x
                        lb = im
                else
                        ub = im
                end
        end
      							# Now a local recovery is done around x
        xloc = (x - starts[lb]-dx/2)/dx			# the local position of x in this grid
        hw = min(w, 5)
        uloc = zeros(2*hw + 1)
        for k=1:2*hw+1
                uloc[k] = sol[mod1(lb + k - hw-1, N)] 	# collect the values from the neighboring cells
        end
        mcoeff = CRM*uloc				# calculate the coefficients using the recovery matrix CRM
        return EvalPoly(mcoeff, xloc)			# evaluate the recovery polynomial
end


# This file performs a convergence analysis for the solver [solver] and saves the resulting plots 
# starting with the filename [basename]. [recon] is the used recovery operator.
function MKConvAna(solver, recon, basename, tend=8.0)
        Narr = floor.(Int, 1.3.^collect(12:22))			# the used number of cells
        Ntest = length(Narr)					# Amount of simulations to carry out.
        rdx = (x[2] - x[1]) / Nref				# distance between two evaluation points 
        Evalps = collect(x[1]:rdx:x[2])[1:Nref] .+ 0.5.*rdx	# the grid of points where the error between reference and simulation is evaluated
        errarr = zeros(Ntest)					# The errors correspond to the cell counts in Narr

        for k=1:Ntest						# This loops iterates over the choosen cell sizes
                println("Running ", k, " of ", Ntest)
                sol = NT(x->sin.(pi*x)/50 .+ 1.0, tend, solver, recon, N=Narr[k])
                println("calculating the norm")
                esum = 0
                for l=1:Nref		# The L1 norm is evaluated summing errors of point evaluation and 
                        esum = esum + abs(refsol(Evalps[l], tend) - EvalSol(sol(tend), Evalps[l]))
                end
                errarr[k] = esum/Nref	# scaling the result accorcing to the number of points
        end

        CairoMakie.activate!()		# enlarge the texts
        fontsize_theme = Theme(fontsize = 28, colormap=:grayC)
        set_theme!(fontsize_theme)
        axis = (xlabel = "Cells", ylabel = "Error", yscale=log10, xscale=log10) # use log axis

					# Plot the errors
        fig = scatter(Narr, errarr, axis = axis, color=:black)
        li = length(Narr)
					# Select a error to hinge the order lines
        mi = li - 5

					# Add order lines with slopes of order 1 to 5 
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1]), errarr[mi]*(Narr[mi]/Narr[end]) ], label = "order 1", color=:black, linestyle=:dot)
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^2, errarr[mi]*(Narr[mi]/Narr[end])^2], label = "order 2", color=:black, linestyle=:dashdot)
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^3, errarr[mi]*(Narr[mi]/Narr[end])^3], label = "order 3", color=:black, linestyle=:dash)
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^4, errarr[mi]*(Narr[mi]/Narr[end])^4], label = "order 4", color=:black)
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^5, errarr[mi]*(Narr[mi]/Narr[end])^5], label = "order 5", color=:black, linestyle=:dot)
        lines!([Narr[1], Narr[end]], [errarr[mi]*(Narr[mi]/Narr[1])^6, errarr[mi]*(Narr[mi]/Narr[end])^6], label = "order 6", color=:black, linestyle=:dashdot)
        axislegend(position=:lb)

	# Save the result
        save(string("pics/convana", w, basename, ".pdf"), fig)
end

# Convergence analysis for Bounding box (Sphere) radius estimator with sphere surface outlier detection and removal
MKConvAna(BurgDIR!, odbbrecon, "outldbbox")

# Convergence analysis for Bounding box (sphere) radius estimator
MKConvAna(BurgDIR!, bbrecon, "bbox")

# Convergence analysis for High Order Edge Sensor based outlier detection with bounding box radius estimate
# MKConvAna(BurgDIR!, hoesbbrecon, "hoesbbox")

# Convergence analysis for High Order Edge Sensor based radius estimate
# MKConvAna(BurgDIR!, hoesrecon, "hoes")

# Convergence analysis for High Order Edge Sensor based outlier detection with variance based radius estimate
# MKConvAna(BurgDIR!, hoesstochrecon, "hoesstoch")

# Convergence analysis for variance based radius estimate
MKConvAna(BurgDIR!, srecon, "stoch")



