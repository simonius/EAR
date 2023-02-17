# This File contains the implementation of the solvers

# Inclusion of needed libraries
using Makie					# For plots
using GLMakie
using CairoMakie
using QuadGK					# Integration of the initial condition
using OrdinaryDiffEq				# Time integration methods


#include("conv.jl")				# convolutions
include("recon.jl")				# recovery operators
include("burgparam.jl")				# global parameters
include("burgflux.jl")				# flux functions
include("derhoes.jl")				# PA operators for derivatives/mean values


# Godunov Scheme for reference solutions
function BurgersGod!(du, u, p, t)
        N = length(u)
        for i = 1:N
                du[i] = (hBurgersGod(u[mod1(i-1, N)], u[i]) - hBurgersGod(u[i], u[mod1(i+1, N)]))/p.dx
        end
end


# Information radius base scheme.
function BurgIR!(du, u, p, t)
        N = length(u)
        dx = p[1]
        dt = p[2]
        rc = p[3]				# the recovery operator
        ur, er = zeros(N), zeros(N)		# space for the recovered values and their radius estimates
        ulocal = zeros(2*w)			# space for a local copy of the conserved variables

	# calculate the recovery and radius estimate
        for k=1:N
                for l = 1:2*w
                        ulocal[l] = u[mod1(k-w+l, N)] # collect the mean values around, respect periodicity
                end
                ur[k], er[k] = rc(ulocal, w) 	# call the recovery procedure
        end

	# use the recovered point values and radius in the scheme.
        for k=1:N
                ulm, umm, urm = u[mod1(k-1, N)], u[mod1(k, N)],u[mod1(k+1, N)]	# left, middle and right mean values
                ulb = min(ur[mod1(k-1, N)] - er[mod1(k-1, N)]) 			# left cell edge lower and upper
                ult = max(ur[mod1(k-1, N)] + er[mod1(k-1, N)]) 			# bounds
                urb = min(ur[mod1(k, N)] - er[mod1(k, N)])   			# right edge lowe and upper bounds
                urt = max(ur[mod1(k, N)] + er[mod1(k, N)])
		# Flux evaluation.
                du[k] = (hBurgDaf(ulb, ult, ulm, umm) - hBurgDaf(urb, urt, umm, urm))/dx
        end
end


UM, SV, VM = svd(Rpoly)
const opnormfac = opnorm(Rpoly, Inf) #maximum(SV)

# Structure of parameters for the solver
mutable struct Params
        dx			# cell size
        dt			# timestep size
        rc			# recovery operator
        c			# used PA operator
end

# Main driver routine for numerical tests.
# u0: 	 	initial condition
# tmax:  	the PDE is solved up to this time
# solver: 	used space discretisation
# rc:		used recovery operator
# N: 		number of points
# CFL:		used grid constant 
function NT(u0, tmax, solver, rc; N=Nstd, CFL=CFLstd)
        tspan = (0, tmax)
        dx = (x[2] - x[1]) / N
        dt = CFL*dx
        xgrid = collect(x[1]:dx:x[2])[1:N] .+ 0.5.*dx

	# Numerical quadrature of the initial condition
        println("integrating initial condition")
        u0ar = zeros(N)
        for k=1:N
                u0ar[k] = quadgk(u0, xgrid[k]-dx/2, xgrid[k]+dx/2, maxevals=10E4)[1] / dx
        end
        println("done")
	# fill farameter struct
        p = Params(dx, dt, rc, findHOES(hoesw)) 
        prob = ODEProblem(solver, u0ar, tspan, p)
	# and solve the defined problem using Hairer Wanner Dorman Price method of order 8
        sol = solve(prob, DP8(), dt=dt, adaptive=false)
        return sol
end

# This function mirrors the driver routine above, but saves time by integrating low order solvers
# like the Godunov solver used for reference solutions using forward Euler steps and big grid constants
function NTEuler(u0, tmax, solver, rc; N=Nstd, CFL=0.5)
        tspan = (0, tmax)
        dx = (x[2] - x[1]) / N
        dt = CFL*dx
        xgrid = collect(x[1]:dx:x[2])[1:N] .+ 0.5.*dx

        # Numerical quadrature of the initial condition
        println("integrating initial condition")
        u0ar = zeros(N)
        for k=1:N
                u0ar[k] = quadgk(u0, xgrid[k]-dx/2, xgrid[k]+dx/2, maxevals=10E4)[1] / dx
        end
        println("done")
        p = Params(dx, dt, rc, findHOES(hoesw))

        prob = ODEProblem(solver, u0ar, tspan, p)

	# up to here as NT, but now using forward Euler steps
        sol = solve(prob, Euler(), dt=dt, adaptive=false)
        return sol
end

# Information radius scheme with redistributed viscosity.
# This redistribution is needed to enforce that the ratio between the jump in the variables and the viscuous
# parts of the flux stays bounded.
function BurgDIR!(du, u, p, t)
        N = length(u)
      	# reconstructed u and error
        ur,er, dur, der = zeros(N), zeros(N), zeros(N), zeros(N)
        # equivalent viscosity
        vis = zeros(N)
        # jumps in the mean values
        jps = zeros(N)
        # buffer for the local values of u
        ulocal = zeros(2*w)
        # maximum possible dissipation for unity radius
        possdis = zeros(N)
	# dissipation that takes place when the reconstructed radii are used directly
        rdis = zeros(N)
        visdiss = zeros(N)

        # discontinuity sensor output
        ds = zeros(N)
        dsloc = zeros(2*hoesw+1)

        # apply discontinuity sensor
        for k=1:N
               for l=0:2*hoesw
                        dsloc[l+1] = u[mod1(k-hoesw+l, N)]
               end
               ds[k] = abs(applyhoes(dsloc, p.c))
        end

        dsloc = zeros(2*w)
        # recosntruct and look at mean jumps
        for k=1:N
                for l = 1:2*w
                        ulocal[l] = u[mod1(k-w+l, N)]
                        dsloc[l] = ds[mod1(k-w+l, N)]
                end
                ur[k], er[k], dur[k], der[k] = p.rc(ulocal, w, dsloc)
                jps[k] = abs(u[mod1(k+1, N)] - u[mod1(k, N)])
        end
	# calculate the equivalent viscosity, i.e. the ratio the radius of information and the jump in the mean values
        vis = er ./ (jps .+ 1.0E-30)

        # calculate actual and possbile dissipation
        for k=1:N
                for l=1:2*w+1
                        bdid = mod1(k+l-w-1, N)
                        rcid = mod1(k+l-w, N)
                        possdis[k] += (jps[bdid])^2*wd[l]
                end
                rdis[k] = jps[k]*er[k]
        end

        # Calculate the needed viscosity
        for k=1:N
                nv = dfac*rdis[k] ./ (possdis[k] + 1.0E-30)
                for l=1:2*w+1
                        bdid = mod1(k+l-w-1, N)
                        visdiss[bdid] += nv*wd[l]
                end
        end

        # calcualte using this viscosity distrubution the 
	# time derivative of the mean values in the cells
        for k=1:N
                du[k] = (f(ur[mod1(k-1, N)]) - f(ur[mod1(k, N)])) / p.dx
                du[k] += visdiss[mod1(k-1, N)] * (u[mod1(k-1, N)] - u[k]) / p.dx
                du[k] += visdiss[k] *(u[mod1(k+1, N)] - u[k]) / p.dx
        end
	# return the viscosity distribution for debugging purposes
        return visdiss
end

