# This file plots severall solutions to test the shock-capturing abilities of the schemes.

# inclusion of the solvers
include("solvers.jl")

# routine that performs severall calculations and saves solution plots.
# solver: the used solver
# recon: the used recovery operator
# basename: this name is used together with the order encoded as digit to save the results.
# Ntest, Nref: Number of cells used for the test and the reference
function MKPics(solver, recon, basename, refsolver=BurgersGod!, Ntest = 50, Nref = 5000)
        println("Producing pictures for filename: ", basename)
        CairoMakie.activate!()
        fontsize_theme = Theme(fontsize = 28, colormap=:grayC)		# increase the fontsize
        set_theme!(fontsize_theme)
        
        axis = (xlabel = "x", ylabel = "u")

        # first testcase: sinus-> shock, saved at t=0.3, t=0.6 and t=1.2
        refsol = NTEuler(x->sin.(pi*x) .+ 0.5, 3.0, refsolver, recon, N=Nref)
        sol =  NT(x->sin.(pi*x) .+ 0.5, 3.0, solver, recon, N=Ntest)

        fig = scatter(locgrid(Ntest),sol(0.3), axis = axis, color=:black)
        lines!(locgrid(Nref), refsol(0.3), color=:black)
        save(string("pics/sin", basename, Ntest, "w", w,"t03.pdf"), fig)

        fig = scatter(locgrid(Ntest),sol(0.6), axis = axis, color=:black)
        lines!(locgrid(Nref), refsol(0.6), color=:black)
        save(string("pics/sin", basename, Ntest, "w", w,"t06.pdf"), fig)

        fig = scatter(locgrid(Ntest),sol(1.2), axis = axis, color=:black)
        lines!(locgrid(Nref), refsol(1.2), color=:black)
        save(string("pics/sin", basename, Ntest, "w", w,"t12.pdf"), fig)

        # Total entropy plots to test for the entropy rate
        teref = zeros(Nref)
        te = zeros(Nref)
        tar = collect(range(0.0, 2.0, length=Nref))
        for k=1:Nref
                teref[k] = sum(U.(refsol(tar[k]))) / Nref
                te[k] = sum(U.(sol(tar[k]))) / Ntest
        end

        # second testcase: Is there a sonic point glitch?
        u1(x) = x < 1 ? -1.0 : 1.0 

        refsol = NTEuler(x->u1.(x), 3.0, refsolver, recon, N=Nref)
        sol =  NT(x->u1.(x), 3.0, solver, recon, N=2*Ntest)
        
        fig = scatter(locgrid(2*Ntest),sol(0.2), axis=(limits=(0.5, 1.5, -1.1, 1.1), xlabel="x", ylabel="u"), color=:black)
        lines!(locgrid(Nref), refsol(0.2), color=:black)
        save(string("pics/rarefaction", basename, Ntest, "w", w,"t02.pdf"), fig)
        return tar, teref, te
end



# The following calls perform the tests for several different recovery operators
# Bounding Box radius indicator
tar, teref, te = MKPics(BurgDIR!, bbrecon, "bbox")
efig = lines(tar, teref, label= "Reference", color=:red, linewidth=3.0)
lines!(tar, te, label= "Bounding sphere", color=:black, linestyle=:dot, linewidth=3.0)

# Bounding Box radius indicator with surface outlier removal
tar, teref, te = MKPics(BurgDIR!, odbbrecon, "outldbbox")
lines!(efig.axis, tar, te, label= "Bound. sp. + surf. discard", color = :black, linestyle=:dash, linewidth=3.0)

# High Order Edge Sensor based radius indicator
# tar, teref, te = MKPics(BurgDIR!, hoesrecon, "hoes")
# lines!(efig.axis, tar, te, label= "PA based ind.", color =:black, linestyle=:solid, linewidth=3.0)

# High Order Edge Sensor outlier removal with variance based radius indicator
# tar, teref, te = MKPics(BurgDIR!, hoesstochrecon, "hoesstoch")
# lines!(efig.axis, tar, te, label= "PA w. + Var. ind.", color =:blue, linestyle=:dot, linewidth=3.0)

# High Order Edge Sensor outlier removal with bounding box radius indicator
# tar, teref, te = MKPics(BurgDIR!, hoesbbrecon, "hoesbbox")
# lines!(efig.axis, tar, te, label= "PA w. + Bound. sp. ind.", color =:blue, linestyle=:dash, linewidth=3.0)

# Variance based radius indicator
tar, teref, te = MKPics(BurgDIR!, srecon, "stoch")
lines!(efig.axis, tar, te, label= "Variance ind.", color =:blue, linestyle=:solid, linewidth=3.0)
axislegend(efig.axis)
save(string("pics/totentro","w", w,".pdf"), efig)


