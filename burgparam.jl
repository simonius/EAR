# This file collects the global parameters for the test cases

# The number of cells, if not overwritten locally in the testcase
Nstd = 100
# The number of evaluations for the calculation of the error norms in the convergence analyses.
Nref = 2^10
# size of the domain
x = (0.0, 2.0)
# and the time intervall
t = (0.0,1.0)
# CFL number, if not overwritten
CFLstd= 0.1
# order of the used base reconstructions
w = 5
# radius of the used hoes stencil
hoesw = ceil(Int, w/2)

# returns the used cell centers
function locgrid(N)
        dx = (x[2] - x[1]) / N
        grid = collect(x[1]:dx:x[2])[1:N] .+ 0.5.*dx
        return grid
end


using DSP                                                       # holds kernel functions (there named windows) used in the convolutions

fd = max(w, 5)                                                  # we are using a convolution diameter of at least 5 cells.
fw = 2*fd+1                                                     # width of the convolution kernel
const wd = hanning(fw)                                          # the used window function saved as a vector of values
const fs = sum(wd)
const dfac = 1/(maximum(wd/fs))                                 # the needed scaling factor to ensude that the result of a convolution with this
                                                                # kernel of a grid function that vanished everywhere but at one point results in
                                                                # a function that has the same value at that the one non-vanishing grid point.

