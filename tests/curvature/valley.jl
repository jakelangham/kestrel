# Test 1 of Scheme_curvature_note.tex: steady, frictionless, supercritical flow
# through the valley b = -20 exp(-x^2/40^2) m, with inflow H0 = 0.5 m and
# u0 = 8 m/s. Runs Kestrel at a sequence of resolutions and reports the
# relative L1 errors in u and H against the reference profiles of each system
# in reference.jl, with the observed orders of convergence.
#
# Usage: julia valley.jl <kestrel binary> <work dir> [geometric factors] [dx...]
#   geometric factors: value passed to 'Geometric factors' (default: not set)
#   dx: grid spacings in m (default: 4 2 1 0.5 0.25)

include("reference.jl")
using Printf
using DelimitedFiles

const A_VALLEY = -20.0
const W_VALLEY = 40.0
const H0 = 0.5
const U0 = 8.0
const XTILESIZE = 40.0
const NXTILES = 12           # domain -240 m to 240 m
const TEND = 200.0
# The flow starts from a uniform cap filling the first interior tile at the
# inflow and is then fed by the Dirichlet boundary. (A uniform layer over the
# whole domain traps a hydraulic jump in the valley, because the fluid there
# lacks the energy to climb out.) The domain-edge tiles hold the boundary
# data, so the cap is kept out of them.
const CAP_X = -XTILESIZE * (NXTILES / 2 - 1.5)
const CAP_RADIUS = XTILESIZE / 2 - 1.0
const NOUT = 20

function valley_input(dx, geofactors, outdir)
   nxpertile = round(Int, XTILESIZE / dx)
   isapprox(nxpertile * dx, XTILESIZE) || error("dx must divide $XTILESIZE")
   geo = geofactors === nothing ? "" : "Geometric factors = $geofactors\n"
   return """
%%ListType

Domain:
Lat = 0
Lon = 0
nXtiles = $NXTILES
nYtiles = 1
nXpertile = $nxpertile
nYpertile = 1
Xtilesize = $XTILESIZE
Boundary Conditions = Dirichlet
Boundary Hn = $H0
Boundary U = $U0

Cap:
  capX = $(CAP_X)
  capHeight = $H0
  capRadius = $(CAP_RADIUS)
  capShape = flat
  capU = $U0
  capConc = 0.0

Parameters:
Drag = None
Erosion = Off
Eddy Viscosity = 0
$(geo)
Solver:
T end = $TEND
cfl = 0.5
max dt = 1.0
Tile Buffer = 3
Height threshold = 1e-6

Output:
N out = $NOUT
directory = $outdir/
Format = txt

Topog:
Type = Function
Topog function = xGauss
Topog params = ($A_VALLEY, $W_VALLEY)
"""
end

# Read x, Hn and u from a Kestrel 1D text output file.
function read_output(file)
   data, header = readdlm(file, ',', header=true)
   names = strip.(vec(header))
   col(n) = Float64.(data[:, findfirst(==(n), names)])
   return col("x"), col("flow_depth"), col("x_velocity")
end

function run_valley(prog, workdir, geofactors, dxs)
   mkpath(workdir)
   bd = GaussBed(A_VALLEY, W_VALLEY)
   x0 = -XTILESIZE * NXTILES / 2
   systems = [("C  (1,-1)", 1, -1), ("B' (1, 0)", 1, 0), ("B  (0, 0)", 0, 0), ("A  (1,+1)", 1, 1)]
   errs = Dict(s[1] => Tuple{Float64,Float64}[] for s in systems)
   for dx in dxs
      name = @sprintf("valley_dx%g", dx)
      infile = joinpath(workdir, "Input_$name.txt")
      write(infile, valley_input(dx, geofactors, name))
      cd(workdir) do
         run(pipeline(`$prog Input_$name.txt`, stdout="$name.log", stderr="$name.log"))
      end
      last = joinpath(workdir, name, @sprintf("%06d.txt", NOUT))
      prev = joinpath(workdir, name, @sprintf("%06d.txt", NOUT - 1))
      x, H, u = read_output(last)
      _, Hp, up = read_output(prev)
      steady = max(rel_l1(u, up), rel_l1(H, Hp))
      @printf("dx = %-5g cells = %-5d change over last %g s: %.2e\n", dx, length(x), TEND / NOUT, steady)
      for (sname, k, sigma) in systems
         uref, Href = sname[1] == 'C' ? valley_exact(bd, x; x0=x0, H0=H0, u0=U0) :
                      steady_profile(bd, x, k, sigma; x0=x0, H0=H0, u0=U0)
         push!(errs[sname], (rel_l1(u, uref), rel_l1(H, Href)))
      end
   end
   println("\nRelative L1 errors against each reference (u, H) and observed orders:")
   for (sname, _, _) in systems
      e = errs[sname]
      @printf("%s:", sname)
      for i in eachindex(e)
         @printf("  %.3e, %.3e", e[i]...)
         if i > 1
            ou = log(e[i-1][1] / e[i][1]) / log(dxs[i-1] / dxs[i])
            oH = log(e[i-1][2] / e[i][2]) / log(dxs[i-1] / dxs[i])
            @printf(" (%.2f, %.2f)", ou, oH)
         end
         i < length(e) && print(" |")
      end
      println()
   end
   return errs
end

if abspath(PROGRAM_FILE) == @__FILE__
   prog = abspath(ARGS[1])
   workdir = ARGS[2]
   geofactors = length(ARGS) >= 3 && ARGS[3] != "default" ? ARGS[3] : nothing
   dxs = length(ARGS) >= 4 ? parse.(Float64, ARGS[4:end]) : [4.0, 2.0, 1.0, 0.5, 0.25]
   run_valley(prog, workdir, geofactors, dxs)
end
