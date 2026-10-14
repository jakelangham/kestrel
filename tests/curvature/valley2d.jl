# Two-dimensional variant of Test 1 of Scheme_curvature_note.tex: the valley
# of valley.jl rotated by theta = 30 degrees and shifted, b = -20
# exp(-(s - s0)^2/40^2) m with s = x cos(theta) + y sin(theta), and inflow
# H0 = 0.5 m, speed u0 = 8 m/s in the s direction. The exact solution depends
# on s only, so this exercises the y fluxes, b_xy and the off-diagonal
# geometric terms.
#
# Kestrel's Dirichlet condition imposes the uniform inflow state on every
# edge. The valley is shifted (s0) so that the inflow (west and south) edges
# see only flat bed, where that state is exact, and the valley leaves the
# domain through the north and east edges, which are supercritical outflow.
#
# The simulation starts from the exact solution, written as a Kestrel initial
# condition file. (Filling the valley from a dry or uniform start can strand
# low-energy fluid in the level trough, where it pools; in 2D the spreading
# fronts lose enough energy numerically for this to happen.) The initial
# condition is made by running Kestrel briefly with a Cube covering the domain,
# to obtain an output file with the right tile layout, and replacing its flow
# fields by the exact solution.
#
# Usage: julia valley2d.jl <kestrel binary> <work dir> [geometric factors] [dx...]
#   dx: grid spacings in m (default: 8 4 2)

include("reference.jl")
using Printf
using DelimitedFiles

const THETA = 30.0
const A_VALLEY = -20.0
const W_VALLEY = 40.0
const S0 = 280.0
const H0 = 0.5
const U0 = 8.0
const RHOW = 1000.0          # Kestrel's default water density
const TILESIZE = 40.0
const NTILES = 20            # x and y from -400 m to 400 m
# Errors are measured where |s - s0| <= 150 m, at least 20 m inside the active
# (non-edge) tiles.
const BOX = 340.0
const SBAND = 150.0
const TEND = 100.0
const NOUT = 10

function valley2d_input(dx, geofactors, outdir; template=false, initfile="")
   n = round(Int, TILESIZE / dx)
   isapprox(n * dx, TILESIZE) || error("dx must divide $TILESIZE")
   c, s = cosd(THETA), sind(THETA)
   geo = geofactors === nothing ? "" : "Geometric factors = $geofactors\n"
   half = TILESIZE * NTILES / 2
   source = template ? """
Cube:
  cubeX = 0.0
  cubeY = 0.0
  cubeLength = $(2half - 2)
  cubeWidth = $(2half - 2)
  cubeHeight = $H0
  cubeShape = flat
  cubeU = $(U0 * c)
  cubeV = $(U0 * s)
  cubeConc = 0.0
""" : ""
   init = template ? "" : "Initial condition = $initfile\n"
   return """
%%ListType

Domain:
Lat = 0
Lon = 0
nXtiles = $NTILES
nYtiles = $NTILES
nXpertile = $n
nYpertile = $n
Xtilesize = $TILESIZE
Ytilesize = $TILESIZE
Boundary Conditions = Dirichlet
Boundary Hn = $H0
Boundary U = $(U0 * c)
Boundary V = $(U0 * s)
$(source)
Parameters:
Drag = None
Erosion = Off
Eddy Viscosity = 0
$(geo)
Solver:
T end = $(template ? 1e-3 : TEND)
cfl = 0.25
max dt = 1.0
Tile Buffer = 2
Height threshold = 1e-6
$(init)
Output:
N out = $(template ? 1 : NOUT)
directory = $outdir/
Format = txt

Topog:
Type = Function
Topog function = xyGauss
Topog params = ($A_VALLEY, $W_VALLEY, $THETA, $S0)
"""
end

# Format a number as Fortran's ES18.10E3 does.
function fortran_es(x)
   m, e = split(@sprintf("%.10e", x), 'e')
   ev = parse(Int, e)
   return @sprintf("%13s", m) * "E" * (ev < 0 ? "-" : "+") * lpad(string(abs(ev)), 3, '0')
end

# Write an initial condition file holding the exact solution, on the tile
# layout of the template output file.
function write_exact_ic(template, icfile)
   bd = GaussBed(A_VALLEY, W_VALLEY)
   c, s = cosd(THETA), sind(THETA)
   lines = readlines(template)
   names = strip.(split(lines[1], ','))
   col(n) = findfirst(==(n), names)
   open(icfile, "w") do io
      println(io, lines[1])
      for line in lines[2:end]
         if isempty(strip(line))
            println(io, line)
            continue
         end
         f = split(line, ',')
         v = parse.(Float64, f[2:end])
         val(n) = v[col(n) - 1]
         sw = c * val("x") + s * val("y") - S0
         us, H = valley_exact(bd, [sw]; x0=-400.0, H0=H0, u0=U0)
         us, H = us[1], H[1]
         gam = sqrt(1 + val("x_slope")^2 + val("y_slope")^2)
         b = val("base_elevation") + val("elevation_change")
         new = Dict("flow_depth" => H, "flow_speed" => us,
                    "x_velocity" => us * c, "y_velocity" => us * s,
                    "density" => RHOW, "solids_fraction" => 0.0,
                    "x_flux" => RHOW * H * us * c, "y_flux" => RHOW * H * us * s,
                    "Hnpsi" => 0.0, "w" => b + H / gam)
         for (k, x) in new
            v[col(k) - 1] = x
         end
         println(io, f[1], ", ", join(fortran_es.(v), ", "))
      end
   end
   cp(template * "_topo", icfile * "_topo"; force=true)
end

function read_output2d(file)
   data, header = readdlm(file, ',', header=true)
   names = strip.(vec(header))
   col(n) = Float64.(data[:, findfirst(==(n), names)])
   return col("x"), col("y"), col("flow_depth"), col("x_velocity"), col("y_velocity")
end

function run_valley2d(prog, workdir, geofactors, dxs)
   mkpath(workdir)
   bd = GaussBed(A_VALLEY, W_VALLEY)
   c, s = cosd(THETA), sind(THETA)
   errs = []
   for dx in dxs
      name = @sprintf("valley2d_dx%g", dx)
      tname = name * "_template"
      icbase = joinpath(abspath(workdir), name * "_ic")
      cd(workdir) do
         write("Input_$tname.txt", valley2d_input(dx, geofactors, tname; template=true))
         run(pipeline(`$prog Input_$tname.txt`, stdout="$tname.log", stderr="$tname.log"))
         write_exact_ic(joinpath(tname, "000000.txt"), icbase * ".txt")
         write("Input_$name.txt", valley2d_input(dx, geofactors, name; initfile=icbase))
         run(pipeline(`$prog Input_$name.txt`, stdout="$name.log", stderr="$name.log"))
      end
      x, y, H, u, v = read_output2d(joinpath(workdir, name, @sprintf("%06d.txt", NOUT)))
      _, _, Hp, up, vp = read_output2d(joinpath(workdir, name, @sprintf("%06d.txt", NOUT - 1)))
      inwin = (abs.(x) .<= BOX) .& (abs.(y) .<= BOX) .& (abs.(c .* x .+ s .* y .- S0) .<= SBAND)
      steady = max(rel_l1(u[inwin], up[inwin]), rel_l1(H[inwin], Hp[inwin]))
      sw = c .* x[inwin] .+ s .* y[inwin] .- S0
      us = c .* u[inwin] .+ s .* v[inwin]          # velocity along s
      un = -s .* u[inwin] .+ c .* v[inwin]         # velocity across s (exactly 0)
      uref, Href = valley_exact(bd, sw; x0=-400.0, H0=H0, u0=U0)
      eu = rel_l1(us, uref)
      eH = rel_l1(H[inwin], Href)
      en = sum(abs.(un)) / sum(abs.(us))
      push!(errs, (eu, eH))
      @printf("dx = %-4g cells measured = %-6d change over last %g s: %.2e  rel L1 error vs exact (C): u_s %.3e  H %.3e  |u_n|/|u_s| %.1e",
              dx, sum(inwin), TEND / NOUT, steady, eu, eH, en)
      if length(errs) > 1
         @printf("  orders (%.2f, %.2f)", log(errs[end-1][1] / eu) / log(dxs[length(errs)-1] / dx),
                 log(errs[end-1][2] / eH) / log(dxs[length(errs)-1] / dx))
      end
      println()
   end
   return errs
end

if abspath(PROGRAM_FILE) == @__FILE__
   prog = abspath(ARGS[1])
   workdir = ARGS[2]
   geofactors = length(ARGS) >= 3 && ARGS[3] != "default" ? ARGS[3] : nothing
   dxs = length(ARGS) >= 4 ? parse.(Float64, ARGS[4:end]) : [8.0, 4.0, 2.0]
   run_valley2d(prog, workdir, geofactors, dxs)
end
