# Test 2 of Scheme_curvature_note.tex: frictionless run-up of a thin layer in
# the bowl b = x^2/(2R), R = 100 m. A parabolic layer (height 0.1 m, radius
# 5 m) is launched from the bottom with speed u0. From the outputs we compute
# the physical mass M = sum(gam Hn dx) (density is constant), the energy
#    E = sum(gam (Hn |u|^2/2 + g Hn b + g Hn^2/(2 gam))) dx   (per unit density)
# with |u|^2 = u^2 gam^2 in 1D, and the mass-weighted height
#    z_b = sum(gam Hn b) / sum(gam Hn).
# For a thin layer, max_t z_b approaches the parcel limit of the note:
#    on (corrected system):       u0^2/(2g)
#    no curvature (model of L23): (R/2)(exp(u0^2/(gR)) - 1)
#    scheme as published:         (R/2) r/(1 - r), r = u0^2/(gR) (unbounded if r >= 1)
#
# Usage: julia bowl.jl <kestrel binary> <work dir> [geometric factors] [r] [dx...]
#   geometric factors: value passed to 'Geometric factors' (default: not set)
#   r = u0^2/(g R) (default 0.81); dx in m (default 1 0.5 0.25)

using Printf
using DelimitedFiles

const R_BOWL = 100.0
const G = 9.81
const H_CAP = 0.1
const L_CAP = 5.0
const XTILESIZE = 20.0
const NXTILES = 100          # domain -1000 m to 1000 m
const TEND = 40.0
const NOUT = 400

parcel_height(r, system) =
   system == :on ? R_BOWL * r / 2 :
   system == :nocurv ? R_BOWL / 2 * (exp(r) - 1) :
   r < 1 ? R_BOWL / 2 * r / (1 - r) : Inf

function bowl_input(dx, r, geofactors, outdir)
   nxpertile = round(Int, XTILESIZE / dx)
   isapprox(nxpertile * dx, XTILESIZE) || error("dx must divide $XTILESIZE")
   u0 = sqrt(r * G * R_BOWL)
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
Boundary Conditions = halt

Cap:
  capX = 0.0
  capHeight = $H_CAP
  capRadius = $L_CAP
  capShape = para
  capU = $u0
  capConc = 0.0

Parameters:
g = $G
Drag = None
Erosion = Off
Eddy Viscosity = 0
$(geo)
Solver:
T end = $TEND
cfl = 0.5
max dt = 0.1
Tile Buffer = 3
Height threshold = 1e-6

Output:
N out = $NOUT
directory = $outdir/
Format = txt

Topog:
Type = Function
Topog function = xParab
Topog params = ($(1 / (2R_BOWL)))
"""
end

# Mass, energy and mass-weighted height from a 1D text output file.
function diagnostics(file, dx)
   data, header = readdlm(file, ',', header=true)
   names = strip.(vec(header))
   col(n) = Float64.(data[:, findfirst(==(n), names)])
   Hn, u, b, bx = col("flow_depth"), col("x_velocity"), col("base_elevation"), col("x_slope")
   gam = sqrt.(1 .+ bx .^ 2)
   M = sum(gam .* Hn) * dx
   E = sum(gam .* (0.5 .* Hn .* u .^ 2 .* gam .^ 2 .+ G .* Hn .* b .+ G .* Hn .^ 2 ./ (2 .* gam))) * dx
   zb = sum(gam .* Hn .* b) / sum(gam .* Hn)
   return M, E, zb, maximum(abs.(col("x")))
end

function run_bowl(prog, workdir, geofactors, r, dxs)
   mkpath(workdir)
   system = geofactors === nothing || lowercase(geofactors) == "on" ? :on :
            lowercase(geofactors) == "no curvature" ? :nocurv : :other
   @printf("u0^2/(gR) = %g, u0 = %.3f m/s; parcel-limit max z_b: on %.1f m, no curvature %.1f m, published scheme %.1f m\n",
           r, sqrt(r * G * R_BOWL), parcel_height(r, :on), parcel_height(r, :nocurv), parcel_height(r, :published))
   results = []
   for dx in dxs
      name = @sprintf("bowl_r%g_dx%g", r, dx)
      write(joinpath(workdir, "Input_$name.txt"), bowl_input(dx, r, geofactors, name))
      ok = true
      cd(workdir) do
         try
            run(pipeline(`$prog Input_$name.txt`, stdout="$name.log", stderr="$name.log"))
         catch
            ok = false
         end
      end
      files = sort(filter(f -> occursin(r"^\d{6}\.txt$", f), readdir(joinpath(workdir, name))))
      M0, E0, _, _ = diagnostics(joinpath(workdir, name, files[1]), dx)
      Mmax = Mmin = M0
      Emax = E0
      zbmax = -Inf
      tzb = 0.0
      for (k, f) in enumerate(files)
         M, E, zb, xext = diagnostics(joinpath(workdir, name, f), dx)
         # Stop once the flow comes within 10 m of the domain boundary.
         xext > XTILESIZE * (NXTILES / 2 - 1) - 10 && break
         Mmax = max(Mmax, M); Mmin = min(Mmin, M)
         Emax = max(Emax, E)
         if zb > zbmax
            zbmax = zb
            tzb = (k - 1) * TEND / NOUT
         end
      end
      bound = E0 / (G * M0)
      @printf("dx = %-5g %s max z_b = %8.3f m (t = %5.2f s)  energy bound E0/(gM) = %7.3f m  max E/E0 = %.4f  M range = %.2e%s\n",
              dx, ok ? "" : "(run stopped early)", zbmax, tzb, bound, Emax / E0, (Mmax - Mmin) / M0,
              system == :other ? "" : @sprintf("  vs parcel limit %+.1f%%", 100 * (zbmax / parcel_height(r, system) - 1)))
      push!(results, (dx=dx, zbmax=zbmax, Emax=Emax / E0, dM=(Mmax - Mmin) / M0, bound=bound))
   end
   return results
end

if abspath(PROGRAM_FILE) == @__FILE__
   prog = abspath(ARGS[1])
   workdir = ARGS[2]
   geofactors = length(ARGS) >= 3 && ARGS[3] != "default" ? ARGS[3] : nothing
   r = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.81
   dxs = length(ARGS) >= 5 ? parse.(Float64, ARGS[5:end]) : [1.0, 0.5, 0.25]
   run_bowl(prog, workdir, geofactors, r, dxs)
end
