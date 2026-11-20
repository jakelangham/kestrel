# Library for the benchmarks and other non-critical checks: tests that measure
# accuracy or behaviour without a firm pass/fail criterion, and so are kept out
# of the test suite in ../tests (run by `make check` and CI), which checks only
# the essential properties of the code.
#
# It includes ../tests/testlib.jl for running Kestrel and reading its results,
# and adds exact solutions (Thacker's oscillations in a bowl, and dam breaks),
# with writers of exact initial conditions, error measures and profiles for
# plotting. make_exact_ic, exact_errors and exact_profile dispatch on the kind
# of solution.

include(joinpath(@__DIR__, "..", "tests", "testlib.jl"))

# Exact solutions of Thacker (1981) for frictionless flow in the bowl
# b = h0 (x^2 + y^2) / a^2 (in 1D, b = h0 x^2 / a^2), with g = 9.81. sol is a
# named tuple with fields kind, h0 and a, and:
#   kind = :planar, d -- the free surface is a plane and the water is a
#      paraboloid (in 1D, a parabola) of radius a. Its centre moves round the
#      circle of radius d, starting at (d, 0) (in 1D, it oscillates along x),
#      with the uniform velocity (-d w sin(w t), d w cos(w t)), w = sqrt(2 g h0) / a.
#   kind = :radial, A -- (2D only) a radially symmetric oscillation that starts
#      at rest with a curved free surface, w = sqrt(8 g h0) / a, 0 < A < 1.
# The functions below give the frequency w, and the free surface elevation and
# the velocity at (x, y, t). The free surface is that of the water where wet,
# and its analytic continuation elsewhere, so the depth is max(surface - b, 0)
# for any bed b. The Input file's topography must be xparab (1D) or xyparab
# (2D) with coefficients h0 / a^2.
const thacker_g = 9.81
function thacker_frequency(sol, dim)
   if sol.kind == :planar
      return sqrt(2 * thacker_g * sol.h0) / sol.a
   elseif sol.kind == :radial && dim == 2
      return sqrt(8 * thacker_g * sol.h0) / sol.a
   end
   error("Unknown Thacker solution $sol in $(dim)D")
end

function thacker_surface(sol, dim, x, y, t)
   h0, a = sol.h0, sol.a
   ωt = thacker_frequency(sol, dim) * t
   if sol.kind == :planar
      xc = sol.d * cos(ωt)
      yc = (dim == 1) ? 0.0 : sol.d * sin(ωt)
      return h0 + (h0 / a^2) * (2 * x * xc - xc^2 + 2 * y * yc - yc^2)
   else
      A = sol.A
      D = 1 - A * cos(ωt)
      r2 = (x^2 + y^2) / a^2
      return h0 * (sqrt(1 - A^2) / D - r2 * (1 - A^2) / D^2 + r2)
   end
end

function thacker_velocity(sol, dim, x, y, t)
   ω = thacker_frequency(sol, dim)
   if sol.kind == :planar
      return (-sol.d * ω * sin(ω * t), (dim == 1) ? 0.0 : sol.d * ω * cos(ω * t))
   else
      c = ω * sol.A * sin(ω * t) / (2 * (1 - sol.A * cos(ω * t)))
      return (c * x, c * y)
   end
end

# Largest distance of the shoreline from the centre of the bowl over time.
function thacker_max_radius(sol)
   if sol.kind == :planar
      return sol.a + sol.d
   else
      return sol.a * ((1 + sol.A) / (1 - sol.A))^0.25
   end
end

# Cell averages of the exact depth and discharge (depth times velocity) over
# the cell centred at (x, y) of size dx by dy, on the bed that Kestrel uses: the
# exact bed at the cell vertices, interpolated linearly (bilinearly in 2D).
function thacker_cell_average(sol, dim, x, y, dx, dy, t)
   bed(x, y) = sol.h0 * (x^2 + y^2) / sol.a^2
   h = 0.0; qx = 0.0; qy = 0.0
   if dim == 1
      n = 16
      bl, br = bed(x - dx / 2, 0.0), bed(x + dx / 2, 0.0)
      for i = 1:n
         s = (i - 0.5) / n
         xs = x + (s - 0.5) * dx
         hs = max(thacker_surface(sol, dim, xs, 0.0, t) - (bl + s * (br - bl)), 0.0)
         u, _ = thacker_velocity(sol, dim, xs, 0.0, t)
         h += hs; qx += hs * u
      end
      return h / n, qx / n, 0.0
   end
   n = 8
   b00, b10 = bed(x - dx / 2, y - dy / 2), bed(x + dx / 2, y - dy / 2)
   b01, b11 = bed(x - dx / 2, y + dy / 2), bed(x + dx / 2, y + dy / 2)
   for j = 1:n, i = 1:n
      s = (i - 0.5) / n; r = (j - 0.5) / n
      xs = x + (s - 0.5) * dx; ys = y + (r - 0.5) * dy
      b = b00 * (1 - s) * (1 - r) + b10 * s * (1 - r) + b01 * (1 - s) * r + b11 * s * r
      hs = max(thacker_surface(sol, dim, xs, ys, t) - b, 0.0)
      u, v = thacker_velocity(sol, dim, xs, ys, t)
      h += hs; qx += hs * u; qy += hs * v
   end
   return h / n^2, qx / n^2, qy / n^2
end

# The value of setting key in the lines of an input file, or nothing.
function input_setting(lines, key)
   for line in strip_comments(lines)
      m = match(Regex("^\\s*$key\\s*=\\s*(.*?)\\s*\$", "i"), line)
      if m !== nothing
         return m.captures[1]
      end
   end
   return nothing
end

# Cell spacings of a result file's cells, from the distinct x (and y) values.
function cell_spacing(v)
   u = sort(unique(v))
   return length(u) > 1 ? minimum(diff(u)) : 0.0
end

# Write an exact initial condition for the test testname to
# $(testname)_ic/000000.txt, which Input_$testname.txt must load with
# 'T start = 0' and 'Initial condition = $(testname)_ic/000000'. The layout of
# the result file (tiles, cells and topography) is taken from the t = 0 output
# of a run of the same input with the initial condition file replaced by the
# block tmplblock (a Cap or Cube covering the flow's extent, so that the tiles
# it covers are active). exact(x, y, dx, dy) gives the cell averages of the
# bed-normal depth Hn and of the discharges Hn u and Hn v at t = 0. The free
# surface is w = b + Hn / gam (Kestrel's ComputeHn inverted), with the
# geometric correction factor gam of the input file's setting. Returns true on
# error.
function write_exact_ic(testname, dim, tmplblock, exact)
   infile = test_inputfile(testname)
   lines = readlines(infile)
   icdir = "$(testname)_ic"
   tmpl = "$(testname)_template"
   gf = input_setting(lines, "Geometric factors")
   geometric = (gf === nothing) || lowercase(gf) != "off"

   tmplfile = test_inputfile(tmpl)
   open(tmplfile, "w") do io
      for line in lines
         if occursin(r"^\s*(T start|Initial condition)\s*="i, line)
            continue
         elseif occursin(r"^\s*directory\s*="i, line)
            println(io, "directory = $tmpl/")
         elseif occursin(r"^\s*T end\s*="i, line)
            println(io, "T end = 1e-9")
         elseif occursin(r"^\s*N out\s*="i, line)
            println(io, "N out = 1")
         else
            println(io, line)
         end
      end
      println(io, "\n", tmplblock)
   end
   rm(tmpl, force = true, recursive = true)
   failed = run_simulation(tmpl)
   rm(tmplfile)
   if failed
      return true
   end

   # Replace the flow in the template's t = 0 output by the exact solution.
   src = joinpath(tmpl, "000000.txt")
   col, _ = read_result_txt(src)
   dx = cell_spacing(col("x"))
   dy = (dim == 1) ? 0.0 : cell_spacing(col("y"))
   rawlines = readlines(src)
   names = strip.(split(rawlines[1], ','))
   idx(name) = findfirst(==(name), names)
   mkpath(icdir)
   open(joinpath(icdir, "000000.txt"), "w") do io
      println(io, rawlines[1])
      for line in rawlines[2:end]
         if isempty(strip(line))
            println(io, line)
            continue
         end
         f = strip.(split(line, ','))
         v = parse.(Float64, f)
         x = v[idx("x")]
         y = (dim == 1) ? 0.0 : v[idx("y")]
         h, qx, qy = exact(x, y, dx, dy)
         u = (h > 0) ? qx / h : 0.0
         vv = (h > 0) ? qy / h : 0.0
         rho = v[idx("density")]
         sx = v[idx("x_slope")]
         sy = (dim == 1) ? 0.0 : v[idx("y_slope")]
         gam = geometric ? sqrt(1 + sx^2 + sy^2) : 1.0
         v[idx("flow_depth")] = h
         v[idx("w")] = v[idx("base_elevation")] + v[idx("elevation_change")] + h / gam
         v[idx("x_velocity")] = u
         v[idx("x_flux")] = rho * qx
         if dim == 2
            v[idx("y_velocity")] = vv
            v[idx("y_flux")] = rho * qy
         end
         v[idx("flow_speed")] = sqrt(u^2 + vv^2)
         print(io, @sprintf("%8d", Int(v[1])))
         for val in v[2:end]
            print(io, ", ", @sprintf("%18.10E", val))
         end
         println(io)
      end
   end
   cp(src * "_topo", joinpath(icdir, "000000.txt_topo"), force = true)
   rm(tmpl, recursive = true)

   return false
end

# Errors of the (.txt) results in directory testname against an exact solution
# at their output times, which are taken as k * (T end) / (N out) for the k-th
# output (T start must be 0). exact(x, y, dx, dy, t) gives the cell averages of
# the exact Hn, Hn u and Hn v; only cells for which inwindow(x, dx, t) is true
# are compared (where the exact solution holds). Errors are measured in the
# 1-norm over those cells, relative to the exact solution:
#   depth: sum |h - h_exact| / sum h_exact,
#   discharge: sum |q - q_exact| / sum |h_exact| U, where q is the discharge
#      vector and U a velocity scale of the solution.
# Returns the largest depth and discharge errors over the output times.
function window_errors(testname, dim, exact, inwindow, U)
   lines = readlines(test_inputfile(testname))
   tend = parse(Float64, input_setting(lines, "T end"))
   nout = parse(Int, input_setting(lines, "N out"))

   ehmax = 0.0
   eqmax = 0.0
   for f in readdir(testname)
      if !is_resultfile_txt(f)
         continue
      end
      k = parse(Int, f[1:end-4])
      t = k * tend / nout
      col, _ = read_result_txt(joinpath(testname, f))
      xs = col("x")
      ys = (dim == 1) ? zeros(length(xs)) : col("y")
      dx = cell_spacing(xs)
      dy = (dim == 1) ? 0.0 : cell_spacing(ys)
      h = col("flow_depth")
      rho = col("density")
      qx = col("x_flux") ./ rho
      qy = (dim == 1) ? zeros(length(xs)) : col("y_flux") ./ rho
      eh = 0.0; eq = 0.0; hsum = 0.0
      for i in eachindex(xs)
         inwindow(xs[i], dx, t) || continue
         he, qxe, qye = exact(xs[i], ys[i], dx, dy, t)
         eh += abs(h[i] - he)
         eq += hypot(qx[i] - qxe, qy[i] - qye)
         hsum += he
      end
      ehmax = max(ehmax, eh / hsum)
      eqmax = max(eqmax, eq / (hsum * U))
   end

   return ehmax, eqmax
end

# Write the last (.txt) output in directory testname along the x axis (in 2D,
# the row of cells nearest y = 0), with the exact solution, to file, for
# plotting. exact and inwindow are as for window_errors: only the cells where
# errors are measured are written, with the exact cell averages. The file is
# comma separated, like Kestrel's output, with a row of column names: x,
# flow_depth, x_discharge, y_discharge, exact_flow_depth, exact_x_discharge and
# exact_y_discharge (discharges are those of window_errors; the y components
# are 0 in 1D). Returns the output time and the y of the row.
function write_window_profile(file, testname, dim, exact, inwindow)
   lines = readlines(test_inputfile(testname))
   tend = parse(Float64, input_setting(lines, "T end"))
   nout = parse(Int, input_setting(lines, "N out"))
   f = maximum(filter(is_resultfile_txt, readdir(testname)))
   t = parse(Int, f[1:end-4]) * tend / nout
   col, _ = read_result_txt(joinpath(testname, f))
   xs = col("x")
   ys = (dim == 1) ? zeros(length(xs)) : col("y")
   dx = cell_spacing(xs)
   dy = (dim == 1) ? 0.0 : cell_spacing(ys)
   h = col("flow_depth")
   rho = col("density")
   qx = col("x_flux") ./ rho
   qy = (dim == 1) ? zeros(length(xs)) : col("y_flux") ./ rho
   y0 = ys[argmin(abs.(ys))]
   cells = sort(filter(i -> ys[i] == y0 && inwindow(xs[i], dx, t), eachindex(xs)), by = i -> xs[i])
   open(file, "w") do io
      println(io, join(["x", "flow_depth", "x_discharge", "y_discharge",
                        "exact_flow_depth", "exact_x_discharge", "exact_y_discharge"], ", "))
      for i in cells
         row = (xs[i], h[i], qx[i], qy[i], exact(xs[i], ys[i], dx, dy, t)...)
         println(io, join([@sprintf("%.10E", v) for v in row], ", "))
      end
   end
   return t, y0
end

# The exact initial condition, errors and profiles (see write_window_profile)
# for any of the solutions here. exact_profile returns the output time and the
# y of the row written, or nothing if the input does not match the solution.
make_exact_ic(testname, dim, sol) =
   (sol.kind == :dambreak) ? make_dambreak_ic(testname, dim, sol) : make_thacker_ic(testname, dim, sol)
exact_comparison(testname, dim, sol) =
   (sol.kind == :dambreak) ? dambreak_comparison(testname, sol) : thacker_comparison(sol, dim)
function exact_errors(testname, dim, sol)
   c = exact_comparison(testname, dim, sol)
   return (c === nothing) ? (NaN, NaN) : window_errors(testname, dim, c...)
end
function exact_profile(file, testname, dim, sol)
   c = exact_comparison(testname, dim, sol)
   return (c === nothing) ? nothing : write_window_profile(file, testname, dim, c[1], c[2])
end

# Write the exact initial condition of the Thacker solution sol (see
# write_exact_ic). Returns true on error.
function make_thacker_ic(testname, dim, sol)
   infile = test_inputfile(testname)
   lines = readlines(infile)

   # Check the bowl in the input file against the solution.
   tp = input_setting(lines, "Topog params")
   coeffs = parse.(Float64, split(strip(tp, ['(', ')', ' ']), ','))
   if any(abs.(coeffs[1:dim] .- sol.h0 / sol.a^2) .> 1e-12 * sol.h0 / sol.a^2)
      println("Topog params in $infile do not match the solution's bowl.")
      return true
   end
   g = input_setting(lines, "g")
   if g !== nothing && parse(Float64, g) != thacker_g
      println("g in $infile must be $thacker_g.")
      return true
   end

   # Template: a level lake over the bowl, out to beyond the furthest
   # shoreline, so that the tiles it will reach are all active.
   R = thacker_max_radius(sol) * 1.05
   level = sol.h0 * (R / sol.a)^2
   tmplblock = (dim == 1) ?
      "Cap:\ncapX = 0.0\ncapHeight = $level\ncapRadius = $R\ncapShape = level\ncapU = 0\ncapConc = 0.0" :
      "Cube:\ncubeX = 0.0\ncubeY = 0.0\ncubeHeight = $level\ncubeLength = $(2R)\ncubeWidth = $(2R)\n" *
      "cubeShape = level\ncubeU = 0\ncubeV = 0\ncubeConc = 0.0"
   return write_exact_ic(testname, dim, tmplblock,
                         (x, y, dx, dy) -> thacker_cell_average(sol, dim, x, y, dx, dy, 0.0))
end

# The comparison with the Thacker solution sol: exact, inwindow and U for
# window_errors, compared everywhere, with U the largest exact speed at the
# shoreline over a period.
function thacker_comparison(sol, dim)
   ω = thacker_frequency(sol, dim)
   R = thacker_max_radius(sol)
   U = maximum(t -> hypot(thacker_velocity(sol, dim, R, 0.0, t)...),
               range(0, 2π / ω, length = 101))
   return (x, y, dx, dy, t) -> thacker_cell_average(sol, dim, x, y, dx, dy, t), (x, dx, t) -> true, U
end

# Dam breaks in 1D: exact solutions of Kestrel's equations for a reservoir of
# depth hl on -L < x < 0, released at t = 0, on a uniform slope with Coulomb
# friction. sol is a named tuple with fields kind = :dambreak and:
#   hl, L -- the reservoir's bed-normal depth and length;
#   u0 -- the reservoir's initial (horizontal) velocity, downhill; > 0 only
#      with a dry bed (hr = 0). Released from rest (u0 = 0), the direction of
#      Coulomb friction is undefined at first, and Kestrel's implicit drag,
#      with a coefficient proportional to 1 / |u|, makes an O(dt) error that
#      persists (the friction doesn't depend on speed, so nothing restores
#      it). Starting the flow in motion avoids that;
#   hr, R -- the bed-normal depth and length of a layer at rest on 0 < x < R
#      ahead of it (hr = 0: a dry bed);
#   s -- the slope, downhill towards +x: the bed is b = -s x;
#   mu -- the Coulomb friction coefficient (0: frictionless);
#   margin -- the distance kept from the disturbances coming from the far
#      ends of the reservoir and the layer when comparing (see below).
# On a uniform slope, Kestrel's 1D equations, in the horizontal coordinate x
# and velocity u, with h the bed-normal depth Hn and gam = sqrt(1 + s^2) (1 if
# geometric factors are off), are
#    h_t + (h u)_x = 0,
#    (h u)_t + (h u^2 + G h^2 / 2)_x = h A,   G = g / gam^3,
#    A = g (s - mu sign(u)) / gam^2.
# (In HydraulicRHS.f90 and Equations.f90: w - b = h / gam, the flux of w is
# divided by gam^2 and the hydrostatic flux g (w - b)^2 / 2 by gam, the
# gravity source is g h s / gam^2, and the Coulomb friction mu (g / gam) h is
# divided by the slope-aligned speed gam |u|.) With cos(theta) = 1 / gam and
# sin(theta) = s / gam, these are the Savage-Hutter equations for flow down an
# incline, with pressure g cos(theta) h^2 / 2 and forcing
# g h (sin(theta) - mu cos(theta)), written in the horizontal coordinate. They
# are the shallow water equations with gravity G and, where u > 0, the uniform
# acceleration A. In the frame xi = x - u0 t - A t^2 / 2, v = u - u0 - A t,
# A and u0 disappear (as in Mangeney, Heinrich & Roche 2000, Pure Appl.
# Geophys. 157, 1081), leaving the dam break of Ritter (1892) for a dry bed,
# and that of Stoker (1957) for a wet one: a rarefaction into the reservoir,
# and a shock into the layer. With mu > 0 this needs u > 0 wherever there is
# flow, which holds if A > 0 (the slope is steeper than the friction angle,
# s > mu), u0 >= 0 and hr = 0 (the wet-bed case is frictionless here). The solution holds until disturbances from the
# far ends of the reservoir (at -L) and the layer (at R) arrive, which travel
# into the undisturbed states at the speed c = sqrt(G h) in the moving frame.
# The errors are measured only where they have not yet arrived, less margin:
# xi > -L + cl t + margin and, if hr > 0, xi < R - cr t - margin.
const dambreak_g = 9.81

function dambreak_constants(sol; geometric = true)
   gam = geometric ? sqrt(1 + sol.s^2) : 1.0
   G = dambreak_g / gam^3
   A = dambreak_g * (sol.s - sol.mu) / gam^2
   cl = sqrt(G * sol.hl)
   cr = sqrt(G * sol.hr)
   hm = um = cm = S = 0.0
   if sol.hr > 0
      # Middle state: the rarefaction's u = 2 (cl - cm) equals the speed behind
      # a shock into the layer at rest.
      f(h) = 2 * (cl - sqrt(G * h)) - (h - sol.hr) * sqrt(G * (h + sol.hr) / (2 * h * sol.hr))
      a, b = sol.hr, sol.hl
      for _ = 1:200
         m = (a + b) / 2
         (f(m) > 0) ? (a = m) : (b = m)
      end
      hm = (a + b) / 2
      cm = sqrt(G * hm)
      um = 2 * (cl - cm)
      S = hm * um / (hm - sol.hr)
   end
   return (; gam, G, A, cl, cr, hm, um, cm, S)
end

# Exact bed-normal depth and velocity at (x, t).
dambreak_shift(sol, k, t) = sol.u0 * t + k.A * t^2 / 2
function dambreak_state(sol, k, x, t)
   if t == 0
      return (x < 0) ? (sol.hl, sol.u0) : (x < sol.R ? sol.hr : 0.0, 0.0)
   end
   ξ = x - dambreak_shift(sol, k, t)
   r = ξ / t
   if r <= -k.cl
      h, v = sol.hl, 0.0
   elseif sol.hr == 0
      h, v = (r < 2 * k.cl) ? ((2 * k.cl - r)^2 / (9 * k.G), 2 * (k.cl + r) / 3) : (0.0, 0.0)
   elseif r <= k.um - k.cm
      h, v = (2 * k.cl - r)^2 / (9 * k.G), 2 * (k.cl + r) / 3
   elseif r <= k.S
      h, v = k.hm, k.um
   else
      h, v = sol.hr, 0.0
   end
   return h, (h > 0) ? v + sol.u0 + k.A * t : 0.0
end

# Points where the exact solution or its derivative is discontinuous.
function dambreak_breaks(sol, k, t)
   t == 0 && return [0.0]
   r = (sol.hr == 0) ? [-k.cl, 2 * k.cl] : [-k.cl, k.um - k.cm, k.S]
   return r .* t .+ dambreak_shift(sol, k, t)
end

# Cell averages of the exact h and h u over the cell centred at x of width dx,
# by 4-point Gauss-Legendre quadrature between the discontinuities.
const gauss4 = ((-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
                (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538))
function dambreak_cell_average(sol, k, x, dx, t)
   a, b = x - dx / 2, x + dx / 2
   pts = sort([a; filter(p -> a < p < b, dambreak_breaks(sol, k, t)); b])
   H = 0.0; Q = 0.0
   for j = 1:length(pts)-1
      m, w = (pts[j] + pts[j+1]) / 2, (pts[j+1] - pts[j]) / 2
      for (z, wt) in gauss4
         h, u = dambreak_state(sol, k, m + w * z, t)
         H += wt * w * h
         Q += wt * w * h * u
      end
   end
   return H / dx, Q / dx, 0.0
end

# Check the input file of a dam break against sol, and return its constants.
function dambreak_check(testname, sol)
   infile = test_inputfile(testname)
   lines = readlines(infile)
   bad(msg) = (println("$infile: $msg"); nothing)
   g = input_setting(lines, "g")
   (g === nothing || parse(Float64, g) == dambreak_g) || return bad("g must be $dambreak_g.")
   tf = lowercase(something(input_setting(lines, "Topog function"), ""))
   if sol.s == 0
      tf == "flat" || return bad("the topography must be flat.")
   else
      tf == "xslope" || return bad("the topography must be xslope.")
      p = parse(Float64, strip(input_setting(lines, "Topog params"), ['(', ')', ' ']))
      p == -sol.s || return bad("Topog params must be $(-sol.s).")
   end
   drag = lowercase(something(input_setting(lines, "Drag"), ""))
   if sol.mu == 0
      # Chezy drag with Co = 1e-8 is what record_reference.sh substitutes for
      # commits without Drag = None.
      chezy0 = drag == "chezy" && parse(Float64, something(input_setting(lines, "Chezy Co"), "1")) <= 1e-8
      drag == "none" || chezy0 || return bad("Drag must be None.")
   else
      drag == "coulomb" || return bad("Drag must be Coulomb.")
      parse(Float64, input_setting(lines, "Coulomb Co")) == sol.mu ||
         return bad("Coulomb Co must be $(sol.mu).")
   end
   gf = input_setting(lines, "Geometric factors")
   k = dambreak_constants(sol; geometric = (gf === nothing) || lowercase(gf) != "off")
   if sol.mu > 0 && (k.A <= 0 || sol.hr > 0)
      return bad("Coulomb friction needs a slope steeper than the friction angle, and a dry bed.")
   end
   (sol.u0 >= 0 && (sol.u0 == 0 || sol.hr == 0)) ||
      return bad("an initial velocity must be downhill, with a dry bed.")
   return k
end

# Write the exact initial condition of the dam break sol (see write_exact_ic).
# Returns true on error.
function make_dambreak_ic(testname, dim, sol)
   dim == 1 || (println("Dam breaks are 1D."); return true)
   k = dambreak_check(testname, sol)
   k === nothing && return true
   xb = (sol.hr > 0) ? sol.R : 0.0
   tmplblock = "Cap:\ncapX = $((xb - sol.L) / 2)\ncapHeight = $(sol.hl)\ncapRadius = $((xb + sol.L) / 2)\n" *
               "capShape = flat\ncapU = 0\ncapConc = 0.0"
   return write_exact_ic(testname, dim, tmplblock,
                         (x, y, dx, dy) -> dambreak_cell_average(sol, k, x, dx, 0.0))
end

# The comparison with the dam break sol: exact, inwindow and U for
# window_errors, compared where the exact solution holds (see above), with U
# the largest exact speed at the end of the run. Nothing if the input file
# does not match sol.
function dambreak_comparison(testname, sol)
   k = dambreak_check(testname, sol)
   k === nothing && return nothing
   lines = readlines(test_inputfile(testname))
   tend = parse(Float64, input_setting(lines, "T end"))
   U = ((sol.hr > 0) ? k.um : 2 * k.cl) + sol.u0 + k.A * tend
   function inwindow(x, dx, t)
      ξa = x - dx / 2 - dambreak_shift(sol, k, t)
      ξb = x + dx / 2 - dambreak_shift(sol, k, t)
      ξa >= -sol.L + k.cl * t + sol.margin || return false
      return sol.hr == 0 || ξb <= sol.R - k.cr * t - sol.margin
   end
   return (x, y, dx, dy, t) -> dambreak_cell_average(sol, k, x, dx, t), inwindow, U
end
