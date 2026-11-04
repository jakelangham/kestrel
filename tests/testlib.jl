using DelimitedFiles

prog = "../src/kestrel"

# Basic checks for conservativity and positivity of the solution output when
# simulating with the input file Input_$testname.txt
# dim = spatial dimension of simulation
function test_flow_consistency(testname, dim)
   if run_simulation(testname)
      return 1
   end

   result = check_conservativity(testname)
   result |= check_positivity(testname, dim)
   result |= check_max_ero_depth(testname, dim)
   
   return result
end

# Run a simulation for which no flow is predicted and verify whether initial and
# final fields are the same.
function test_no_flow(testname, dim)
   if run_simulation(testname)
      return 1
   end

   result = check_no_flow(testname, dim)

   return result
end

# Take two simulations, testname1 & testname2 and check if they produce the same
# output. This can be used for checking that simulations are independent of
# settings that shouldn't matter, like tile layout.
function test_identical_simulations(testname, infile1, infile2, tol, dim)
   if run_simulation(infile1) || run_simulation(infile2)
      return 1
   end

   result = check_identical_simulations(infile1, infile2, tol, dim)

   return result
end

# Check that a simulation gives identical results regardless of the number of
# OpenMP threads used to run it. The simulation for testname is run with 1
# thread and twice with nthreads threads (the repeat is to catch
# nondeterministic results). All outputs must match exactly: within each
# parallel loop, every tile is computed using the same arithmetic as in serial,
# so any difference indicates a bug (e.g. a data race).
function test_thread_independence(testname, dim, nthreads)
   dir_serial = "$(testname)_threads1"
   dir_par = "$(testname)_threads$(nthreads)"
   dir_repeat = "$(testname)_threads$(nthreads)_repeat"

   if run_simulation_with_threads(testname, 1, dir_serial) ||
      run_simulation_with_threads(testname, nthreads, dir_par) ||
      run_simulation_with_threads(testname, nthreads, dir_repeat)
      return 1
   end

   result = check_identical_outputs(dir_serial, dir_par)
   result |= check_identical_outputs(dir_par, dir_repeat)

   return result
end

# This test runs the same simulation twice. In the second instance, an interrupt
# is called and the simulation is restarted. If both outputs are the same,
# return 0.
function test_netcdf_restart(testname, dim)
   # The test itself is farmed off to an external shell script ("$testname.sh"),
   # since this is easier and we don't need to do any numerical work on the
   # result.
   if run_test_script(testname)
      return 1
   else
      return 0
   end
end

# Run-up of a thin layer launched from the bottom of a 1D bowl. The maximum,
# over the outputs, of the mass-weighted height of the layer,
#    z_b = sum(gam Hn b) / sum(gam Hn),   gam = sqrt(1 + b_x^2),
# must be within rtol (relative) of zmax, the height reached by a point mass
# sliding on the bed, and the mass sum(gam Hn) must be constant. For a
# frictionless flow with the curvature force (check_energy = true), the energy
#    E = sum(gam Hn |u|^2/2 + g gam Hn b + g Hn^2/2),   |u|^2 = u^2 gam^2
# (per unit density and cell width) must not rise by more than 1%.
function test_curvature_runup(testname, zmax, rtol, check_energy)
   if run_simulation(testname)
      return 1
   end

   g = input_value(test_inputfile(testname), "g", 9.81)
   zbmax = -Inf
   M = Float64[]
   E = Float64[]
   for f in sort(filter(is_resultfile_txt, readdir(testname)))
      data, header = readdlm("$testname/$f", ',', header = true)
      col(name) = Float64.(data[:, findfirst(==(name), strip.(vec(header)))])
      Hn, u, b, bx = col("flow_depth"), col("x_velocity"), col("base_elevation"), col("x_slope")
      gam = sqrt.(1 .+ bx .^ 2)
      push!(M, sum(gam .* Hn))
      push!(E, sum(gam .* Hn .* (0.5 .* u .^ 2 .* gam .^ 2 .+ g .* b) .+ 0.5 .* g .* Hn .^ 2))
      zbmax = max(zbmax, sum(gam .* Hn .* b) / M[end])
   end

   result = 0
   println("Maximum mass-weighted height: $zbmax m (point mass: $zmax m)")
   if abs(zbmax / zmax - 1) > rtol
      println("This differs from the point mass by more than $(100 * rtol)%")
      result = 1
   end
   dM = (maximum(M) - minimum(M)) / M[1]
   println("Relative change in mass:      $dM")
   if dM > 1e-8
      result = 1
   end
   if check_energy
      dE = maximum(E) / E[1] - 1
      println("Largest relative rise in energy: $dE")
      if dE > 0.01
         result = 1
      end
   end

   return result
end

# Thin layer moving across a bank in 2D. At the last output, the direction of
# the mass-weighted mean horizontal velocity (degrees from the x axis) must be
# within angtol of angle, and the centroid of the layer within postol (m) of
# centroid, these being the values for point masses sliding on the bed. (Each
# cell is weighted by its mass gam Hn, with gam = sqrt(1 + b_x^2 + b_y^2).) The
# volume must also be conserved.
function test_curvature_bank(testname, angle, angtol, centroid, postol)
   if run_simulation(testname)
      return 1
   end

   f = sort(filter(is_resultfile_txt, readdir(testname)))[end]
   data, header = readdlm("$testname/$f", ',', header = true)
   col(name) = Float64.(data[:, findfirst(==(name), strip.(vec(header)))])
   x, y, Hn, u, v = col("x"), col("y"), col("flow_depth"), col("x_velocity"), col("y_velocity")
   m = sqrt.(1 .+ col("x_slope") .^ 2 .+ col("y_slope") .^ 2) .* Hn
   dir = atand(sum(m .* v), sum(m .* u))
   dist = hypot(sum(m .* x) / sum(m) - centroid[1], sum(m .* y) / sum(m) - centroid[2])

   result = 0
   println("Direction of the mean velocity: $dir deg (point masses: $angle deg)")
   if abs(dir - angle) > angtol
      println("This differs from the point masses by more than $angtol deg")
      result = 1
   end
   println("Distance of the centroid from the point masses': $dist m")
   if dist > postol
      println("This is more than $postol m")
      result = 1
   end
   result |= check_volume_change(testname, 1e-8)

   return result
end

# Layer slumping against a steep wall. The simulation must reach its end time
# (within a wall-clock limit, since a collapsing time step would otherwise
# stall it) and the speed of the fluid, in cells deeper than 0.1 mm, must not
# exceed umax. If check_volume is set, the total (flow + bed) volume must be
# conserved.
function test_curvature_wall(testname, umax, check_volume)
   if run_simulation_with_timeout(testname, 120.0)
      return 1
   end

   speed = 0.0
   for f in filter(is_resultfile_txt, readdir(testname))
      data, header = readdlm("$testname/$f", ',', header = true)
      col(name) = Float64.(data[:, findfirst(==(name), strip.(vec(header)))])
      wet = col("flow_depth") .> 1e-4
      any(wet) && (speed = max(speed, maximum(abs.(col("x_velocity")[wet]))))
   end

   result = 0
   println("Largest speed: $speed m/s")
   if speed > umax
      println("This is more than $umax m/s")
      result = 1
   end
   if check_volume
      result |= check_conservativity(testname)
   end

   return result
end

# Run $prog on the input file corresponding to the test with name testname and
# log the output in "$testname.log". Return a helpful error if the command
# fails.
function run_simulation(testname)
   infile = test_inputfile(testname)
   logfile = "$testname.log"

   error = false
   t = @async begin
      try
         run(pipeline(`$prog $infile`, stdout = logfile, stderr = logfile));
      catch
         println("Simulation exited unexpectedly.")
         println("The command was $prog $infile.")
         println("Any output has been logged in $logfile.")
         error = true
      end
   end
   wait_for_task(t)

   return error
end

# Like run_simulation, but stop the simulation, and return an error, if it runs
# for longer than timeout (s) of wall-clock time.
function run_simulation_with_timeout(testname, timeout)
   infile = test_inputfile(testname)
   logfile = "$testname.log"

   p = run(pipeline(`$prog $infile`, stdout = logfile, stderr = logfile), wait = false)
   t0 = time()
   while process_running(p) && time() - t0 < timeout
      sleep(0.1)
   end
   if process_running(p)
      kill(p)
      println("Simulation did not finish within $timeout s and was stopped.")
      println("Any output has been logged in $logfile.")
      return true
   end
   if !success(p)
      println("Simulation exited unexpectedly.")
      println("The command was $prog $infile.")
      println("Any output has been logged in $logfile.")
      return true
   end

   return false
end

# Like run_simulation, but with nthreads OpenMP threads and writing the output to
# outdir. This is done via a temporary copy of the input file (with its output
# directory replaced), which is removed afterwards.
function run_simulation_with_threads(testname, nthreads, outdir)
   infile = test_inputfile(testname)
   tmpfile = test_inputfile(outdir)
   logfile = "$outdir.log"

   replaced = false
   open(tmpfile, "w") do io
      for line in readlines(infile)
         if occursin(r"^\s*directory\s*="i, line)
            println(io, "directory = $outdir/")
            replaced = true
         else
            println(io, line)
         end
      end
   end
   if !replaced
      println("No output directory found in $infile.")
      rm(tmpfile)
      return true
   end
   rm(outdir, force = true, recursive = true)

   cmd = addenv(`$prog $tmpfile`, "OMP_NUM_THREADS" => string(nthreads),
                "OMP_WAIT_POLICY" => "passive")
   error = false
   t = @async begin
      try
         run(pipeline(cmd, stdout = logfile, stderr = logfile));
      catch
         println("Simulation exited unexpectedly.")
         println("The command was OMP_NUM_THREADS=$nthreads $prog $tmpfile.")
         println("Any output has been logged in $logfile.")
         error = true
      end
   end
   wait_for_task(t)
   rm(tmpfile)

   return error
end

# Like run_simulation, but for cases where the test is deferred to a bespoke
# shell script, assumed to be called "$testname.sh".
function run_test_script(testname)
   script = "$testname.sh"
   logfile = "$testname.log"

   error = false
   t = @async begin
      try
         run(pipeline(`./$script`, stdout = logfile, stderr = logfile));
      catch
         println("External test script exited unexpectedly.")
         println("The command was $script.")
         println("Any output has been logged in $logfile.")
         error = true
      end
   end
   wait_for_task(t)

   return error
end

# Check that simulation (.txt) outputs in directory dir conserved their total
# volume and solids volume by scrutinising the output of Volume.txt.
function check_conservativity(dir)
   vfile = "$dir/Volume.txt"
   infile = test_inputfile(dir)

   # Tolerance for relative error. We typically lose 3 - 4 digits of precision
   # by conducting arithmetic with values over roughly this many orders of
   # magntide, so this seems like an acceptable choice.
   rtol = 1e-10

   # Retrieve any volume input by flux sources (total, solids)
   flux_vol, flux_sol = total_flux_sources(infile)

   # Get volume data from Volume.txt using calls to awk. There are native Julia
   # ways to do this but they require libraries so less portable.
   #
   # Flow volumes (at start/end of simulation)
   vol_t0 = parse(Float64, replace(readchomp(`awk 'NR==2 {print $2}' $vfile`), [','] => ""))
   vol_tn = parse(Float64, replace(readchomp(`awk 'END {print $2}' $vfile`), [','] => ""))
   # Bed volumes
   b_t0 = parse(Float64, replace(readchomp(`awk 'NR==2 {print $3}' $vfile`), [','] => ""))
   b_tn = parse(Float64, replace(readchomp(`awk 'END {print $3}' $vfile`), [','] => ""))
   # Solids masses
   solm_t0 = parse(Float64, replace(readchomp(`awk 'NR==2 {print $6}' $vfile`), [','] => ""))
   solm_tn = parse(Float64, replace(readchomp(`awk 'END {print $6}' $vfile`), [','] => ""))
   # Bed solids masses
   bsm_t0 = parse(Float64, replace(readchomp(`awk 'NR==2 {print $7}' $vfile`), [','] => ""))
   bsm_tn = parse(Float64, replace(readchomp(`awk 'END {print $7}' $vfile`), [','] => ""))

   expected_vol = flux_vol + vol_t0 + b_t0
   final_vol = vol_tn + b_tn
   vol_err = abs((expected_vol - final_vol) / expected_vol)
   println("Relative error in volume:     $vol_err")

   # Volume.txt outputs solids mass, so to get volume must divide by density.
   ρs = solids_density(infile)

   expected_solids_vol = flux_sol + (solm_t0 + bsm_t0) / ρs
   final_solids_vol = (solm_tn + bsm_tn) / ρs
   if (expected_solids_vol == 0.0)
      # In case initial condition is dilute
      solid_err = abs(expected_solids_vol - final_solids_vol)
      println("Absolute error in solid mass: $solid_err")
   else
      solid_err = abs((expected_solids_vol - final_solids_vol) / expected_solids_vol)
      println("Relative error in solid mass: $solid_err")
   end

   if (vol_err > rtol || solid_err > rtol)
      return 1
   else
      return 0
   end
end

# Hunt for the flux sources specified as time series in the input file. Then
# integrate over each source between the limits of the simulation (t start/end),
# sum them up and return both the total volume flux and solid volume flux.
function total_flux_sources(infile)
   search_sources = `grep -i "Source" $infile`
   source_lines = []
   try
      source_lines = lowercase.(readlines(search_sources))
   catch
      return 0, 0 # no sources found
   end

   source_lines = strip_comments(source_lines)

   # Split up into blocks
   blocks = []
   block_start = 1
   block_end = 0
   for i = 1:length(source_lines)
      if occursin("source:", source_lines[i]) && block_end > 0
         push!(blocks, source_lines[block_start:block_end])
         block_start = i
      end
      block_end += 1
   end
   push!(blocks, source_lines[block_start:block_end])

   # Get start and end time
   tstart, tend = simulation_interval(infile)

   # For each source block, read its time series data.
   t = []; Q = []; ψ = []
   Qtot = 0; Qψtot = 0
   for b in blocks
      for i = 1:length(b)
         if occursin("sourcetime", b[i])
            t = extract_time_series(b[i])
         elseif occursin("sourceflux", b[i])
            Q = extract_time_series(b[i])
         elseif occursin("sourceconc", b[i])
            ψ = extract_time_series(b[i])
         end
      end
      # Integrate and add to total flux
      Q, Qψ = integrate_source_time_series(t, Q, ψ, tstart, tend)
      Qtot += Q
      Qψtot += Qψ
   end

   return Qtot, Qψtot
end

# Assuming a string of the form "BLAH = (f1, f2, f3, etc)", extract its time
# series to an array of Float64's.
function extract_time_series(str)
   return parse.(Float64, split(strip(split(str, "=")[2])[2:end-1], ","))
end

# Retrieve the simulation interval [tstart, tend] from the given input file.
function simulation_interval(infile)
   tstart_d = 0 # default value
   tstart = tstart_d; tend = 0

   try
      tstart_str = readlines(`grep -i "t start" $infile`)
      tstart = parse(Float64, strip(split(tstart_str[1], "=")[2]))
   catch
   end

   # Assume t end present, since otherwise simulation wouldn't have run.
   tend_str = readlines(`grep -i "t end" $infile`)
   tend = parse(Float64, strip(split(tend_str[1], "=")[2]))

   return tstart, tend
end

# Given timeseries data t = (t1, t2, ...) and variable var = (v1, v2, ...),
# integrate var between t = tstart and t = tend, assuming piecewise linear.
function integrate_source_time_series(t, Q, ψ, tstart, tend)
   total_flux = 0; solids_flux = 0

   # Series can be one point, which is a special case
   if (length(t) == 1 && t[1] < tend)
      duration = min(tend - tstart, tend - t[1])
      total_flux = Q[1] * duration
      solids_flux = ψ[1] * total_flux
   end

   # Loop through each interval in time series. There is a bit of logic to find
   # the bounds on integration so that they don't lie outside [tstart,tend].
   for i = 1:(length(t)-1)
      if (t[i+1] < tstart) || (t[i] > tend)
         continue # outside simulation interval
      end

      dQdt = (Q[i+1] - Q[i]) / (t[i+1] - t[i])
      dψdt = (ψ[i+1] - ψ[i]) / (t[i+1] - t[i])
      t_lower = t[i]; t_upper = t[i+1]
      Q_lower = Q[i]; Q_upper = Q[i+1]
      ψ_lower = ψ[i]; ψ_upper = ψ[i+1]

      if (t[i] < tstart)
         t_lower = tstart
         Q_lower = Q[i] + dQdt * (tstart - t[i])
         ψ_lower = ψ[i] + dψdt * (tstart - t[i])
      end
      if (t[i+1] > tend)
         t_upper = tend
         Q_upper = Q[i] + dQdt * (tend - t[i])
         ψ_upper = ψ[i] + dψdt * (tend - t[i])
      end
      # Trapezium rule for Q
      Δt = t_upper - t_lower
      total_flux += Δt * (Q_lower + Q_upper) / 2
      # Qψ is quadratic formed by multiplication of two linear pieces, so may be
      # determined explicitly too
      solids_flux += (Δt/6) * 
              (Q_lower * ψ_upper + Q_upper * ψ_lower +
          2 * (Q_lower * ψ_lower + Q_upper * ψ_upper))
   end

   return total_flux, solids_flux
end

# Check that all (.txt) solutions in directory dir have positive flow depth
# everywhere.
function check_positivity(dir, dim)
   for f in readdir(dir)
      if is_resultfile_txt(f)
         Hnf = Hn_data_field(dim)
         awkcmd = `awk -F ',' '(FNR > 1) && (NF > 0) && ($'$Hnf' < -1e-14) \
                   { print $'$Hnf' }' ./$dir/$f`
         if !isempty(read(awkcmd))
            println("Solution depth field is negative somewhere in file",
                    " ./$dir/$f")
            return 1
         end
      end
   end

   println("Solution depth is everywhere positive")
   return 0
end

# Check that none of the (.txt) solutions in directory dir exceed the maximum
# erosion depth.
function check_max_ero_depth(dir, dim)
   infile = test_inputfile(dir)

   # Get the maximum erosion depth
   ero_dep = erosion_depth(infile)
   #ero_dep = run_settings(infile, "erosion depth", Float64, 0.0)

   # Check each file
   for f in readdir(dir)
      if is_resultfile_txt(f)
         btf = bt_data_field(dim)
         awkcmd = `awk -F ',' '(FNR > 1) && (NF > 0) && ($'$btf' < -'$ero_dep') \
                   { print $'$btf' }' ./$dir/$f`
         if !isempty(read(awkcmd))
            println("Erosion exceeds maximum depth somewhere in file ./$dir/$f")
            return 1
         end
      end
   end

   println("Erosion nowhere exceeds maximum depth")
   return 0
end

# Check each result file produced by test testname against "000000.txt". Return
# 1 (fail) if depth column differs for any file, 0 otherwise.
function check_no_flow(testname, dim)
   result = 0
   f0 = "000000.txt"

   cd(testname)
   for fn in readdir()
      if is_resultfile_txt(fn) && fn != f0
         Hnf = Hn_data_field(dim)
         # Store columns 2 and Hnf from f0 in a, then check if they are to be
         # found in fn. i.e. verify that height data from two files is equal
         checkcmd = `awk -F ',' '(NR==FNR) { a[$2$'$Hnf']; next } \
                                !($2$'$Hnf' in a) { print $2$'$Hnf' }' $f0 $fn`
         if !isempty(read(checkcmd))
            result = 1
            println("Solution depth field in $testname/$fn differs from ",
                     "initial condition")
            break
         end
      end
   end
   cd("..")

   if result == 0
      println("Initial condition and subsequent depth fields agree")
   end
   return result
end

# Check if two simulations with (.txt) result files in dirs testname1 and
# testname2 are identical, up to tolerance tol. We use a maximum norm to make
# this comparison.
function check_identical_simulations(testname1, testname2, tol, dim)
   result = 0 # return value

   # Compare Hn, u, (v in 2D), Hnψ & bt. Cells are matched by their
   # coordinates. A cell output by only one simulation (because its tile was
   # inactive in the other) is compared against a dry cell.
   #
   # The text output holds 11 significant figures, so values that differ only
   # by round-off can differ by one unit in the last printed digit, i.e. by
   # up to about 1e-10 relative to their magnitude. The tolerance tol is
   # therefore relative: for each field, the difference must not exceed tol
   # times the largest magnitude of that field in the two outputs.
   fields = [("Hn", Hn_data_field(dim)), ("u", u_data_field(dim))]
   if dim == 2
      push!(fields, ("v", v_data_field(dim)))
   end
   push!(fields, ("Hnψ", Hnψ_data_field(dim)), ("bt", bt_data_field(dim)))
   coords = (dim == 1) ? [2] : [2, 3]

   dir2_contents = readdir(testname2)
   for f in readdir(testname1)
      if is_resultfile_txt(f)
         if !(f in dir2_contents)
            result = 1
            println("Directories contain nonequal sets of solution data")
            break
         end

         f1 = "$testname1/$f"
         f2 = "$testname2/$f"

         soln1 = readdlm(f1, ',', skipstart = 1)
         soln2 = readdlm(f2, ',', skipstart = 1)

         # Pair up the rows of the two outputs by cell coordinates (rounded,
         # in case the last printed digit of a coordinate differs).
         key(soln, i) = Tuple(round.(Float64.(soln[i, coords]), digits = 6))
         rows2 = Dict(key(soln2, i) => i for i in 1:size(soln2, 1))
         pairs = Tuple{Int,Int}[]
         matched2 = falses(size(soln2, 1))
         for i in 1:size(soln1, 1)
            j = get(rows2, key(soln1, i), 0)
            push!(pairs, (i, j))
            j > 0 && (matched2[j] = true)
         end
         append!(pairs, [(0, j) for j in findall(.!matched2)])

         # Field value of row i (0 = a dry cell absent from that output).
         val(soln, i, c) = (i == 0) ? 0.0 : Float64(soln[i, c])

         for (name, c) in fields
            scale = max(maximum(abs, Float64.(soln1[:, c])),
                        maximum(abs, Float64.(soln2[:, c])), floatmin(Float64))
            maxdiff, worst = 0.0, (0, 0)
            for (i, j) in pairs
               # A cell missing from one output is only compared by its
               # depth, since the other fields of a dry cell are arbitrary.
               (i == 0 || j == 0) && name != "Hn" && continue
               d = abs(val(soln1, i, c) - val(soln2, j, c))
               if d > maxdiff
                  maxdiff, worst = d, (i, j)
               end
            end
            if maxdiff > tol * scale
               i, j = worst
               where = (i > 0) ? soln1[i, coords] : soln2[j, coords]
               println("Solution field in $f1 differs from $f2: $name differs by ",
                       maxdiff, " (", maxdiff / scale, " of its largest magnitude) at ",
                       Tuple(where))
               result = 1
            end
         end
         result == 1 && break
      end
   end

   return result
end

# Check that the output directories dir1 and dir2 contain the same files with
# exactly the same contents, apart from RunInfo.txt (which records the date and
# input file name).
function check_identical_outputs(dir1, dir2)
   files1 = sort(filter(f -> f != "RunInfo.txt", readdir(dir1)))
   files2 = sort(filter(f -> f != "RunInfo.txt", readdir(dir2)))
   if files1 != files2
      println("Directories $dir1 and $dir2 contain different sets of files")
      return 1
   end
   for f in files1
      if read("$dir1/$f") != read("$dir2/$f")
         println("$dir1/$f differs from $dir2/$f")
         return 1
      end
   end

   return 0
end

# Retrieve the solids density ρs from the input file.
function solids_density(infile)
   ρs = 2000.0 # Default value
   ρs_str = ""

   # See if there is a user-specified value
   try
      ρs_str = readlines(`grep -i "rhos" $infile`)[1]
   catch
      return ρs
   end

   ρs = parse(Float64, strip(split(ρs_str, "=")[2]))
   return ρs
end

# Clearly these funcs can be generalised...
function erosion_depth(infile)
   ero_dep = 0.0 # Default value
   ero_dep_str = ""

   try
      ero_dep_str = readlines(`grep -i "erosion depth" $infile`)[1]
   catch
      return ero_dep
   end

   ero_dep = parse(Float64, strip(split(ero_dep_str, "=")[2]))
   return ero_dep
end

# Retrieve the numerical value of the setting key (e.g. "g") from the input
# file, or default if it is not set.
function input_value(infile, key, default)
   for line in readlines(infile)
      parts = split(split(line, "%")[1], "=")
      if length(parts) == 2 && lowercase(strip(parts[1])) == lowercase(key)
         return parse(Float64, strip(parts[2]))
      end
   end

   return default
end

# Check that the flow volume in Volume.txt, for the simulation output in dir,
# changed by at most rtol (relative).
function check_volume_change(dir, rtol)
   rows = filter(!isempty, strip.(readlines("$dir/Volume.txt")))[2:end]
   volume(row) = parse(Float64, strip(split(row, ",")[2]))
   err = abs(volume(rows[end]) / volume(rows[1]) - 1)
   println("Relative change in volume:    $err")

   return (err > rtol) ? 1 : 0
end

function test_inputfile(testname)
   return "Input_$testname.txt"
end

# Check if the string s names a .txt result file.
function is_resultfile_txt(s)
   if length(s) <= 4
      return false
   end
   if s[end-3:end] != ".txt"
      return false
   end
   for c in s[1:end-4]
      if !isnumeric(c)
         return false
      end
   end

   return true
end

# Return the .txt file output field for different fields, which depends on 
# whether the simulation is 1 or 2D.
function Hn_data_field(dim)
   return (dim == 1) ? 3 : 6
end
function u_data_field(dim)
   return (dim == 1) ? 5 : 8
end
function v_data_field(dim)
   return (dim == 1) ? nothing : 9
end
function Hnψ_data_field(dim)
   return (dim == 1) ? 9 : 14
end
function bt_data_field(dim)
   return (dim == 1) ? 12 : 17
end

# Remove anything after first '%' for array of strings
function strip_comments(strings)
   stripped = []
   for i = 1:length(strings)
      s = split(strings[i], "%")[1]
      if !isempty(s)
         push!(stripped, s)
      end
   end

   return stripped
end

# Verify if the executable prog links to the NetCDF libraries.
function links_to_netcdf(prog)
    lib_printer_exec = "ldd" # default library checker
    if Sys.isapple()
        lib_printer_exec = ["otool", "-L"]
    end

    cmd = `$lib_printer_exec $prog`
    cmd_output = try
        read(cmd, String)
    catch
        println("Warning: unable to look for netcdf support with command: $cmd")
        println("Continuing without testing netcdf features.")
        return false
    end

    return match(r"netcdf", cmd_output) !== nothing
end

# Check whether the executable prog was built with OpenMP, which it reports on
# start up (before it fails to find the given input file).
function built_with_openmp(prog)
   cmd_output = try
      read(pipeline(ignorestatus(`$prog no_such_input_file.txt`), stderr = devnull), String)
   catch
      println("Warning: unable to check whether $prog was built with OpenMP.")
      println("Continuing without testing OpenMP features.")
      return false
   end

   return occursin("OpenMP thread", cmd_output)
end

# Print a pending animation while waiting for task t to complete.
function wait_for_task(t)
   while !istaskdone(t)
      print("/\r"); sleep(0.1)
      print("-\r"); sleep(0.1)
      print("\\\r"); sleep(0.1)
      print("|\r"); sleep(0.1)
   end
end
