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

   # Check that Hn, u, Hnψ & bt are all equal
   Hnf = Hn_data_field(dim)
   uf = u_data_field(dim)
   Hnψf = Hnψ_data_field(dim)
   btf = bt_data_field(dim)

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

         # Compute the maximum difference of Hn,u,Hnψ for the two solutions,
         # sorted by their Hn value, assumed to have enough digits of precision
         # to be unique. (Though not completely rigorous, this avoids pain
         # screening out any the empty cells that are not common between f1 and
         # f2.)
         idx = [Hnf,uf,Hnψf,btf]
         fields1 = sortslices(soln1[:, idx], dims = 1)
         fields2 = sortslices(soln2[:, idx], dims = 1)
         nfields = min(size(fields1)[1], size(fields2)[1])
         diff = maximum(abs.(fields1[end-nfields+1:end] - fields2[end-nfields+1:end]))

         if (diff > tol)
            result = 1
            println("Solution field in $f1 differs from $f2")
            break
         end
      end
   end

   return result
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
        lib_printer_exec = "otool -L"
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

# Print a pending animation while waiting for task t to complete.
function wait_for_task(t)
   while !istaskdone(t)
      print("/\r"); sleep(0.1)
      print("-\r"); sleep(0.1)
      print("\\\r"); sleep(0.1)
      print("|\r"); sleep(0.1)
   end
end
