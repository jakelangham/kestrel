include("testlib.jl")
using Printf

tests_1d = [
   ("flat_depositional", 1)
   ("cap_dilute", 1)
   ("cap_conc", 1)
   ("cap_morpho", 1)
   ("flux_hydro", 1)
   ("flux_edwards2019", 1)
   ("flux_morpho", 1)
]

tests_2d = [
   ("flat_depositional_2d", 2)
   ("cap_dilute_2d", 2)
   ("cap_conc_2d", 2)
   ("cap_morpho_2d", 2)
   ("flux_hydro_2d", 2)
   ("flux_edwards2019_2d", 2)
   ("flux_morpho_2d", 2)
   ("flux_single_pt", 2)
]

tests_noflow = [
   ("lake_at_rest_hydro", 1)
   ("lake_at_rest_morpho", 1)
   ("lake_at_rest_hydro_2d", 2)
   ("lake_at_rest_morpho_2d", 2)
]

tests_identical = [
   # These two tests check for independence of the tiling layout by comparing
   # simulations with different tile sizes. The static test is designed so
   # that the flow stays largely within the tiles added at the start, while
   # the dynamic test features a flow that enters inactive regions.
   # Comparisons are made up to some tolerance (4th arg), relative to the
   # largest magnitude of each field (see check_identical_simulations). Both
   # tests are bitwise independent of the tile layout, so their output must
   # be identical.
   ("tile_indep_static", "tile_indep_static_100m", "tile_indep_static_50m", 0.0, 2)
   ("tile_indep_dynamic", "tile_indep_dynamic_100m", "tile_indep_dynamic_20m", 0.0, 2)
   # As tile_indep_static, on a curved bed (a rotated Gaussian valley), which
   # exercises the curvature terms across tile boundaries. It is also bitwise
   # independent of the tile layout.
   ("tile_indep_curved", "tile_indep_curved_100m", "tile_indep_curved_50m", 0.0, 2)
]

# Tests of the curvature of the bed ('Geometric factors' and 'Curvature
# friction'). The flows are thin enough that each part of them should move
# like a point mass sliding on the bed, so their reference values are those of
# a point mass (see the test functions in testlib.jl).
tests_curvature = [
   # Frictionless run-up in the bowl b = x^2/(2R), R = 100 m, launched from
   # the bottom with speed u0, u0^2/(gR) = 0.81. Maximum height: with the
   # curvature force, u0^2/(2g) = 40.5 m; without it ('no curvature', the
   # model of Langham et al. 2023), (R/2)(exp(u0^2/(gR)) - 1) = 62.395 m.
   # (2nd arg: reference height, 3rd: relative tolerance, 4th: whether the
   # energy must be conserved.)
   ("curvature_bowl_on", test_curvature_runup, 40.5, 0.01, true)
   ("curvature_bowl_nocurv", test_curvature_runup, 62.395, 0.01, false)
   # As above, with Coulomb friction (mu = 0.2) and the centripetal term in
   # the normal stress ('Curvature friction = On'): 21.692 m, against
   # 26.061 m without the term.
   ("curvature_bowl_friction_on", test_curvature_runup, 21.692, 0.01, false)
   # Frictionless layer crossing a lateral bank that fades out along the flow,
   # so that b_xy != 0: direction of the mean velocity (deg) and centroid (m)
   # at the end, with tolerances. Without the curvature force the direction
   # is 1.35 deg off and the centroid 15 mm.
   ("curvature_bank_fade_on", test_curvature_bank, 10.852, 0.25, (0.38812, -0.02319), 0.003)
   # Slump against a steep wall that the grid does not resolve: the time step
   # must not collapse through runaway velocities in thin cells (2nd arg:
   # largest speed allowed, m/s, in cells deeper than 0.1 mm; 3rd: whether the
   # total volume must be conserved).
   ("curvature_wall_on", test_curvature_wall, 10.0, false)
   # As above, for a debris-like flow with drag on an erodible bed.
   ("curvature_wall_morpho_on", test_curvature_wall, 10.0, true)
]

# Check if NetCDF support enabled and define extra test if so.
tests_netcdf = []
with_netcdf = links_to_netcdf(prog)
if with_netcdf
   tests_netcdf = [
      ("netcdf_restart", 2)
   ]
end

# Check if built with OpenMP and define thread independence tests if so. Each
# simulation is run with 1 thread and (twice) with nthreads threads, and the
# outputs must be identical. These cover 1D/2D, hydraulic/morphodynamic flows
# and (in tile_indep_dynamic_20m) many small tiles being added during the flow.
tests_threads = []
with_openmp = built_with_openmp(prog)
if with_openmp
   nthreads = 4
   tests_threads = [
      ("flux_morpho", 1, nthreads)
      ("flux_hydro_2d", 2, nthreads)
      ("flux_morpho_2d", 2, nthreads)
      ("tile_indep_dynamic_20m", 2, nthreads)
   ]
end

function run_test(i, name, testfunc, args...)
   println("Test #$i: ($name)")
   res = testfunc(args...)
   if res == 1
      printstyled("FAIL\n"; color=:red, bold=true)
   else
      printstyled("PASS\n"; color=:green)
   end
   return 1 - res
end

function run_1d()
   println("Running 1D tests...")
   numpassed = 0
   for i = 1:length(tests_1d)
      args = tests_1d[i]
      testname = tests_1d[i][1]
      testfunc = test_flow_consistency
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_2d()
   println("\nRunning 2D tests...")
   numpassed = 0
   for i = 1:length(tests_2d)
      testname = tests_2d[i][1]
      testfunc = test_flow_consistency
      args = tests_2d[i]
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_noflow()
   println("\nStatic flow tests...")
   numpassed = 0
   for i = 1:length(tests_noflow)
      testname = tests_noflow[i][1]
      testfunc = test_no_flow
      args = tests_noflow[i]
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_identical()
   println("\nTile independence tests...")
   numpassed = 0
   for i = 1:length(tests_identical)
      testname = tests_identical[i][1]
      testfunc = test_identical_simulations
      args = tests_identical[i]
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_curvature()
   println("\nCurvature tests...")
   numpassed = 0
   for i = 1:length(tests_curvature)
      testname = tests_curvature[i][1]
      testfunc = tests_curvature[i][2]
      args = (testname, tests_curvature[i][3:end]...)
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_netcdf()
   println("\nNetCDF tests...")
   numpassed = 0
   for i = 1:length(tests_netcdf)
      testname = tests_netcdf[i][1]
      testfunc = test_netcdf_restart
      args = tests_netcdf[i]
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_threads()
   println("\nThread independence tests...")
   numpassed = 0
   for i = 1:length(tests_threads)
      testname = tests_threads[i][1]
      testfunc = test_thread_independence
      args = tests_threads[i]
      numpassed += run_test(i, testname, testfunc, args...)
   end
   return numpassed
end

function run_all()

   numtests = length(tests_1d) + length(tests_2d) + length(tests_noflow) +
         length(tests_identical) + length(tests_curvature) + length(tests_netcdf) +
         length(tests_threads)
   numpassed = 0

   print("Running full test suite. Some tests (espcially 2D ones) may take a",
         " while\n\n")
   
   numpassed += run_1d()
   
   numpassed += run_2d()
   
   numpassed += run_noflow()

   numpassed += run_identical()

   numpassed += run_curvature()

   if with_netcdf
      numpassed += run_netcdf()
   end

   if with_openmp
      numpassed += run_threads()
   else
      println("\nKestrel is not built with OpenMP: skipping thread independence tests.")
   end

   report_pass_rate_and_exit(numpassed, numtests)
end

function report_pass_rate_and_exit(numpassed, numtests)
   ratio = 100 * numpassed / numtests
   @printf "Results: %i/%i passed" numpassed numtests
   @printf "%s (%.1f%c)\n" (numpassed == numtests ? "!" : "") ratio '%'

   if numpassed == numtests
      exit()
   else
      exit(1)
   end
end

# Print a list of available test categories.
function print_options()
   println("Options are 'all', '1d', '2d', 'noflow', 'identical', 'curvature', 'netcdf' and 'threads'.")
end

tests = (if isempty(ARGS); ["all"]; else ARGS; end)
if "help" in tests
   print_options()
elseif "all" in tests
   run_all()
else
   numtests = 0
   numpassed = 0
   if "1d" in tests
      numtests += length(tests_1d)
      numpassed += run_1d()
   end
   if "2d" in tests
      numtests += length(tests_2d)
      numpassed += run_2d()
   end
   if "noflow" in tests
      numtests += length(tests_noflow)
      numpassed += run_noflow()
   end
   if "identical" in tests
      numtests += length(tests_identical)
      numpassed += run_identical()
   end
   if "curvature" in tests
      numtests += length(tests_curvature)
      numpassed += run_curvature()
   end
   if "netcdf" in tests
      if with_netcdf
         numtests += length(tests_netcdf)
         numpassed += run_netcdf()
      else
         printstyled("FAIL Kestrel is not built with NetCDF support\n";
                     color = :red, bold = true)
      end
   end
   if "threads" in tests
      if with_openmp
         numtests += length(tests_threads)
         numpassed += run_threads()
      else
         printstyled("FAIL Kestrel is not built with OpenMP support\n";
                     color = :red, bold = true)
      end
   end

   if numtests >= 1
      report_pass_rate_and_exit(numpassed, numtests)
   else
      println("No valid tests selected.")
      print_options()
   end
end

exit(1) # default is to return error
