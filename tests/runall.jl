include("testlib.jl")
using Printf

tests_1d = [
   ("flat_depositional", 1)
   ("cap_dilute", 1)
   ("cap_conc", 1)
   ("cap_morpho", 1)
   ("flux_hydro", 1)
   ("flux_morpho", 1)
]

tests_2d = [
   ("flat_depositional_2d", 2)
   ("cap_dilute_2d", 2)
   ("cap_conc_2d", 2)
   ("cap_morpho_2d", 2)
   ("flux_hydro_2d", 2)
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
   # that all tiles are added at the start with none added during the flow,
   # while the dynamic test features a flow that enters inactive regions.
   # Comparisons are made up to some tolerance (4th arg). We don't insist on
   # perfect equivalence because changing the tile layout inevitably reorders
   # some arithmetic in the method. The dynamic test accepts a weaker
   # tolerance because tiles are added only when Hn >
   # RunParams%heightThreshold so it is possible in some cases for very thin
   # flows to run into inactive regions. This is not expected to
   # fundamentally matter, but may potentially to larger than machine eps
   # deviations.
   ("tile_indep_static", "tile_indep_static_100m", "tile_indep_static_50m", 1e-13, 2)
   ("tile_indep_dynamic", "tile_indep_dynamic_100m", "tile_indep_dynamic_20m", 1e-11, 2)
]

# Check if NetCDF support enabled and define extra test if so.
tests_netcdf = []
with_netcdf = links_to_netcdf(prog)
if with_netcdf
   tests_netcdf = [
      ("netcdf_restart", 2)
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

function run_all()

   numtests = length(tests_1d) + length(tests_2d) + length(tests_noflow) +
         length(tests_identical) + length(tests_netcdf)
   numpassed = 0

   print("Running full test suite. Some tests (espcially 2D ones) may take a",
         " while\n\n")
   
   numpassed += run_1d()
   
   numpassed += run_2d()
   
   numpassed += run_noflow()

   numpassed += run_identical()

   if with_netcdf
      numpassed += run_netcdf()
   end
  
   ratio = 100 * numpassed / numtests
   @printf "Results: %i/%i passed" numpassed numtests
   @printf "%s (%.1f%c)\n" (numpassed == numtests ? "!" : "") ratio '%'
end

tests = ( if isempty(ARGS); ["all"]; else ARGS; end)
if "all" in tests
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
   if "netcdf" in tests
      if with_netcdf
         numtests += length(tests_netcdf)
         numpassed += run_netcdf()
      else
         printstyled("FAIL Kestrel is not built with NetCDF support\n"; 
                     color = :red, bold = true)
      end
   end

   if numtests >= 1
      ratio = 100 * numpassed / numtests
      @printf "Results: %i/%i passed" numpassed numtests
      @printf "%s (%.1f%c)\n" (numpassed == numtests ? "!" : "") ratio '%'
   else
      println("No tests selected.  Options are 'all', '1d', '2d', 'noflow', 'identical' and 'netcdf'.")
   end

end
