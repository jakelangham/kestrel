# Fast curvature verification tests for runall.jl (see Scheme_curvature_note.tex
# and the full studies in valley.jl, valley2d.jl and bowl.jl). Each test function
# returns 0 on success and 1 on failure.

module Valley
include("valley.jl")
end

module Bowl
include("bowl.jl")
end

# Steady valley flow at dx = 1 m. The error against the reference for the
# selected system must be below tol, and the errors against the references of
# the other systems at least sep*tol, so that the test identifies the system.
# (The no-curvature system and the model of Langham et al. 2023 differ only in
# the volume they conserve, so their profiles are close.)
function test_curvature_valley(testname, prog, geofactors, tol, sep)
   systems = Dict("on" => "C  (1,-1)", "no curvature" => "B' (1, 0)")
   errs = redirect_stdout(devnull) do
      Valley.run_valley(abspath(prog), testname, geofactors, [1.0])
   end
   result = 0
   for (sname, e) in errs
      err = maximum(e[1])
      if sname == systems[geofactors]
         if err > tol
            println("Error against $sname is $err, more than $tol")
            result = 1
         end
      elseif err < sep * tol
         println("Error against $sname is only $err, less than $(sep * tol)")
         result = 1
      end
   end
   return result
end

# Frictionless run-up in a bowl at dx = 1 m with u0^2/(gR) = 0.81. The maximum
# mass-weighted height must be within ztol (relative) of the parcel limit, the
# physical mass must be constant to output precision and, for the corrected
# system, the maximum energy must not exceed E0 (1 + Etol).
function test_curvature_bowl(testname, prog, geofactors, ztol, Etol)
   system = geofactors == "on" ? :on : :nocurv
   res = redirect_stdout(devnull) do
      Bowl.run_bowl(abspath(prog), testname, geofactors, 0.81, [1.0])
   end[1]
   target = Bowl.parcel_height(0.81, system)
   result = 0
   if abs(res.zbmax / target - 1) > ztol
      println("Maximum height $(res.zbmax) m differs from $target m by more than $(100ztol)%")
      result = 1
   end
   if res.dM > 1e-8
      println("Physical mass changed by $(res.dM) (relative)")
      result = 1
   end
   if system == :on && res.Emax > 1 + Etol
      println("Energy rose to $(res.Emax) times its initial value")
      result = 1
   end
   return result
end
