# Accuracy benchmarks: Kestrel's errors against exact solutions, and their
# convergence under grid refinement. Unlike the tests in ../tests, these have no
# pass/fail tolerances: they report the errors, the observed orders of
# convergence, and the change from the errors of a reference commit (by
# default the head of main) on the same grids.
#
# Run from this directory:
#    julia runbench.jl [selection ...] [--levels N] [--reference REF] [--save FILE]
#                      [--markdown FILE]
#
#    selection            'all' (the default), '1d', '2d' or the name of a case
#                         (see `cases` below); several may be given, as in
#                         ../tests/runall.jl. 'help' (or --help) describes
#                         the options and the cases.
#    --levels N           run N grid levels per case: the input file's grid, then
#                         each refinement by 2 (default: each case's own number)
#    --reference REF      what to compare with: a git ref (default: main), a
#                         reference file (as written by --save), or 'none'. A
#                         git ref's errors are kept in references/<hash>.txt; if
#                         that file does not exist, record_reference.sh builds
#                         the commit and records it first (several minutes).
#    --save FILE          write this run's errors to FILE, in the format read
#                         by --reference
#    --markdown FILE      also write the results table to FILE as Markdown
#
# Each run is reported as it finishes, with its status: OK, WORSE (an error
# more than 10% above its reference value), BETTER (both errors more than 10%
# below it), NEW (no reference) or FAIL (the run
# failed or gave a negative depth). A summary table and the number of runs
# that are not WORSE or FAIL follow. The exit status is 1 if any run is WORSE
# or FAIL, and 0 otherwise. Each run uses Input_<case>.txt with its
# cells per tile multiplied by 2^level; its output goes to <case>_n<cells>/.
# Set OMP_NUM_THREADS for OpenMP builds; the results do not depend on it.

include("benchlib.jl")

# The cases: (name, dim, solution, default number of levels). The solutions
# are in benchlib.jl, and must match the input files Input_<name>.txt:
# - Thacker (1981) oscillations in the bowl h0 (x^2 + y^2) / a^2 (see
#   thacker_frequency);
# - dam breaks, onto a dry bed (Ritter), a wet bed (Stoker) and, with Coulomb
#   friction, a dry slope (Mangeney et al.), released from rest and started
#   in motion (see dambreak_constants).
const thacker_bowl = (h0 = 10.0, a = 3000.0)
const dambreak = (kind = :dambreak, hl = 10.0, L = 1000.0, u0 = 0.0, margin = 100.0)
const cases = [
   ("thacker_planar", 1, (kind = :planar, thacker_bowl..., d = 500.0), 4)
   ("thacker_planar_2d", 2, (kind = :planar, thacker_bowl..., d = 500.0), 3)
   ("thacker_radial_2d", 2, (kind = :radial, thacker_bowl..., A = 0.5), 3)
   ("dambreak_dry", 1, (dambreak..., hr = 0.0, R = 0.0, s = 0.0, mu = 0.0), 4)
   ("dambreak_wet", 1, (dambreak..., hr = 2.0, R = 800.0, s = 0.0, mu = 0.0), 4)
   ("dambreak_coulomb_slope", 1, (dambreak..., hr = 0.0, R = 0.0, s = 0.5, mu = 0.3), 4)
   ("dambreak_coulomb_slope_moving", 1, (dambreak..., u0 = 10.0, hr = 0.0, R = 0.0, s = 0.5, mu = 0.3), 4)
]

# One-line summaries of the cases, for 'julia runbench.jl help'.
const case_summaries = Dict(
   "thacker_planar" => "Thacker (1981): frictionless oscillation in a parabolic bowl, planar free surface. Smooth flow with moving shorelines.",
   "thacker_planar_2d" => "Thacker (1981): the lake goes round a paraboloid bowl, planar free surface. Smooth flow with moving shorelines.",
   "thacker_radial_2d" => "Thacker (1981): radially symmetric oscillation in a paraboloid bowl, curved free surface.",
   "dambreak_dry" => "Ritter (1892): dam break onto a dry flat bed, frictionless. A rarefaction and a sharp wet/dry front.",
   "dambreak_wet" => "Stoker (1957): dam break onto a 2 m layer, frictionless. A rarefaction and a shock.",
   "dambreak_coulomb_slope" => "Mangeney et al. (2000): dam break down a slope with Coulomb friction, released from rest. Includes the start-up error of the friction.",
   "dambreak_coulomb_slope_moving" => "As dambreak_coulomb_slope, but the reservoir starts moving at 10 m/s. The front and rarefaction, without the start-up error.",
)

const reference_dir = "references"
const worse_ratio = 1.1

# Write the input file for case at cells cells per tile, and return its name.
function refined_input(case, dim, cells)
   name = "$(case)_n$cells"
   open(test_inputfile(name), "w") do io
      for line in readlines(test_inputfile(case))
         if occursin(r"^\s*nXpertile\s*="i, line)
            println(io, "nXpertile = $cells")
         elseif dim == 2 && occursin(r"^\s*nYpertile\s*="i, line)
            println(io, "nYpertile = $cells")
         elseif occursin(r"^\s*directory\s*="i, line)
            println(io, "directory = $name/")
         elseif occursin(r"^\s*Initial condition\s*="i, line)
            println(io, "Initial condition = $(name)_ic/000000")
         else
            println(io, line)
         end
      end
   end
   return name
end

# The reference to compare with, from the --reference argument spec: returns
# the file holding its errors (or nothing for 'none') and a description of it.
# A git ref is resolved to its commit, whose errors are kept in
# references/<hash>.txt and recorded by record_reference.sh if that is missing
# or lacks any of the runs in needed (a list of (case, cells per tile)), as when
# cases have been added since.
function resolve_reference(spec, needed)
   spec == "none" && return nothing, "none"
   if isfile(spec)
      return spec, reference_header(spec)
   end
   hash = try
      strip(read(pipeline(`git rev-parse --short --verify --quiet "$spec^{commit}"`), String))
   catch
      ""
   end
   if isempty(hash)
      println("--reference $spec is neither a file nor a git ref.")
      exit(1)
   end
   file = joinpath(reference_dir, "$hash.txt")
   missing_runs = isfile(file) ? setdiff(needed, keys(read_reference(file))) : needed
   if !isempty(missing_runs)
      println(isfile(file) ? "The reference for $spec ($hash) lacks $(join(unique(first.(missing_runs)), ", ")): re-recording it with record_reference.sh." :
              "No reference recorded for $spec ($hash): recording it with record_reference.sh.")
      println("This builds $hash and runs every benchmark with it, and takes several minutes.\n")
      run(`./record_reference.sh $spec`)
      isfile(file) || (println("record_reference.sh did not write $file."); exit(1))
      println()
   end
   return file, reference_header(file)
end

# Description of a reference file from its header: the commit it was recorded
# at, and any note on how (such as a change of drag).
function reference_header(file)
   desc = "unknown commit"; notes = String[]
   for line in readlines(file)
      startswith(line, "#") || break
      m = match(r"at commit (.+?)\.?$", line)
      if m !== nothing
         desc = "commit " * m.captures[1]
      elseif startswith(line, "# Note: ")
         push!(notes, line[9:end])
      end
   end
   return join([desc; notes], ", ")
end

# Reference errors in file, as a Dict (case, cells) => (dx, depth error,
# discharge error); empty for no file.
function read_reference(file)
   ref = Dict{Tuple{String,Int},Tuple{Float64,Float64,Float64}}()
   file === nothing && return ref
   for line in readlines(file)
      if startswith(line, "#") || isempty(strip(line))
         continue
      end
      f = strip.(split(line, ','))
      ref[(f[1], parse(Int, f[2]))] = (parse(Float64, f[3]), parse(Float64, f[4]),
                                       parse(Float64, f[5]))
   end
   return ref
end

function git_revision()
   try
      rev = strip(read(`git rev-parse --short HEAD`, String))
      dirty = !isempty(strip(read(`git status --porcelain --untracked-files=no -- ../src`, String)))
      return dirty ? "$rev (with uncommitted changes to src/)" : rev
   catch
      return "unknown"
   end
end

# Write the errors in errs (as read_reference returns) to file.
function write_reference(file, errs)
   open(file, "w") do io
      println(io, "# Errors of the benchmarks in runbench.jl, recorded by")
      println(io, "# 'julia runbench.jl --save' at commit $(git_revision()).")
      println(io, "# case, cells per tile, dx (m), depth error, discharge error")
      for k in sort(collect(keys(errs)))
         dx, eh, eq = errs[k]
         @printf(io, "%s, %d, %g, %.4e, %.4e\n", k[1], k[2], dx, eh, eq)
      end
   end
end

function parse_args(args)
   names = String[]; levels = nothing; refspec = "main"; save = nothing; markdown = nothing
   i = 1
   while i <= length(args)
      a = args[i]
      if a == "--levels"
         levels = parse(Int, args[i+1]); i += 1
      elseif a == "--reference"
         refspec = args[i+1]; i += 1
      elseif a == "--save"
         save = args[i+1]; i += 1
      elseif a == "--markdown"
         markdown = args[i+1]; i += 1
      elseif a in ("--help", "-h")
         push!(names, "help")
      elseif startswith(a, "--")
         error("Unknown option $a")
      else
         push!(names, a)
      end
      i += 1
   end
   isempty(names) && push!(names, "all")
   if "help" in names
      print_help()
      exit()
   end
   selected = String[]
   for (case, dim, _, _) in cases
      if "all" in names || "$(dim)d" in names || case in names
         push!(selected, case)
      end
   end
   bad = setdiff(names, ["all", "1d", "2d"], [c[1] for c in cases])
   if !isempty(bad) || isempty(selected)
      println(isempty(bad) ? "No valid benchmarks selected." :
              "Unknown selection $(join(bad, ", ")).")
      print_options()
      exit(1)
   end
   return selected, levels, refspec, save, markdown
end

# Print a list of the available selections.
function print_options()
   println("Options are 'all', '1d', '2d' and the cases ",
           join(["'$(c[1])'" for c in cases], ", "), ".")
   println("'help' describes the cases. See the header of runbench.jl and README.md for more.")
end

# Describe the options and the cases, with the grids each runs by default.
function print_help()
   println("""
   usage: julia runbench.jl [selection ...] [--levels N] [--reference REF] [--save FILE] [--markdown FILE]

     selection          'all' (the default), '1d', '2d', or case names (below)
     --levels N         run N grids per case, halving dx each time (default: the case's own number)
     --reference REF    compare with a git ref (default: main), a file from --save, or 'none'
     --save FILE        write this run's errors to FILE
     --markdown FILE    also write the summary table to FILE

   Cases (grid spacings at their default levels):""")
   for (case, dim, _, nlev) in cases
      lines = readlines(test_inputfile(case))
      dx0 = parse(Float64, input_setting(lines, "Xtilesize")) / parse(Int, input_setting(lines, "nXpertile"))
      grids = join([@sprintf("%g", dx0 / 2^l) for l = 0:nlev-1], ", ")
      println("\n  $case  ($(dim)D; dx = $grids m)")
      println("    ", get(case_summaries, case, ""))
   end
   println("\nSee README.md for the solutions, the error measures and the reference.")
end

order(e0, e1) = (isnan(e0) || isnan(e1) || e1 <= 0) ? NaN : log2(e0 / e1)
fmt(x, f) = isnan(x) ? "-" : @sprintf("%s", f(x))

# Status of a run, and the colour it is printed in: FAIL if the run failed or
# produced a negative depth, WORSE if an error is more than worse_ratio times
# its reference, BETTER if both are less than 1 / worse_ratio times (a prompt
# to re-record the reference), NEW if there is no reference, and OK otherwise.
function run_status(err, r, rh, rq)
   err && return "FAIL"
   r === nothing && return "NEW"
   (rh > worse_ratio || rq > worse_ratio) && return "WORSE"
   (rh < 1 / worse_ratio && rq < 1 / worse_ratio) && return "BETTER"
   return "OK"
end
const status_style = Dict("FAIL" => (color = :red, bold = true),
                          "WORSE" => (color = :red, bold = true),
                          "BETTER" => (color = :cyan, bold = false),
                          "NEW" => (color = :yellow, bold = false),
                          "OK" => (color = :green, bold = false))

function main(args)
   if !isfile("runbench.jl")
      println("Run runbench.jl from the benchmarks directory.")
      exit(1)
   end
   names, levels, refspec, save, markdown = parse_args(args)
   needed = Tuple{String,Int}[]
   for (case, dim, sol, nlev) in cases
      case in names || continue
      cells0 = parse(Int, input_setting(readlines(test_inputfile(case)), "nXpertile"))
      append!(needed, [(case, cells0 * 2^lev) for lev = 0:(something(levels, nlev) - 1)])
   end
   reffile, refdesc = resolve_reference(refspec, needed)
   ref = read_reference(reffile)
   println("Reference: $refspec, $refdesc; its errors on the same grid as each run.\n")
   rows = []
   errs = Dict{Tuple{String,Int},Tuple{Float64,Float64,Float64}}()

   for (case, dim, sol, nlev) in cases
      case in names || continue
      lines = readlines(test_inputfile(case))
      cells0 = parse(Int, input_setting(lines, "nXpertile"))
      tilesize = parse(Float64, input_setting(lines, "Xtilesize"))
      nrun = something(levels, nlev)
      println(isempty(rows) ? "" : "\n", "Running $case ($(dim)D, $nrun grid$(nrun == 1 ? "" : "s"))...")
      eh0 = NaN; eq0 = NaN
      for lev = 0:(nrun - 1)
         cells = cells0 * 2^lev
         dx = tilesize / cells
         name = refined_input(case, dim, cells)
         println("Run #$(lev + 1): ($name, dx = $(@sprintf("%g", dx)) m)")
         # Only the simulation is timed: the first call of each Julia function
         # compiles it, which would make the first run look the slowest.
         err = make_exact_ic(name, dim, sol)
         t0 = time()
         err = err || run_simulation(name)
         wall = time() - t0
         eh = NaN; eq = NaN
         if !err
            err = check_positivity(name, dim) != 0
            eh, eq = exact_errors(name, dim, sol)
         end
         rm(test_inputfile(name), force = true)
         r = get(ref, (case, cells), nothing)
         rh = (r === nothing) ? NaN : eh / r[2]
         rq = (r === nothing) ? NaN : eq / r[3]
         status = run_status(err, r, rh, rq)
         if !err
            oh = order(eh0, eh); oq = order(eq0, eq)
            for (what, e, o, k) in (("Depth", eh, oh, 2), ("Discharge", eq, oq, 3))
               @printf("%s error %.3e, order %s; reference %s", what, e,
                       fmt(o, x -> @sprintf("%.2f", x)),
                       r === nothing ? "none\n" : @sprintf("%.3e, ratio %.4f\n", r[k], e / r[k]))
            end
            @printf("Simulation time %.1f s\n", wall)
         end
         printstyled(status, "\n"; status_style[status]...)
         push!(rows, (case, dx, eh, order(eh0, eh), eq, order(eq0, eq), rh, rq, wall, status))
         err || (errs[(case, cells)] = (dx, eh, eq))
         eh0, eq0 = eh, eq
      end
   end

   hdr = ["case", "dx (m)", "depth error", "order", "discharge error", "order",
          "depth / ref", "discharge / ref", "sim time (s)", "status"]
   rowcells(r) = [r[1], @sprintf("%g", r[2]), fmt(r[3], x -> @sprintf("%.3e", x)),
               fmt(r[4], x -> @sprintf("%.2f", x)), fmt(r[5], x -> @sprintf("%.3e", x)),
               fmt(r[6], x -> @sprintf("%.2f", x)), fmt(r[7], x -> @sprintf("%.4f", x)),
               fmt(r[8], x -> @sprintf("%.4f", x)), @sprintf("%.1f", r[9]), r[10]]
   table = vcat([hdr], [rowcells(r) for r in rows])
   width = [maximum(length(t[j]) for t in table) for j in eachindex(hdr)]
   println("\nSummary (ratios to the reference $refspec):")
   for t in table
      println(join([lpad(t[j], width[j]) for j in eachindex(t)], "  "))
   end
   if markdown !== nothing
      open(markdown, "w") do io
         println(io, "Kestrel benchmarks at commit $(git_revision()). Errors are relative 1-norms,")
         println(io, "maximised over output times; ratios are to the errors on the same grid of the reference $refspec ($refdesc).\n")
         println(io, "| ", join(hdr, " | "), " |")
         println(io, "|", repeat("---|", length(hdr)))
         for t in table[2:end]
            println(io, "| ", join(t, " | "), " |")
         end
      end
   end
   if save !== nothing
      write_reference(save, errs)
      println("\nWrote the errors to $save.")
   end
   nbad = count(r -> r[10] in ("FAIL", "WORSE"), rows)
   n = length(rows)
   @printf("\nResults: %i/%i within %.0f%% of the reference or better%s (%.1f%%)\n",
           n - nbad, n, 100 * (worse_ratio - 1), nbad == 0 ? "!" : "", 100 * (n - nbad) / n)
   for st in ("FAIL", "WORSE", "BETTER", "NEW")
      k = count(r -> r[10] == st, rows)
      k > 0 && printstyled("$k $st\n"; status_style[st]...)
   end
   exit(nbad == 0 ? 0 : 1)
end

main(ARGS)
