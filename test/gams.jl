using GAMSFiles, QuadDIRECT, ProgressMeter, Compat

struct ProblemInfo
    smoothness::String
    convexity::String
    nvars::Int
    minval::Float64
end

gamsdir = joinpath(@__DIR__, "gams")
if !isdir(gamsdir)
    mkdir(gamsdir)
end
# Download the problems studied by L. M. Rios and N. V. Sahinidis, "Derivative-free optimization:
# a review of algorithms and comparison of software implementations."
# Journal of Global Optimization,  July 2013, Volume 56, Issue 3, pp 1247–1293.
site = "http://thales.cheme.cmu.edu/dfo/comparison"
gamsdirs = ("convex_nonsmooth_gams",
            "convex_smooth_gams",
            "nonconvex_nonsmooth_gams",
            "nonconvex_smooth_gams")
cd(gamsdir) do
    for dir in gamsdirs
        subdir = joinpath(gamsdir, dir)
        if !isdir(subdir)
            pathfile = subdir*".zip"
            download(joinpath(site, dir*".zip"), pathfile)
            mkdir(subdir)
            run(`unzip $pathfile -d $subdir`)
        end
    end
end
problemdatafile = joinpath(gamsdir, "problemdata.zip")
if !isfile(problemdatafile)
    download(joinpath(site, "problemdata.zip"), problemdatafile)
    cd(gamsdir) do
        mkdir("problemdata")
        run(`unzip $problemdatafile -d problemdata`)
    end
end
answerfile = joinpath(gamsdir, "answers.tsv")
if !isfile(answerfile)
    error("must parse the answer file using tips from https://eagereyes.org/data/scrape-tables-using-google-docs using URL ", joinpath(site, "comp.html"))
end
answers_array = readdlm(answerfile, '\t')
@assert(answers_array[1,2] == "problem")
@assert(answers_array[1,3] == "smoothness")
@assert(answers_array[1,4] == "convexity")
@assert(answers_array[1,5] == "variables")
@assert(answers_array[1,7] == "solution")
answers = Dict{String,ProblemInfo}()
for i = 2:size(answers_array, 1)
    prob = answers_array[i,2]
    info = ProblemInfo(answers_array[i,3], answers_array[i,4], answers_array[i,5], answers_array[i,7])
    answers[prob] = info
end

function read_problemdata(file)
    lower, upper, x0 = open(file) do io
        values = split(readstring(io))
        k = 1
        nvars = parse(Int, values[1])
        lower = Vector{Float64}(uninitialized, nvars)
        upper = similar(lower)
        x0 = similar(lower)
        for i = 1:nvars
            lower[i] = parse(Float64, values[k+=1])
        end
        for i = 1:nvars
            upper[i] = parse(Float64, values[k+=1])
        end
        for i = 1:nvars
            x0[i] = parse(Float64, values[k+=1])
        end
        @assert(k == length(values))
        return lower, upper, x0
    end
    return lower, upper, x0
end

stopping_criterion(v) = max(1.01*v, v+0.01)
stopping_criterion(info::ProblemInfo) = stopping_criterion(info.minval)

# Debugging utility: run a single GAMS file
function gams1(filename::AbstractString; kwargs...)
    dirname, fname = splitdir(filename)
    basename, ext = splitext(fname)
    lower, upper, x0 = read_problemdata(joinpath(gamsdir, "problemdata", basename)*".problem.data")
    gams = parsegams(filename)
    modex, axs = parsegams(Module, basename, gams)
    mod = eval(modex)
    f = getfield(mod, :objective)
    fwrap = LoggedFunction(f)
    root, x0 = @eval analyze($fwrap, $x0, $lower, $upper; rtol=0, fvalue=stopping_criterion($(answers[basename])), $(kwargs...))
    return root, x0, fwrap
end

failures = String[]
successes = Dict{String,Any}()

# Because of https://github.com/JuliaLang/julia/issues/25927
excludes = ["nonconvex_smooth_gams/nnls.gms",
            "nonconvex_smooth_gams/chebyqad.gms"]

mode = "optimize"   # any choice other than "optimize" tests that the objective can be evaluated

terminating = false
for dir in gamsdirs
    terminating && break
    subdir = joinpath(gamsdir, dir)
    cd(subdir) do
        files = readdir()
        @showprogress 1 "Working on $subdir" for file in files
            endswith(file, ".gms") || (println("skipping ", file); continue)
            fullfile = joinpath(subdir, file)
            exclude = false
            for excl in excludes
                if endswith(fullfile, excl)
                    exclude = true
                    println("excluding ", fullfile)
                end
            end
            exclude && continue
            basename, ext = splitext(file)
            lower, upper, x0 = read_problemdata(joinpath(gamsdir, "problemdata", basename)*".problem.data")
            local x0g
            try
                modex, axs = parsegams(Module, file)
                x0g = axs[1]  # initialization provided by the file
                if isdefined(Main, Symbol(basename))
                    mod = getfield(Main, Symbol(basename))
                else
                    mod = eval(modex)
                end
                f = getfield(mod, :objective)
                # # FIXME: we need a better way of setting the splits
                if mode == "optimize"
                    fwrap = LoggedFunction(f)
                    root, x0 = @eval analyze($fwrap, $x0, $lower, $upper; rtol=0, fvalue=stopping_criterion($(answers[basename])))
                    successes[file] = fwrap.values
                else
                    val = Base.invokelatest(f, x0)
                    if !isfinite(val)
                        # Also try any initialization present in the GAMS file
                        val = Base.invokelatest(f, x0g)
                        if isfinite(val)
                            println(fullfile, ": finite result from file initialization but not problemdata")
                        else
                            println("\nnot finite on ", fullfile)
                            println()
                            push!(failures, fullfile)
                        end
                    end
                    if isfinite(val)
                        successes[fullfile] = val
                    end
                end
            catch ex
                if ex isa InterruptException
                    terminating = true
                    break
                end
                push!(failures, fullfile)
                println("\nfailed on ", fullfile)
                println("problemdata axes ", Compat.axes(x0), ", and from the file it's ", Compat.axes(x0g))
                showerror(STDOUT, ex)
                println()
            end
        end
    end
end
