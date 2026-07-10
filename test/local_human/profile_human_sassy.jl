#!/usr/bin/env julia

using Pkg

const ROOT_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const PROFILE_ENV = get(ENV, "CHOPOFF_PROFILE_ENV", "/home/rstudio/livemount/kornel_dev/temp_upload/profiletools")
const PROFILE_MODE = lowercase(strip(get(ENV, "CHOPOFF_PROFILE_MODE", "baseline")))
const NEED_PPROF = PROFILE_MODE in ("cpu", "allocs", "all")
const NEED_STATHTML = PROFILE_MODE in ("cpu", "all")

if NEED_PPROF || NEED_STATHTML
    Pkg.activate(PROFILE_ENV; io = devnull)
    NEED_PPROF && @eval using PProf
    NEED_STATHTML && @eval using StatProfilerHTML
end

Pkg.activate(ROOT_DIR; io = devnull)
using CHOPOFF
using BioSequences
using Dates
using Profile
using Printf

const DEFAULT_GENOME = "/home/rstudio/livemount/Bio_data/references/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
const DEFAULT_GUIDES = joinpath(@__DIR__, "data", "guides_for_tests.txt")
const SAMPLE_GENOME = joinpath(ROOT_DIR, "test", "sample_data", "genome", "semirandom.fa")
const SAMPLE_GUIDES = joinpath(ROOT_DIR, "test", "sample_data", "guides.txt")

function parse_bool_env(name::String, default::Bool)
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    raw in ("1", "true", "yes", "y") && return true
    raw in ("0", "false", "no", "n") && return false
    error("Invalid boolean in environment variable $name: '$raw'")
end

function parse_int_env(name::String, default::Int)
    raw = strip(get(ENV, name, string(default)))
    isempty(raw) && return default
    return parse(Int, raw)
end

function parse_float_env(name::String, default::Float64)
    raw = strip(get(ENV, name, string(default)))
    isempty(raw) && return default
    return parse(Float64, raw)
end

function parse_threads_env(name::String, default::String)
    raw = strip(get(ENV, name, default))
    return [parse(Int, strip(x)) for x in split(raw, ',') if !isempty(strip(x))]
end

function require_file(path::String)
    isfile(path) || error("Required file not found: $path")
    return path
end

function validate_genome(path::String)
    require_file(path)
    if !endswith(lowercase(path), ".2bit")
        require_file(path * ".fai")
    end
    return path
end

function load_guides(path::String, motif::Motif; limit::Int = 0)
    require_file(path)
    raw = String[]
    for line in eachline(path)
        s = uppercase(strip(line))
        isempty(s) && continue
        startswith(s, "#") && continue
        push!(raw, s)
        limit > 0 && length(raw) >= limit && break
    end
    isempty(raw) && error("No guides found in $path")

    expected_len = length_noPAM(motif)
    bad = findall(!=(expected_len), length.(raw))
    if !isempty(bad)
        first_bad = bad[1]
        error(
            "Guide length mismatch for motif $(motif.alias): expected $expected_len, " *
            "got $(length(raw[first_bad])) at guide $first_bad ($(raw[first_bad]))",
        )
    end
    return LongDNA{4}.(raw)
end

function count_result_rows(path::String)
    isfile(path) || return 0
    rows = 0
    for _ in eachline(path)
        rows += 1
    end
    return max(rows - 1, 0)
end

function write_kv(path::String, pairs)
    open(path, "w") do io
        for (k, v) in pairs
            println(io, k, ",", v)
        end
    end
end

function warmup_compile(distance::Int, use_avx512::Bool, force_safe_minima::Bool, algorithm::Symbol = :sassy)
    isfile(SAMPLE_GENOME) || return
    isfile(SAMPLE_GENOME * ".fai") || return
    isfile(SAMPLE_GUIDES) || return

    motif = Motif("Cas9"; distance = min(distance, 3))
    guides = load_guides(SAMPLE_GUIDES, motif; limit = 1)
    tmp = mktempdir(prefix = "chopoff_profile_warmup_")
    try
        if algorithm == :sassy
            search_sassy(
                guides,
                SAMPLE_GENOME,
                motif,
                joinpath(tmp, "warmup.csv");
                distance = min(distance, 3),
                early_stopping = fill(10, min(distance, 3) + 1),
                use_avx512 = use_avx512,
                force_safe_minima = force_safe_minima,
            )
        elseif algorithm == :prefixHashScan
            CHOPOFF.search_prefixHashScan(
                guides,
                SAMPLE_GENOME,
                motif,
                joinpath(tmp, "warmup.csv");
                distance = min(distance, 3),
                early_stopping = fill(10, min(distance, 3) + 1),
                query_variant = :auto,
            )
        end
    finally
        rm(tmp; recursive = true, force = true)
    end
end

Base.@kwdef struct ProfileConfig
    genome::String
    guides_path::String
    motif_name::String
    distance::Int
    guide_limit::Int
    use_avx512::Bool
    force_safe_minima::Bool
    algorithm::Symbol
    output_parent::String
    run_label::String
    sample_rate::Float64
    profile_delay::Float64
    profile_n::Int
end

function load_config(mode::String)
    genome = validate_genome(String(strip(get(ENV, "CHOPOFF_PROFILE_GENOME", DEFAULT_GENOME))))
    guides_path = String(strip(get(ENV, "CHOPOFF_PROFILE_GUIDES", DEFAULT_GUIDES)))
    motif_name = String(strip(get(ENV, "CHOPOFF_PROFILE_MOTIF", "Cas9")))
    distance = parse_int_env("CHOPOFF_PROFILE_DISTANCE", 3)
    guide_limit = parse_int_env("CHOPOFF_PROFILE_GUIDE_LIMIT", 0)
    use_avx512 = parse_bool_env("CHOPOFF_PROFILE_USE_AVX512", true)
    force_safe_minima = parse_bool_env("CHOPOFF_PROFILE_FORCE_SAFE_MINIMA", false)
    algorithm = Symbol(strip(get(ENV, "CHOPOFF_PROFILE_ALGORITHM", "sassy")))
    algorithm in (:sassy, :prefixHashScan) || error("CHOPOFF_PROFILE_ALGORITHM must be sassy or prefixHashScan")
    sample_rate = parse_float_env("CHOPOFF_PROFILE_ALLOC_SAMPLE_RATE", 0.01)
    profile_delay = parse_float_env("CHOPOFF_PROFILE_DELAY", 0.001)
    profile_n = parse_int_env("CHOPOFF_PROFILE_BUFFER", 10_000_000)

    stamp = Dates.format(now(Dates.UTC), dateformat"yyyymmdd_HHMMSS")
    output_parent = String(strip(get(ENV, "CHOPOFF_PROFILE_OUTPUT_PARENT", joinpath(@__DIR__, "outputs", "profile_" * stamp))))
    run_label = String(strip(get(ENV, "CHOPOFF_PROFILE_RUN_LABEL", mode)))

    return ProfileConfig(
        genome = genome,
        guides_path = guides_path,
        motif_name = motif_name,
        distance = distance,
        guide_limit = guide_limit,
        use_avx512 = use_avx512,
        force_safe_minima = force_safe_minima,
        algorithm = algorithm,
        output_parent = output_parent,
        run_label = run_label,
        sample_rate = sample_rate,
        profile_delay = profile_delay,
        profile_n = profile_n,
    )
end

function prepare_run(cfg::ProfileConfig, mode::String)
    run_dir = joinpath(cfg.output_parent, cfg.run_label)
    mkpath(run_dir)
    mkpath(joinpath(run_dir, "profiles"))

    motif = Motif(cfg.motif_name; distance = cfg.distance)
    guides = load_guides(cfg.guides_path, motif; limit = cfg.guide_limit)

    println("CHOPOFF human profile")
    println("mode: ", mode)
    println("time_utc: ", Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println("run_dir: ", run_dir)
    println("genome: ", cfg.genome)
    println("guides: ", cfg.guides_path)
    println("guide_count: ", length(guides))
    println("motif: ", cfg.motif_name)
    println("distance: ", cfg.distance)
    println("threads: ", Threads.nthreads())
    println("use_avx512: ", cfg.use_avx512)
    println("force_safe_minima: ", cfg.force_safe_minima)
    println("algorithm: ", cfg.algorithm)

    return run_dir, motif, guides
end

function run_search(cfg::ProfileConfig, motif::Motif, guides::Vector{LongDNA{4}}, out_csv::String)
    if cfg.algorithm == :sassy
        search_sassy(
            guides,
            cfg.genome,
            motif,
            out_csv;
            distance = cfg.distance,
            early_stopping = fill(1_000_000, cfg.distance + 1),
            use_avx512 = cfg.use_avx512,
            force_safe_minima = cfg.force_safe_minima,
        )
    elseif cfg.algorithm == :prefixHashScan
        CHOPOFF.search_prefixHashScan(
            guides,
            cfg.genome,
            motif,
            out_csv;
            distance = cfg.distance,
            early_stopping = fill(1_000_000, cfg.distance + 1),
            query_variant = :auto,
        )
    end
    return count_result_rows(out_csv)
end

function write_run_summary(path::String, cfg::ProfileConfig, mode::String, elapsed::Float64, rows::Int)
    write_kv(path, [
        "timestamp" => Dates.format(now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"),
        "mode" => mode,
        "genome" => cfg.genome,
        "guides" => cfg.guides_path,
        "guide_limit" => cfg.guide_limit,
        "motif" => cfg.motif_name,
        "distance" => cfg.distance,
        "threads" => Threads.nthreads(),
        "use_avx512" => cfg.use_avx512,
        "force_safe_minima" => cfg.force_safe_minima,
        "algorithm" => cfg.algorithm,
        "elapsed_s" => elapsed,
        "rows" => rows,
    ])
end

function run_baseline(cfg::ProfileConfig; label::String = "baseline")
    run_dir, motif, guides = prepare_run(cfg, label)
    warmup_compile(cfg.distance, cfg.use_avx512, cfg.force_safe_minima, cfg.algorithm)
    out_csv = joinpath(run_dir, "$(cfg.algorithm).csv")
    rows_ref = Ref(0)
    elapsed = @elapsed rows_ref[] = run_search(cfg, motif, guides, out_csv)
    write_run_summary(joinpath(run_dir, "summary.csv"), cfg, label, elapsed, rows_ref[])
    @printf("baseline elapsed: %.3fs | rows=%d\n", elapsed, rows_ref[])
    println("summary: ", joinpath(run_dir, "summary.csv"))
    return run_dir, elapsed, rows_ref[]
end

function save_pprof_top(pb_path::String, top_path::String)
    try
        run(pipeline(`$(PProf.pprof_jll.pprof()) -top $pb_path`, stdout = top_path))
    catch err
        @warn "Unable to write pprof top output" pb_path top_path exception = (err, catch_backtrace())
    end
end

function with_label(cfg::ProfileConfig, label::String)
    return ProfileConfig(
        genome = cfg.genome,
        guides_path = cfg.guides_path,
        motif_name = cfg.motif_name,
        distance = cfg.distance,
        guide_limit = cfg.guide_limit,
        use_avx512 = cfg.use_avx512,
        force_safe_minima = cfg.force_safe_minima,
        algorithm = cfg.algorithm,
        output_parent = cfg.output_parent,
        run_label = label,
        sample_rate = cfg.sample_rate,
        profile_delay = cfg.profile_delay,
        profile_n = cfg.profile_n,
    )
end

function run_cpu_profile(cfg::ProfileConfig)
    run_dir, motif, guides = prepare_run(cfg, "cpu")
    warmup_compile(cfg.distance, cfg.use_avx512, cfg.force_safe_minima, cfg.algorithm)
    out_csv = joinpath(run_dir, "$(cfg.algorithm).csv")
    rows_ref = Ref(0)

    Profile.init(n = cfg.profile_n, delay = cfg.profile_delay)
    Profile.clear()
    elapsed = @elapsed Profile.@profile rows_ref[] = run_search(cfg, motif, guides, out_csv)

    profiles_dir = joinpath(run_dir, "profiles")
    flat_path = joinpath(profiles_dir, "cpu_flat.txt")
    open(flat_path, "w") do io
        Profile.print(io; format = :flat, sortedby = :count, maxdepth = 80)
    end

    pb_path = PProf.pprof(out = joinpath(profiles_dir, "cpu.pb.gz"), web = false, from_c = false)
    save_pprof_top(pb_path, joinpath(profiles_dir, "cpu_top.txt"))
    statprofilehtml(path = joinpath(profiles_dir, "statprof"))

    write_run_summary(joinpath(run_dir, "summary.csv"), cfg, "cpu", elapsed, rows_ref[])
    @printf("cpu profile elapsed: %.3fs | rows=%d\n", elapsed, rows_ref[])
    println("pprof: ", pb_path)
    println("flat: ", flat_path)
    println("html: ", joinpath(profiles_dir, "statprof", "index.html"))
    return run_dir, elapsed, rows_ref[]
end

function run_alloc_profile(cfg::ProfileConfig)
    run_dir, motif, guides = prepare_run(cfg, "allocs")
    warmup_compile(cfg.distance, cfg.use_avx512, cfg.force_safe_minima, cfg.algorithm)
    out_csv = joinpath(run_dir, "$(cfg.algorithm).csv")
    rows_ref = Ref(0)

    Profile.Allocs.clear()
    elapsed = @elapsed Profile.Allocs.@profile sample_rate = cfg.sample_rate rows_ref[] = run_search(cfg, motif, guides, out_csv)

    profiles_dir = joinpath(run_dir, "profiles")
    pb_path = PProf.Allocs.pprof(out = joinpath(profiles_dir, "allocs.pb.gz"), web = false)
    save_pprof_top(pb_path, joinpath(profiles_dir, "allocs_top.txt"))

    write_run_summary(joinpath(run_dir, "summary.csv"), cfg, "allocs", elapsed, rows_ref[])
    @printf("alloc profile elapsed: %.3fs | rows=%d | sample_rate=%.4f\n", elapsed, rows_ref[], cfg.sample_rate)
    println("alloc pprof: ", pb_path)
    return run_dir, elapsed, rows_ref[]
end

function run_scaling_driver(cfg::ProfileConfig)
    threads = parse_threads_env("CHOPOFF_PROFILE_SCALING_THREADS", "1,2,4,8,16")
    parent = cfg.output_parent
    mkpath(parent)
    summary_path = joinpath(parent, "scaling_commands.txt")
    julia = get(ENV, "CHOPOFF_PROFILE_JULIA", joinpath(dirname(Sys.BINDIR), "bin", "julia"))
    if !isfile(julia)
        julia = joinpath(Sys.BINDIR, "julia")
    end

    open(summary_path, "w") do io
        println(io, "threads,command")
        for t in threads
            env = copy(ENV)
            env["JULIA_NUM_THREADS"] = string(t)
            env["CHOPOFF_PROFILE_MODE"] = "baseline"
            env["CHOPOFF_PROFILE_OUTPUT_PARENT"] = parent
            env["CHOPOFF_PROFILE_RUN_LABEL"] = "threads_$(t)"
            cmd = `$(julia) --project=$(ROOT_DIR) $(abspath(@__FILE__))`
            println(io, t, ",", join(cmd.exec, " "))
            println("running thread-scaling baseline with JULIA_NUM_THREADS=", t)
            run(setenv(cmd, env))
        end
    end
    println("scaling command log: ", summary_path)
end

function main()
    mode = PROFILE_MODE
    cfg = load_config(mode)

    if mode == "baseline"
        run_baseline(cfg)
    elseif mode == "cpu"
        run_cpu_profile(cfg)
    elseif mode == "allocs"
        run_alloc_profile(cfg)
    elseif mode == "all"
        run_baseline(with_label(cfg, "baseline"))
        run_cpu_profile(with_label(cfg, "cpu"))
        run_alloc_profile(with_label(cfg, "allocs"))
    elseif mode == "scaling"
        run_scaling_driver(cfg)
    else
        error("Unknown CHOPOFF_PROFILE_MODE=$mode. Use baseline, cpu, allocs, all, or scaling.")
    end
end

main()
