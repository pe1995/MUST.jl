#= M3D spectrum synthesis for dispatch models =#
using Pkg; Pkg.activate(".")
using MUST 
using TSO
using ArgParse

struct ChemicalComposition
    abundances ::Dict
end

# argument parsing for arrays
ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:Number} = parse.(A, split(x, ",", keepempty=false))
ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:AbstractString} = split(x, ",", keepempty=false)
ArgParse.parse_item(::Type{ChemicalComposition}, x::AbstractString) = begin
    elements_abundances = split(x, ",", keepempty=false)
    elements = [Symbol(split(a, "_", keepempty=false) |> first |> string) for a in elements_abundances]
    abundances = [parse(Float64, split(a, "_", keepempty=false) |> last) for a in elements_abundances]

    ChemicalComposition(Dict(e=>a for (e, a) in zip(elements, abundances)))
end

# available command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--run", "-r"
        help = "DISPATCH run to compute spectra for."
        arg_type = String
        required=true
    "--snapshot", "-n"
        help = "Number of the snapshot within `run`. Converts if it not converted yet."
        arg_type = Int
        required=true
    "--multi_threads", "-t"
        help = "Number of threads to use for M3D."
        arg_type = Int
        default=20
    "--lambda_start", "-s"
        help = "Wavelength window start wavelength [Å]."
        arg_type = Float64
        default = 6100.0
    "--lambda_end", "-e"
        help = "Wavelength window end wavelength [Å]."
        arg_type = Float64
        default = 6790.0
    "--lambda_step", "-d"
        help = "Wavelength window end wavelength [Å]."
        arg_type = Float64
        default = 0.5
    "--linelists"
        help = "List of linelists. Choose `all` to use the default.
        Also available: CHonly (turn molecules except CH off)"
        arg_type = String
        default="all"
    "--feh"
        help = "Metallicity [Fe/H]."
        arg_type = Float64
        default = 0.0
    "--alpha"
        help = "Alpha enhancement [alpha/Fe]."
        arg_type = Float64
        default = 0.0
    "--composition"
        help = "Chemical composition [X/Fe]. 
        List different elements (case insensitive) with underscore e.g.: C_3.0,N_0.4,O_0.5"
        arg_type = ChemicalComposition
        default = ChemicalComposition(Dict())
    "--horizontal_resolution", "-x"
        help = "Resample cube to this resolution in the horizontal."
        arg_type = Int
        default=5
    "--vertical_resolution", "-z"
        help = "Resample cube to this resolution in the vertical. Set `-1` to use  `size(box, 3) * 2 - 1`."
        arg_type = Int
        default=-1
    "--atom"
        help = "Atom to use."
        arg_type = String
        default=""
    "--datafolder"
        help = "Dispatch data directory."
        arg_type = String
        default="data"
    "--absolutepath"
        help = "Snapshot data directory in `datafolder` is given as absolute path, instead of relative."
        action = :store_true
    "--dims"
        help = "Number of vertical atmosphere splits."
        arg_type = Int
        default=32
    "--nnu"
        help = "Number of frequency splits."
        arg_type = Int
        default=32
    "--save_resolved"
        help = "Save resolved spectrum."
        action = :store_true
    "--long_scheme"
        help = "Scheme for long characteristics."
        arg_type = String
        default="lobatto"
    "--short_scheme"
        help = "Scheme for short characteristics."
        arg_type = String
        default="radau"
    "--absdat"
        help = "Absdat file."
        arg_type = String
        default="./input_multi3d/TS_absdat.dat"
    "--abund"
        help = "Default abund file. Will be modified by `--composition`."
        arg_type = String
        default="./input_multi3d/abund_magg"
    "--NLTE"
        help = "Turn on NLTE."
        action = :store_true
    "--move"
        help = "Move results to DISPATCH after computation has finished."
        action = :store_true
    "--extension"
        help = "Additional, optional name extension. No '_' added automatically."
        default = ""
        arg_type = String
end

# code setup
begin
    # read and parse command line arguments
    arguments = parse_args(ARGS, s)

    # import dispatch and M3D
    MUST.@import_m3dis "../../../Multi3D"
    MUST.@import_dispatch "../../../dispatch2"
end

# simulation details
begin
    # wavelengh coverage
	λs = arguments["lambda_start"]
	λe = arguments["lambda_end"]
	Δλ = arguments["lambda_step"]
    nλ = (λe - λs) / Δλ
    window = MUST.@sprintf "lam_%i-%i" λs λe

    # simulation
    name = arguments["run"]
    if arguments["absolutepath"]
        run = joinpath("$(arguments["datafolder"])/$(name)")
    else
        run = @in_dispatch("$(arguments["datafolder"])/$(name)")
    end
    isnap = arguments["snapshot"]

    extension = arguments["extension"]
    @info "Spectrum synthesis for model $(name) ($(extension)) in spectral window $(window) Å."
end

# convert snapshot if not already done, convert to M3D aswell (full res.)
begin
	snaps_ready = MUST.list_snapshots(MUST.converted_snapshots(run))

	if !(isnap in snaps_ready)
		@info "Converting snapshot $(name), Number: $(isnap)"
		snapshotBox(isnap, folder=run, to_multi=true, legacy=false)		
		@info "Snapshot $(name), Number: $(isnap) converted."
	end

	b, bτ = pick_snapshot(run, isnap)
end

# chemical details
begin
    linelists = if arguments["linelists"] == "all" 
        String[
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list",
            "input_multi3d/LINE-LISTS/25500-200000_cut-4/atom_25500-200000.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/Hlinedata",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/12CH_multi.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/13CH_multi.list",
		    #"input_multi3d/LINE-LISTS/combined_molecules/H2O_multi.list",
            MUST.relative_path.(@in_m3dis("."), MUST.glob("*.list", @in_m3dis("input_multi3d/LINE-LISTS/combined_molecules/most_relevant/")))...
        ]
    elseif arguments["linelists"] == "CHonly" 
        String[
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/vald_2490-25540.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/1000-2490-vald.list",
            "input_multi3d/LINE-LISTS/25500-200000_cut-4/atom_25500-200000.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/12CH_multi.list",
            "input_multi3d/LINE-LISTS/ADDITIONAL-LISTS/13CH_multi.list",
        ]
    else
        arguments["linelists"]
    end
    maskL = isfile.([@in_m3dis(l) for l in linelists])
    if any(.!maskL)
        @warn "Linelists could not be found at: $(linelists[.!maskL])"
    end
    linelists = linelists[maskL]

    FeH = arguments["feh"]
    α = arguments["alpha"]
    composition = arguments["composition"]
	abund_file = MUST.abund_abundances(;
		α=α, 
		composition.abundances...,
		default=arguments["abund"]
	)
    
    cs = join(["[$(k)/Fe]=$(v)" for (k, v) in composition.abundances], ",")
    cstring = "[Fe/H]=$(FeH), [α/Fe]=$(α), with: $(cs)"
    @info "Chemical composition: $(cstring)"
end

# run M3D
begin
    nz = if arguments["vertical_resolution"] == -1
        size(b, 3) * 2 - 1
    else
        arguments["vertical_resolution"]
    end

    spectrum_namelist = Dict(
        :model_folder=>run,
        :linelist=>nothing,
        :absmet=>nothing,
        :linelist_params=>(:line_lists=>linelists,),
        :atom_params=>(:atom_file=>arguments["atom"],),
        :spectrum_params=>(:daa=>Δλ, :aa_blue=>λs, :aa_red=>λe, :in_log=>false),
        :atmos_params=>(
            :dims=>arguments["dims"], 
            :atmos_format=>"must",
            :use_density=>true, 
            :use_ne=>false,
            :FeH=>FeH,
            :nx=>arguments["horizontal_resolution"],
            :ny=>arguments["horizontal_resolution"],
            :nz=>nz
        ),
        :m3d_params=>(
            :save_resolved=>arguments["save_resolved"],
            :long_scheme=>arguments["long_scheme"],
            :short_scheme=>arguments["short_scheme"],
            :n_nu=>arguments["nnu"],
        ),
        :composition_params=>(
            :absdat_file=>arguments["absdat"],
            :abund_file=>abund_file
        )
    )

    m3dis_kwargs = Dict(
        :threads=>arguments["multi_threads"]
    )

    result = MUST.spectrum(
        "m3dis_$(isnap)"; 
        name=name*"_"*window*extension, 
        NLTE=arguments["NLTE"], 
        slurm=false, 
        namelist_kwargs=spectrum_namelist,
        m3dis_kwargs=m3dis_kwargs
    )

    newp = @in_m3dis("data/m3dis_$(isnap)_$(name)_$(window)$(extension)")
    if arguments["move"]
        mvp = joinpath(run, "spectra_sn$(isnap)_$(window)$(extension)")
        mv(newp, mvp, force=true)
        @info "M3D run saved at $(mvp)."
    else
        @info "M3D run saved at $(newp)."
    end
end

