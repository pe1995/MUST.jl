using Pkg; Pkg.activate(".")
using MUST 
using TSO
using ArgParse

# argument parsing for arrays
ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:Number} = parse.(A, split(x, ",", keepempty=false))
ArgParse.parse_item(::Type{Vector{A}}, x::AbstractString) where {A<:AbstractString} = split(x, ",", keepempty=false)

# available command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--modelFolder"
        help = "Where to store new models + binned opacities."
        arg_type = String
        default = "interpolatedModels/M1"
    "--version"
        help = "Version to be added to the EoS table."
        default="v1.0"
        arg_type = String
    "--grid"
        help = "Grid in which to interpolate (Stagger, MARCS)."
        arg_type = String
        default = "Stagger"
    "--teff", "-t"
        help = "List of effective temperatures."
        arg_type = Vector{Float64}
        required = true
    "--logg", "-l"
        help = "List of log10 surface gravities."
        arg_type = Vector{Float64}
        required = true
    "--feh", "-f"
        help = "List of metallicities [Fe/H]."
        arg_type = Vector{Float64}
        required = true
    "--eos"
        help = "EoS to use. Can be set to `closest` to pick the closest EoS in the grid."
        arg_type = Vector{String}
        default = ["closest"]
    "--tau_bottom"
        help = "For interpolation: Bottom bounday in terms of optical depth."
        arg_type = Float64
        default = 6
    "--adiabatic_extrapolation"
        help = """Extrapolate adiabatically below the bottom boundary `tau_bottom` after the model has been averaged.
        The extrapolation is carried out until the geometrical bottom boundary as specified in the grid has been reached.
        """
        action = :store_true
    # General parameters
    "--patch_size"
        help = "Points per patch."
        arg_type = Float64
        default = 14
    "--tau_up"
        help = "Upper boundary (target) of simulation domain."
        arg_type = Float64
        default = -5.0
    "--tau_down"
        help = "Lower boundary (target) of simulation domain."
        arg_type = Float64
        default = 6.0
    "--tau_ee0"
        help = "Surface energy during Newton phase picked from this height."
        arg_type = Float64
        default = -2.5
    "--tau_zee0"
        help = "Location of the Newton cooling surface."
        arg_type = Float64
        default = -1.75
    "--tau_rho0"
        help = "Density normalization height."
        arg_type = Float64
        default = 0.0
    "--scale_resolution"
        help = "Down or upsampling of simulation domain."
        arg_type = Float64
        default = 0.75
    # Namelist specific parameters
    "--newton_time"
        help = "End time of Newton cooling."
        arg_type = Float64
        default = 100.0
    "--friction_time"
        help = "End time of artifical friction."
        arg_type = Float64
        default = 120.0
    "--newton_decay_scale"
        help = "Decay scale of Newtong cooling contribution after Newton time."
        arg_type = Float64
        default = 20.0
    "--friction_decay_scale"
        help = "Decay scale of friction after friction time."
        arg_type = Float64
        default = 20.0
    "--courant_hd"
        help = "Hydrodynamic courant factor."
        arg_type = Float64
        default = 0.25
    "--courant_rt"
        help = "Radiative transfer courant factor."
        arg_type = Float64
        default = 0.5
    "--out_time"
        help = "Spacing beetween snapshots."
        arg_type = Float64
        default = 1.0
    "--end_time"
        help = "Simulation duration."
        arg_type = Float64
        default = 400.0
    "--aux_params"
        help = "Aux data to save during dispatch run (`flux`, `qr`, `dt_rt`)."
        arg_type = Vector{String}
        default = ["flux"]
    "--upper_bc"
        help = "Upper boundary condition."
        arg_type = Int
        default = 7
    "--lower_bc"
        help = "Upper boundary condition."
        arg_type = Int
        default = 5
    "--radiation_bc"
        help = "Upper RT boundary condition."
        arg_type = Int
        default = 3
    "--smallr"
        help = "Lower density limit."
        arg_type = Float64
        default = 1.0e-9
    "--dt_small"
        help = "Lower timestep limit."
        arg_type = Float64
        default = 1.0e-10
    "--rt_freq"
        help = "Upper limit on RT timestep relative to the HD."
        arg_type = Float64
        default = 0.0
    "--cdtd"
        help = "Additional courant coundition for RT diffusion."
        arg_type = Float64
        default = 1.0
end

# read and parse command line arguments
arguments = parse_args(ARGS, s)

# include dispatch path and untility functions
MUST.@import_dispatch("../../../dispatch2/")
iniCond = MUST.ingredients("initial_condition.jl")

# Dispatch setup from command line
begin
    patch_size = arguments["patch_size"]                 
    τ_up = arguments["tau_up"]                    
    τ_surf = 0.0                    
    τ_down = arguments["tau_down"]                 
    τ_ee0 = arguments["tau_ee0"]                 
    τ_eemin = τ_up                  
    τ_zee0 = arguments["tau_zee0"]                  
    τ_rho0 =  arguments["tau_rho0"]                   
    dxdz_max = 3.0                  
    scale_resolution = arguments["scale_resolution"]        
    namelist_kwargs = Dict(         
        :newton_time=>arguments["newton_time"],        
        :friction_time=>arguments["friction_time"],     
        :newton_decay_scale=>arguments["newton_decay_scale"],  
        :friction_decay_scale=>arguments["friction_decay_scale"],  
        :courant_target=>0.25,
        :courant_hd=>arguments["courant_hd"],
        :courant_rt=>arguments["courant_rt"],
        :newton_params=>(            
            :on=>true,              
            :delay_rt=>true
        ),
        :io_params=>(
            :out_time=>arguments["out_time"],
            :end_time=>arguments["end_time"]
        ),
        :aux_params=>(
            :select=>arguments["aux_params"],
        ),
        :boundary_params=>(
            :upper_bc=>arguments["upper_bc"],
            :lower_bc=>arguments["lower_bc"],
            :radiation_bc=>arguments["radiation_bc"],
            :smallr=>arguments["smallr"]
        ),
        :an_params=>(
            :smallr=>arguments["smallr"],
            :dlnr_limit=>0.7
        ),
        :patch_params=>(
            :grace=>0.01,
            :nt=>5,
            :dt_small=>arguments["dt_small"]
        ),
        :sc_rt_params=>(
            :rt_grace=>0.01,
            :rt_freq=>arguments["rt_freq"],
            :cdtd=>arguments["cdtd"]
        ),
        :dispatcher0_params=>(
            :retry_stalled=>60,
        )
    )
end

# interpolate the average models
begin
    modelFolder = arguments["modelFolder"]
	avModels = joinpath(modelFolder, "av_models")
    !isdir(avModels) && mkpath(avModels)

    # check parameter input
    @assert length(arguments["teff"])==length(arguments["logg"])==length(arguments["feh"])
    nparas = length(arguments["teff"])
    paras = zeros(nparas, 3)
    paras[:, 1] .= arguments["teff"]
    paras[:, 2] .= arguments["logg"]
    paras[:, 3] .= arguments["feh"]

    # pick EoS
    eos = [e == "closest" ? :closest : e for e in arguments["eos"]]
    eos = if length(eos) == 1
        [first(eos) for _ in axes(paras, 1)]
    else
        @assert size(paras, 1) == length(eos)
        eos
    end

    ingrid = if arguments["grid"] == "Stagger"
        iniCond.staggergrid
    elseif arguments["grid"] == "MARCS"
        iniCond.marcsgrid
    else
        error("Given input grid doesn't exist.")
    end

    # interpolate in grid
    grid = iniCond.initialmodel(
        paras;
        grid=ingrid,
        eos=eos, 
        savedir=avModels, 
        τbottom=arguments["tau_bottom"], 
        common_size=3000,
        adiabatic_extrapolation=arguments["adiabatic_extrapolation"]
    )
end

# binning opacities
begin
    iniCond.prepare(
	    grid,
	    name_extension=modelFolder,
        version=arguments["version"],
        patch_size=patch_size,
        τ_up=τ_up,
        τ_surf=τ_surf,
        τ_down=τ_down,
        τ_ee0=τ_ee0,
        τ_eemin=τ_eemin,
        τ_zee0=τ_zee0,
        τ_rho0=τ_rho0,
        dxdz_max=dxdz_max,
        scale_resolution=scale_resolution,
        namelist_kwargs=namelist_kwargs
    )
end