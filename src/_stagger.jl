#=========== Stagger module from Mikolaj ===========#
const stagger_endian = Ref("big-endian")
mutable struct StaggerMesh
    mx::Int64
    dxm::Vector{Float32}
    dxmdn::Vector{Float32}
    x::Vector{Float32}
    xmdn::Vector{Float32}
    dxidxup::Vector{Float32}
    dxidxdn::Vector{Float32}
    my::Int64
    dym::Vector{Float32}
    dymdn::Vector{Float32}
    y::Vector{Float32}
    ymdn::Vector{Float32}
    dyidyup::Vector{Float32}
    dyidydn::Vector{Float32}
    mz::Int64
    dzm::Vector{Float32}
    dzmdn::Vector{Float32}
    z::Vector{Float32}
    zmdn::Vector{Float32}
    dzidzup::Vector{Float32}
    dzidzdn::Vector{Float32}
    n::Int64

    function StaggerMesh(mesh_file::String)
        f = FortranFile(mesh_file, "r", convert=stagger_endian[]);
        dims = read(f,(Int32,3))
        # -- x direction
        mx = dims[1]
        my = dims[2]
        mz = dims[3]
        x = read(f,(Float32,(mx,6)))
        y = read(f,(Float32,(my,6)))
        z = read(f,(Float32,(mz,6)))
        new(mx,
            x[:,1],x[:,2],x[:,3],x[:,4],x[:,5],x[:,6],
            my,
            y[:,1],y[:,2],y[:,3],y[:,4],y[:,5],y[:,6],
            mz,
            z[:,1],z[:,2],z[:,3],z[:,4],z[:,5],z[:,6],
            mx*my*mz);
    end;
    function StaggerMesh(mesh::StaggerMesh, cutlb::Int64 = 0, cutub::Int64 = 0);
        new_my = size(mesh.y[1+cutlb:end-cutub])[1]
        new(mesh.mx,mesh.dxm,mesh.dxmdn,mesh.x,mesh.xmdn,mesh.dxidxup,mesh.dxidxdn,
            new_my,
            mesh.dym[1+cutlb:end-cutub],
            mesh.dymdn[1+cutlb:end-cutub],
            mesh.y[1+cutlb:end-cutub],
            mesh.ymdn[1+cutlb:end-cutub],
            mesh.dyidyup[1+cutlb:end-cutub],
            mesh.dyidydn[1+cutlb:end-cutub],
            mesh.mz,mesh.dzm,mesh.dzmdn,mesh.z,mesh.zmdn,mesh.dzidzup,mesh.dzidzdn,
            mesh.mx*new_my*mesh.mz);
    end;
end

mutable struct StaggerEOS

    mrho::Int32
    iupdte::Int32
    nvar::Int32
    mbox::Int32
    eosxmin::Float32
    dbox::Float32
    ul::Float32
    ut::Float32
    ur::Float32
    eps::Float32
    tff::Float32
    grv::Float32
    abnd::Float32

    # -- real arrays
    tmean::Vector{Float32}
    tamp::Vector{Float32}
    rhm::Vector{Float32} # rho axis
    xcorr::Vector{Float32}
    thmin::Vector{Float32}
    thmax::Vector{Float32}
    dth::Vector{Float32}
    eemin::Vector{Float32} # eemin/max
    eemax::Vector{Float32}
    deetab::Vector{Float32}

    # -- read int arrays
    itab::Vector{Int32}
    mtab::Vector{Int32}

    # --
    mtable::Int32
    tab::Vector{Float32}

    function StaggerEOS(eos_file::String)

        # open (12,file=tablefile,form='unformatted',status='old',readonly)
        # read (12) mrho,iupdte,nvar,mbox,eosxmin,dbox,ul,ut,ur,eps,tff,grv,abnd
        # allocate (tmean(mrho),tamp(mrho),rhm(mrho),xcorr(mrho),thmin(mrho),thmax(mrho))
        # allocate (dth(mrho),eemin(mrho),eemax(mrho),deetab(mrho),itab(mrho),mtab(mrho))
        # read (12) tmean,tamp,rhm,xcorr,thmin,thmax,dth,eemin,eemax,deetab,itab,mtab
        # read (12) mtable
        # allocate (tab(mtable))
        # read (12) (tab(i), i=1,mtable)
        # if (master) print *, 'init_eos: read   ', mtable, ' table values'
        # if (master) print *, 'init_eos: rho limits   =',rhm(1),rhm(mrho)
        # if (master) print *, 'init_eos: ee  limits   =',eemin(1),eemax(1)
        # if (master) print *, 'init_eos: T0_ideal, ee =',T0_ideal,ee0
        # close (12)

        f = FortranFile(eos_file);
        # mrho, iupdte, nvar, mbox = read(f,(Int32,4))
        # eosxmin,dbox,ul,ut,ur,eps,tff,grv,abnd = read(f,(Float32,9))
        t0, t1 = read(f,(Int32,4),(Float32,9))
        t2, t3 = read(f,(Float32,(t0[1],10)),(Int32,(t0[1],2)))

        mtable = read(f,Int32)
        tab    = read(f,(Float32,mtable))

        new(t0[1],t0[2],t0[3],t0[4],
            t1[1],t1[2],t1[3],t1[4],t1[5],t1[6],t1[7],t1[8],t1[9],
            t2[:,1],t2[:,2],t2[:,3],t2[:,4],t2[:,5],t2[:,6],t2[:,7],t2[:,8],t2[:,9],t2[:,10],
            t3[:,1],t3[:,2],
            mtable,
            tab)
    end;
end


mutable struct StaggerSnap
    mesh_file ::String
    dat_file  ::String
    aux_file  ::String
    mesh      ::StaggerMesh
    data      ::Dict{String,Array{Float32,3}}
    order     ::Vector{String}
end


function br_squeeze( A :: AbstractArray )
    keepdims = Tuple(i for i in size(A) if i != 1);
    return reshape( A, keepdims );
end

function br_heatmap_xz(A :: AbstractArray, M :: StaggerMesh; kwargs...)
    return heatmap(M.x,M.z,br_squeeze(A)',yflip=true; kwargs...)
end

function br_heatmap_yz(A :: AbstractArray, M :: StaggerMesh; kwargs...)
    return heatmap(M.y,M.z,br_squeeze(A)',yflip=true; kwargs...)
end

function br_heatmap_xy(A :: AbstractArray, M :: StaggerMesh; kwargs...)
    return heatmap(M.x,M.y,br_squeeze(A); kwargs...)
end

function br_arr_ffile(file_name::String, mesh::StaggerMesh; rpos::Int)
    f = FortranFile(file_name,"r",access="direct",recl=mesh.n*4, convert=stagger_endian[])
    var = read(f,rec=rpos,(Float32,(mesh.mx,mesh.my,mesh.mz)));
    close(f)
    return var
end

function StaggerMesh2BifrostMesh(snap::StaggerSnap; cutlb::Int64 = 0, cutub::Int64 = 0)

    M = snap.mesh
    my  = size(M.y[1+cutlb:end-cutub])[1]

    open(bifrost_name(snap.mesh_file),"w") do io
        println(io,@sprintf "%d" M.mz)
        println(io,join([@sprintf "%e" x for x in M.z], " "))
        println(io,join([@sprintf "%e" x for x in M.zmdn], " "))
        println(io,join([@sprintf "%e" x for x in M.dzidzup], " "))
        println(io,join([@sprintf "%e" x for x in M.dzidzdn], " "))
        println(io,@sprintf "%d" M.mx)
        println(io,join([@sprintf "%e" x for x in M.x], " "))
        println(io,join([@sprintf "%e" x for x in M.xmdn], " "))
        println(io,join([@sprintf "%e" x for x in M.dxidxup], " "))
        println(io,join([@sprintf "%e" x for x in M.dxidxdn], " "))
        # This dimension can be trimed
        println(io,@sprintf "%d"  my)
        println(io,join([@sprintf "%e" x for x in M.y[1+cutlb:end-cutub]], " "))
        println(io,join([@sprintf "%e" x for x in M.ymdn[1+cutlb:end-cutub]], " "))
        println(io,join([@sprintf "%e" x for x in M.dyidyup[1+cutlb:end-cutub]], " "))
        println(io,join([@sprintf "%e" x for x in M.dyidydn[1+cutlb:end-cutub]], " "))
    end
end

function StaggerSnap2BifrostSnap(file_name::String, mesh::StaggerMesh; do_mhd::Bool = false, cutlb::Int64 = 0, cutub::Int64 = 0)
    # create new mesh if cuts are needed
    if (cutlb > 0 || cutub > 0)
        fmesh = StaggerMesh(mesh,cutlb,cutub);
        StaggerMesh2BifrostMesh(mesh,cutlb,cutub);
    else
        fmesh = mesh;
    end
    # order of vars in Legacy Stagger: r,px,py,pz,e,d,Bx,By,Bz
    f = FortranFile("bifrost.snap","w",access="direct",recl=fmesh.n*4);
    tmp = br_arr_ffile(file_name, mesh, rpos=1) # r
    write(f, 1, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])   # permute: x,y,z --> z,x,y
    tmp = br_arr_ffile(file_name, mesh, rpos=4) # legacy pz -> px
    write(f, 2, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
    tmp = br_arr_ffile(file_name, mesh, rpos=2) # legacy px -> py
    write(f, 3, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
    tmp = br_arr_ffile(file_name, mesh, rpos=3) # legacy py -> pz
    write(f, 4, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
    tmp = br_arr_ffile(file_name, mesh, rpos=5) # e
    write(f, 5, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
    if do_mhd
        tmp = br_arr_ffile(file_name, mesh, rpos=9) # legacy bz -> bx
        write(f, 6, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
        tmp = br_arr_ffile(file_name, mesh, rpos=7) # legacy bz -> by
        write(f, 7, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
        tmp = br_arr_ffile(file_name, mesh, rpos=8) # legacy by -> bz
        write(f, 8, permutedims(tmp,[3,1,2])[:,:,1+cutlb:end-cutub])
    end
end

function StaggerSnap2BifrostSnap(snap::StaggerSnap; do_mhd::Bool = false, cutlb::Int64 = 0, cutub::Int64 = 0)
    # create new mesh if cuts are needed
    if (cutlb > 0 || cutub > 0)
        fmesh = StaggerMesh(snap.mesh,cutlb,cutub);
        #StaggerMesh2BifrostMesh(snap.mesh;cutlb,cutub);
    else
        fmesh = snap.mesh;
    end

    # order of vars in Legacy Stagger: r,px,py,pz,e,d,Bx,By,Bz
    f = FortranFile(bifrost_name(snap.dat_file), "w", access="direct", recl=fmesh.n*4);

    tmp = snap.data["r"]
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=1)   # permute: x,y,z --> z,x,y
    tmp = snap.data["pz"] # legacy pz -> px
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=2)
    tmp = snap.data["px"] # legacy px -> py
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=3)
    tmp = snap.data["py"] # legacy py -> pz
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=4)
    tmp = snap.data["e"] # e
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=5)

    if do_mhd
        tmp = snap.data["bz"] # legacy bz -> bx
        write(f, tmp[:,:,1+cutlb:end-cutub]; rec=6)
        tmp = snap.data["bx"] # legacy bx -> by
        write(f, tmp[:,:,1+cutlb:end-cutub]; rec=7)
        tmp = snap.data["by"] # legacy by -> bz
        write(f, tmp[:,:,1+cutlb:end-cutub]; rec=8)
    end

    close(f)
end

function StaggerAux2BifrostAux(snap::StaggerSnap; cutlb::Int64 = 0, cutub::Int64 = 0)
    # create new mesh if cuts are needed
    if (cutlb > 0 || cutub > 0)
        fmesh = StaggerMesh(snap.mesh,cutlb,cutub);
        #StaggerMesh2BifrostMesh(mesh,cutlb,cutub);
    else
        fmesh = snap.mesh;
    end

    # order of vars in Legacy Stagger: r,px,py,pz,e,d,Bx,By,Bz
    f = FortranFile(bifrost_name(snap.aux_file), "w", access="direct", recl=fmesh.n*4);

    tmp = exp.(snap.data["lpp"])
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=1)   

    tmp = exp.(snap.data["ltemp"])
    write(f, tmp[:,:,1+cutlb:end-cutub]; rec=2);

    close(f)

    nothing
end

function BifrostSnap(filename, folder="../input_data/")
    mesh_file = joinpath(folder,filename*".msh")
    dat_file  = joinpath(folder,filename*".dat")
    aux_file  = joinpath(folder,filename*".aux")
    mesh      = StaggerMesh(mesh_file)
    result    = Dict{String,Array{Float32,3}}()
    order    = ["r","px","py","pz","e","bx","by","bz"]
    for (i,para) in enumerate(order)
        try
            result[para] = br_arr_ffile(dat_file,mesh, rpos=i)
            #reverse!(result[para]; dims=3)
        catch e
            if isa(e, EOFError)
                @warn "$(para) not present."
            else
                error("Reading Failed.")
            end
        end
    end
    result["ee"] = result["e"] ./ result["r"] 

    aux_order = ["lpp", "ltemp"]
    for (i,para) in enumerate(aux_order)
        try
            result[para] = br_arr_ffile(aux_file, mesh, rpos=i)
            #reverse!(result[para]; dims=3)
        catch e
            if isa(e, EOFError)
                @warn "$(para) not present."
            else
                throw(e)
            end
        end
    end

    StaggerSnap(mesh_file, dat_file, aux_file, mesh, result, order)
end

"""
Read a Stagger Snapshot from .msh and .dat files.
"""
function StaggerSnap(filename, folder="../input_data/")
    mesh_file = joinpath(folder,filename*".msh")
    dat_file  = joinpath(folder,filename*".dat")
    aux_file  = joinpath(folder,filename*".aux")

    mesh_file = if !isfile(mesh_file)
        f = filename[1:findlast('_', filename)-1]
        joinpath(folder, f*".msh")
    else
        mesh_file
    end

    mesh      = StaggerMesh(mesh_file)
    result    = Dict{String,Array{Float32,3}}()
    order    = ["r","px","py","pz","e","d","bx","by","bz"]
    for (i,para) in enumerate(order)
        try
            result[para] = permutedims(br_arr_ffile(dat_file,mesh, rpos=i),[3,1,2])[:,:,2:end-1]
            #reverse!(result[para]; dims=3)
        catch e
            if isa(e, EOFError)
                #@warn "$(para) not present."
                nothing
            else
                error("Reading Failed.")
            end
        end
    end
    result["ee"] = result["e"] ./ result["r"] 

    aux_order = ["lpp", "lross", "ltemp", "lne", "lplanck", "ltau"]
    for (i,para) in enumerate(aux_order)
        try
            result[para] = permutedims(br_arr_ffile(aux_file, mesh, rpos=i),[3,1,2])[:,:,2:end-1]
            #reverse!(result[para]; dims=3)
        catch e
            if isa(e, EOFError)
                #@warn "$(para) not present."
                nothing
            else
                throw(e)
            end
        end
    end

    StaggerSnap(mesh_file, dat_file, aux_file, StaggerMesh(mesh, 1, 1), result, order)
end

"""
Apply the scaling in converters to snap (This is not ideal yet, should make it more flexible)
"""
function normalize!(snap, converter::AtmosUnits, paras...; reverse=false)
	operation(a,b) = reverse ? a ./ b : a .* b
	for p in paras 
		p == "r" ? 
		snap.data["r"] = operation(snap.data["r"],converter.d) : nothing
		p == "ee" ? 
		snap.data["ee"] = operation(snap.data["ee"],converter.ee) : nothing
		p == "e" ? 
		snap.data["e"] = operation(snap.data["e"],converter.e) : nothing
	end
end

"""
Lookup para in the given EOS using use_paras as variables.
"""
function get_from_eos(snap::StaggerSnap,eos,para="T",use_paras=["r","e"]; use_log=[true,true])
    snap_size = size(snap.data[snap.order[1]])
    result    = similar(snap.data[snap.order[1]])
    result   .= 0.0
    do_log    = [i ? x->log(x) : identity for i in use_log]
    @inbounds for k in 1:snap_size[3]
        @inbounds for j in 1:snap_size[2]
            @inbounds for i in 1:snap_size[1]
                try
                    result[i,j,k] = eos.lookup(para,
                                            do_log[1](snap.data[use_paras[1]][i,j,k]),
                                            do_log[2](snap.data[use_paras[2]][i,j,k]))
                catch 
                    result[i,j,k] = NaN32
                end
            end
        end
    end
    result
end

"""
Container for Python stagger snap (usefull for dispatching)
"""
struct StaggerLegacySnap
    snap::PyCall.PyObject
end

bifrost_name(stagger_name::String) = begin
    i_ext     = first(findlast(".", stagger_name))
    path, ext = stagger_name[1:i_ext-1], stagger_name[i_ext+1:end]

    path*"_bifrost.$(ext)"
end

function getindex(s::StaggerSnap, key::Symbol)
    return if key == :px
        s.data["pz"]
    elseif key == :py
        s.data["px"]
    elseif key == :pz
        s.data["py"] 
    elseif key == :x
        s.mesh.z
    elseif key == :y
        s.mesh.x
    elseif key == :z
        s.mesh.y
    elseif key == :T
        exp.(s.data["ltemp"])
    elseif key in Symbol.(["pp", "ross", "temp", "ne", "planck", "tau"])
        exp.(s.data["l$(key)"])
    else
        s.data[String(key)]
    end
end
