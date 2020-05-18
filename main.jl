using Format


include("input.jl")
include("read_pdb.jl")
include("dir_chk.jl")
include("psf_gen.jl")


cfg = read_config("config.in")
onemol = collect(read_pdb(i) for i in cfg.subName)

println(onemol)


#generate input flows
folder_check("flow")
generate_flow(onemol, cfg)

