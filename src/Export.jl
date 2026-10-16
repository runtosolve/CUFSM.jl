module Export


# to_MAT is implemented in ext/CUFSMMATExt.jl, which Julia loads automatically
# once MAT is loaded alongside CUFSM (`using CUFSM, MAT`). Keeping MAT a weak
# dependency means projects that never write .mat files do not install
# MAT -> HDF5 -> MPI.
function to_MAT end

to_MAT(args...) = error("CUFSM.Export.to_MAT needs MAT.jl: add it to your environment and run `using MAT`.")

end  #module
