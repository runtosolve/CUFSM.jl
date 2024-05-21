module Export


using MAT


function to_MAT(file_path, file_name, elem, lengths, node, prop)

    m_all_input = Matrix{Any}(undef, (1, size(lengths)[1]))

    for i in eachindex(lengths)
        m_all_input[i] = 1
    end

    matwrite(joinpath(file_path, file_name), Dict(
        "BC" => "S-S",
        "clas" => Matrix{Float64}(undef, 0, 0),
        "constraints" => 0,
        "curve" => Matrix{Float64}(undef, 0, 0),
        "elem" => elem,
        "GBTcon" => Dict{String, Any}("local"=>0.0, "ospace"=>1.0, "other"=>0.0, "couple"=>1.0, "orth"=>2.0, "dist"=>0.0, "glob"=>0.0),
        "lengths" => lengths,
        "m_all" =>  m_all_input,
        "node" => node,
        "prop" => prop,
        "shapes" =>  Matrix{Float64}(undef, 0, 0),
        "springs" => 0
    ); compress = false)

end

end  #module 