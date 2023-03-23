module Show

using LinearAlgebra, CairoMakie, LinesCurvesNodes



function define_element_global_dof(elem, node)

    global_dof_elem = Vector{Vector{Int}}(undef, size(elem, 1))

    num_nodes = size(node, 1)
    for i in eachindex(elem[:,1])

        node_i_index = Int(elem[i, 2])
        node_j_index = Int(elem[i, 3])

        global_dof_node_i = [node_i_index*2 - 1, node_i_index*2, 2*num_nodes + node_i_index*2 - 1, 2*num_nodes + node_i_index*2] 
        global_dof_node_j = [node_j_index*2 - 1, node_j_index*2, 2*num_nodes + node_j_index*2 - 1, 2*num_nodes + node_j_index*2] 
        global_dof_elem[i] = [global_dof_node_i; global_dof_node_j]  #u, v, w, θ

    end

    return global_dof_elem

end

function get_element_deformed_shape(local_u_elem, L, Γ, n)

    x = range(0.0, L, n)

    dof = [2, 3, 5, 6]
    q1, q2, q3, q4 = local_u_elem[dof]
    w_xy = LinesCurvesNodes.beam_shape_function.(q1, q2, q3, q4, L, x)

    #local element deformation
    δ = [zeros(Float64, 2) for i in eachindex(x)]

    for i in eachindex(x)
        δ[i][2] = w_xy[i]
    end

    #global element deformation
    Δ = [Γ'[1:2,1:2] * δ[i] for i in eachindex(x)]

    return δ, Δ, x

end



function define_cross_section_deformation(local_u_elem, L, Γ, n)

    δ = Vector{Vector{Vector{Float64}}}(undef, length(L))
    Δ = Vector{Vector{Vector{Float64}}}(undef, length(L))
    x = Vector{Vector{Float64}}(undef, length(L))

    for i in eachindex(L)

        δ[i], Δ[i], x[i] = get_element_deformed_shape(local_u_elem[i], L[i], Γ[i], n[i])

    end

    return δ, Δ, x

end


function cross_section_mode_shape!(ax, elem, node, mode, n, scale, attributes)

    #element lengths
    L = [norm([node[Int(elem[i, 3]), 2], node[Int(elem[i, 3]), 3]] - [node[Int(elem[i, 2]), 2], node[Int(elem[i, 2]), 3]]) for i in eachindex(elem[:,1])]

    #element rotation matrices
    Γ = Vector{Matrix{Float64}}(undef, length(L))
    for i in eachindex(L)
        Γ[i], angles = LinesCurvesNodes.define_rotation_matrix([node[Int(elem[i, 2]), 2], node[Int(elem[i, 2]), 3], 0.0], [node[Int(elem[i, 3]), 2], node[Int(elem[i, 3]), 3], 0.0], 0.0) 
    end
    # Γ_2D = [Γ[i][[1, 2, 6, 7, 8, 12],[1, 2, 6, 7, 8, 12]] for i in eachindex(L)]  #reduce to 2D rotation matrix, u1, v1, θ1, u1, v2, θ2

    global_dof_elem = define_element_global_dof(elem, node)

    global_u_elem = [mode[global_dof_elem[i][[1, 3, 4, 5, 7, 8]]] for i in eachindex(global_dof_elem)]  #u1 w1 θ1 u1 w2 θ2 

    local_u_elem = [Γ[i][[1, 2, 6, 7, 8, 12],[1, 2, 6, 7, 8, 12]] * global_u_elem[i] for i in eachindex(global_u_elem)]

    δ, Δ, x = define_cross_section_deformation(local_u_elem, L, Γ, n)

    cross_section_coords = define_cross_section_geometry(elem, node, Γ, x)

    [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], scale, attributes) for i in eachindex(Δ)];
   
end;


function define_cross_section_geometry(elem, node, Γ, x)

    cross_section_coords = Vector{Vector{Vector{Float64}}}(undef, length(x))
    for i in eachindex(x)
        node_i = [node[Int(elem[i, 2]), 2], node[Int(elem[i, 2]), 3], 0.0]
        cross_section_coords_3D = LinesCurvesNodes.discretized_element_global_coords(node_i, Γ[i], x[i])
        cross_section_coords[i] = [cross_section_coords_3D[i][1:2] for i in eachindex(cross_section_coords_3D)]

    end

    return cross_section_coords

end


function show_element_deformed_shape!(ax, element_XY, Δ, scale, attributes)

    X = [element_XY[i][1] for i in eachindex(element_XY)]
    Y = [element_XY[i][2] for i in eachindex(element_XY)]

    ΔX = [Δ[i][1] for i in eachindex(element_XY)]
    ΔY = [Δ[i][2] for i in eachindex(element_XY)]

    for i=1:(length(X)-1)

        scatterlines!(ax, [X[i] + scale[1] * ΔX[i], X[i+1] + scale[1] * ΔX[i+1]], [Y[i] + scale[2] * ΔY[i], Y[i+1] + scale[2] * ΔY[i+1]], linestyle=attributes.linestyle, color=attributes.color, linewidth=attributes.linewidth, marker = attributes.marker, markersize=attributes.markersize)

    end

end

end #module 