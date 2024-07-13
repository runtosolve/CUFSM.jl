module Show

using LinearAlgebra, CairoMakie, LinesCurvesNodes


struct SignatureCurveOptions

    drawing_size::Vector{Float64}
    backgroundcolor::Symbol
    linecolor::Symbol
    linewidth::Float64
    fontsize::Float64
    markersize::Float64

end

struct ModeShapeOptions

    drawing_size::Vector{Float64}
    thickness_scale::Float64
    backgroundcolor::Symbol
    linecolor::Symbol
    linestyle::Symbol
    joinstyle::Symbol

end


using ..Tools 

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
        δ[i][1] = local_u_elem[1]
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


function find_maximum_element_deformed_dimensions(element_XY, Δ, scale)

    X = [element_XY[i][1] for i in eachindex(element_XY)]
    Y = [element_XY[i][2] for i in eachindex(element_XY)]

    ΔX = [Δ[i][1] for i in eachindex(element_XY)] * scale[1]
    ΔY = [Δ[i][2] for i in eachindex(element_XY)] * scale[2]

    max_X = maximum(X .+ ΔX)
    min_X = minimum(X .+ ΔX)
    max_Y = maximum(Y .+ ΔY)
    min_Y = minimum(Y .+ ΔY)

    # ΔX_deformed = max_X - min_X
    # ΔY_deformed = max_Y - min_Y

    return max_X, min_X, max_Y, min_Y

end


function cross_section_mode_shape_info(elem, node, mode, n, scale)

    #element lengths
    L = [norm([node[Int(elem[i, 3]), 2], node[Int(elem[i, 3]), 3]] - [node[Int(elem[i, 2]), 2], node[Int(elem[i, 2]), 3]]) for i in eachindex(elem[:,1])]

    #element rotation matrices
    Γ = Vector{Matrix{Float64}}(undef, length(L))
    for i in eachindex(L)
        Γ[i] = LinesCurvesNodes.define_rotation_matrix([node[Int(elem[i, 2]), 2], node[Int(elem[i, 2]), 3]], [node[Int(elem[i, 3]), 2], node[Int(elem[i, 3]), 3]]) 
    end
    # Γ_2D = [Γ[i][[1, 2, 6, 7, 8, 12],[1, 2, 6, 7, 8, 12]] for i in eachindex(L)]  #reduce to 2D rotation matrix, u1, v1, θ1, u1, v2, θ2

    global_dof_elem = define_element_global_dof(elem, node)

    global_u_elem = [mode[global_dof_elem[i][[1, 3, 4, 5, 7, 8]]]  for i in eachindex(global_dof_elem)]  #u1 w1 θ1 u1 w2 θ2 

    local_u_elem = [Γ[i] * global_u_elem[i] for i in eachindex(global_u_elem)]

    δ, Δ, x = define_cross_section_deformation(local_u_elem, L, Γ, n)

    cross_section_coords = define_cross_section_geometry(elem, node, Γ, x)

    max_X = Vector{Float64}(undef, size(Δ, 1))
    min_X = Vector{Float64}(undef, size(Δ, 1))
    max_Y = Vector{Float64}(undef, size(Δ, 1))
    min_Y = Vector{Float64}(undef, size(Δ, 1))
    # ΔY_deformed = Vector{Float64}(undef, size(Δ, 1))
    for i in eachindex(Δ)
        max_X[i], min_X[i], max_Y[i], min_Y[i] = find_maximum_element_deformed_dimensions(cross_section_coords[i], Δ[i], scale)
    end

    X_range_deformed = maximum(max_X) - minimum(min_X)
    Y_range_deformed = maximum(max_Y) - minimum(min_Y)
    
    figure_max_dim_range = (X_range_deformed, Y_range_deformed)

    return cross_section_coords, Δ, figure_max_dim_range

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


function show_element_deformed_shape!(ax, element_XY, Δ, scale, attributes, line_thickness)

    X = [element_XY[i][1] for i in eachindex(element_XY)]
    Y = [element_XY[i][2] for i in eachindex(element_XY)]

    ΔX = [Δ[i][1] for i in eachindex(element_XY)]
    ΔY = [Δ[i][2] for i in eachindex(element_XY)]


    lines!(ax, X .+ scale[1] * ΔX, Y .+ scale[2] * ΔY, linestyle=attributes.linestyle, color=attributes.color, linewidth=line_thickness, joinstyle=attributes.joinstyle, linecap=attributes.linecap)

    # for i=1:(length(X)-1)

    #     scatterlines!(ax, [X[i] + scale[1] * ΔX[i], X[i+1] + scale[1] * ΔX[i+1]], [Y[i] + scale[2] * ΔY[i], Y[i+1] + scale[2] * ΔY[i+1]], linestyle=attributes.linestyle, color=attributes.color, linewidth=line_thickness, marker = attributes.marker, markersize=attributes.markersize)

    # end

end

function signature_curve(model, eig, options)

    Pcr = Tools.get_load_factor(model, eig)
    figure = Figure(size = (options.drawing_size[1], options.drawing_size[2]), backgroundcolor=options.backgroundcolor, fontsize = options.fontsize)
    ax = Axis(figure[1, 1], backgroundcolor=options.backgroundcolor)
    [scatterlines!(model.lengths, Pcr, color=options.linecolor, markersize=options.markersize, linewidth = options.linewidth) for i in eachindex(Pcr)];
    
    return ax, figure 

end

function signature_curve(model, eig, options, figure)

    Pcr = Tools.get_load_factor(model, eig)
    # figure = Figure(size = (options.drawing_size[1], options.drawing_size[2]), backgroundcolor=options.backgroundcolor, fontsize = options.fontsize)
    ax = Axis(figure[1, 1], backgroundcolor=options.backgroundcolor)
    [scatterlines!(model.lengths, Pcr, color=options.linecolor, markersize=options.markersize, linewidth = options.linewidth) for i in eachindex(Pcr)];
    
    return ax, figure 

end

#function signature_curve(model, eig, scale; figure = Figure())
# function signature_curve(model, eig, scale, figure)

#     Pcr = Tools.get_load_factor(model, eig)
#     # figure = Figure(resolution = (4.0*72, 4.0*72))
#     ax = Axis(figure[1, 1])
#     [scatterlines!(model.lengths .* scale[1], Pcr .* scale[2], color=:blue) for i in eachindex(Pcr)];
    
#     return ax, figure 

# end


function minimum_mode_shape(model, eig, t, deformation_scale, options)
    

    Pcr = Tools.get_load_factor(model, eig)
    mode_index = argmin(Pcr)
    mode = model.shapes[mode_index][:, eig]

    n = fill(10, length(t))

    cross_section_coords, Δ, figure_max_dims = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)

    XY_discretized =Float64[]
    t_discretized = Float64[]
    Δ_discretized = Float64[]

    for i in eachindex(cross_section_coords)

        if i != size(cross_section_coords)[1]
            XY_segment = [cross_section_coords[i][j] for j in eachindex(cross_section_coords[i])][1:end-1]
        else
            XY_segment = [cross_section_coords[i][j] for j in eachindex(cross_section_coords[i])] #last segment, add end node 
        end

        XY_discretized = vcat(XY_discretized, XY_segment)

        t_discretized = vcat(t_discretized, fill(t[i], n[i]-1))

        if i != size(cross_section_coords)[1]
            Δ_segment = [Δ[i][j] for j in eachindex(Δ[i])][1:end-1]
        else
            Δ_segment = [Δ[i][j] for j in eachindex(Δ[i])]  #last segment, add end node 
        end

        Δ_discretized = vcat(Δ_discretized, Δ_segment)

    end

    X_discretized = [XY_discretized[i][1] for i in eachindex(XY_discretized)]
    Y_discretized = [XY_discretized[i][2] for i in eachindex(XY_discretized)]

    ΔX_discretized = [Δ_discretized[i][1] for i in eachindex(Δ_discretized)]
    ΔY_discretized = [Δ_discretized[i][2] for i in eachindex(Δ_discretized)]

    linesegment_ranges, t_segments = LinesCurvesNodes.find_linesegments(t_discretized)
    coords_as_linesegments = LinesCurvesNodes.combine_points_into_linesegments(linesegment_ranges, X_discretized, Y_discretized)
    Δ_as_linesegments = LinesCurvesNodes.combine_points_into_linesegments(linesegment_ranges, ΔX_discretized, ΔY_discretized)


    # figure_span_x = figure_max_dims[1]
    # figure_span_y = figure_max_dims[2]

    # figure = Figure(size = (figure_span_x*72, figure_span_y*72) .* drawing_scale)
    figure = Figure(size = (options.drawing_size[1], options.drawing_size[2]), backgroundcolor=options.backgroundcolor)
    # ax = Axis(figure[1, 1], aspect = figure_span_x/figure_span_y)

    ax = Axis(figure[1, 1], backgroundcolor = options.backgroundcolor)


    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
    # thickness_scale = maximum(t) * 72 * drawing_scale
    # linewidths = t ./ maximum(t) * thickness_scale

    # thickness_scale = maximum(t_segments) * 72 * drawing_scale
    linewidths = t_segments .* options.thickness_scale

    attributes = (color=options.linecolor, linestyle=options.linestyle, joinstyle=options.joinstyle, linecap=:butt)

    # linesegment_ranges = LinesCurvesNodes.find_line_segments(linewidths)

    # [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], deformation_scale, attributes, attributes.linewidth[i]) for i in eachindex(Δ)];

    [show_element_deformed_shape!(ax, coords_as_linesegments[i], Δ_as_linesegments[i], deformation_scale, attributes, linewidths[i]) for i in eachindex(Δ_as_linesegments)];

    num_nodes = size(model.node)[1]
    num_elem = size(model.elem)[1]

    # if length(X_discretized) == length(t_discretized)  #closed section
    if num_nodes == num_elem

        coords_last_linesegment = Vector{Vector{Float64}}(undef, 3)
        Δ_last_linesegment = Vector{Vector{Float64}}(undef, 3)

        coords_last_linesegment[1] = coords_as_linesegments[end][end]
        coords_last_linesegment[2] = coords_as_linesegments[1][1]
        coords_last_linesegment[3] = coords_as_linesegments[1][2]

        Δ_last_linesegment[1] = Δ_as_linesegments[end][end]
        Δ_last_linesegment[2] = Δ_as_linesegments[1][1]
        Δ_last_linesegment[3] = Δ_as_linesegments[1][2]

        attributes = (color=options.linecolor, linestyle=options.linestyle, joinstyle=options.joinstyle, linecap=:round) #just for connecting element


        show_element_deformed_shape!(ax, coords_last_linesegment, Δ_last_linesegment, deformation_scale, attributes, linewidths[end])

    end
    
    return ax, figure

end



function minimum_mode_shape(model, eig, t, deformation_scale, options, figure)
    

    Pcr = Tools.get_load_factor(model, eig)
    mode_index = argmin(Pcr)
    mode = model.shapes[mode_index][:, eig]

    n = fill(10, length(t))

    cross_section_coords, Δ, figure_max_dims = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)

    XY_discretized =Float64[]
    t_discretized = Float64[]
    Δ_discretized = Float64[]

    for i in eachindex(cross_section_coords)

        if i != size(cross_section_coords)[1]
            XY_segment = [cross_section_coords[i][j] for j in eachindex(cross_section_coords[i])][1:end-1]
        else
            XY_segment = [cross_section_coords[i][j] for j in eachindex(cross_section_coords[i])] #last segment, add end node 
        end

        XY_discretized = vcat(XY_discretized, XY_segment)

        t_discretized = vcat(t_discretized, fill(t[i], n[i]-1))

        if i != size(cross_section_coords)[1]
            Δ_segment = [Δ[i][j] for j in eachindex(Δ[i])][1:end-1]
        else
            Δ_segment = [Δ[i][j] for j in eachindex(Δ[i])]  #last segment, add end node 
        end

        Δ_discretized = vcat(Δ_discretized, Δ_segment)

    end

    X_discretized = [XY_discretized[i][1] for i in eachindex(XY_discretized)]
    Y_discretized = [XY_discretized[i][2] for i in eachindex(XY_discretized)]

    ΔX_discretized = [Δ_discretized[i][1] for i in eachindex(Δ_discretized)]
    ΔY_discretized = [Δ_discretized[i][2] for i in eachindex(Δ_discretized)]

    linesegment_ranges, t_segments = LinesCurvesNodes.find_linesegments(t_discretized)
    coords_as_linesegments = LinesCurvesNodes.combine_points_into_linesegments(linesegment_ranges, X_discretized, Y_discretized)
    Δ_as_linesegments = LinesCurvesNodes.combine_points_into_linesegments(linesegment_ranges, ΔX_discretized, ΔY_discretized)


    # figure_span_x = figure_max_dims[1]
    # figure_span_y = figure_max_dims[2]

    # figure = Figure(size = (figure_span_x*72, figure_span_y*72) .* drawing_scale)
    # figure = Figure(size = (options.drawing_size[1], options.drawing_size[2]), backgroundcolor=options.backgroundcolor)
    # ax = Axis(figure[1, 1], aspect = figure_span_x/figure_span_y)

    ax = Axis(figure[1, 1], backgroundcolor = options.backgroundcolor)


    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
    # thickness_scale = maximum(t) * 72 * drawing_scale
    # linewidths = t ./ maximum(t) * thickness_scale

    # thickness_scale = maximum(t_segments) * 72 * drawing_scale
    linewidths = t_segments .* options.thickness_scale

    attributes = (color=options.linecolor, linestyle=options.linestyle, joinstyle=options.joinstyle)

    # linesegment_ranges = LinesCurvesNodes.find_line_segments(linewidths)

    # [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], deformation_scale, attributes, attributes.linewidth[i]) for i in eachindex(Δ)];

    [show_element_deformed_shape!(ax, coords_as_linesegments[i], Δ_as_linesegments[i], deformation_scale, attributes, linewidths[i]) for i in eachindex(Δ_as_linesegments)];

    num_nodes = size(model.node)[1]
    num_elem = size(model.elem)[1]

    # if length(X_discretized) == length(t_discretized)  #closed section
    if num_nodes == num_elem

        coords_last_linesegment = Vector{Vector{Float64}}(undef, 3)
        Δ_last_linesegment = Vector{Vector{Float64}}(undef, 3)

        coords_last_linesegment[1] = coords_as_linesegments[end][end]
        coords_last_linesegment[2] = coords_as_linesegments[1][1]
        coords_last_linesegment[3] = coords_as_linesegments[1][2]

        Δ_last_linesegment[1] = Δ_as_linesegments[end][end]
        Δ_last_linesegment[2] = Δ_as_linesegments[1][1]
        Δ_last_linesegment[3] = Δ_as_linesegments[1][2]

        show_element_deformed_shape!(ax, coords_last_linesegment, Δ_last_linesegment, deformation_scale, attributes, linewidths[end])

    end
    
    return ax, figure

end




# function minimum_mode_shape(model, eig, t, deformation_scale, drawing_scale, options, figure)
    
#     # x = model.node[:, 2]
#     # y = model.node[:, 3]
#     # Δx = abs(maximum(x) - minimum(x))
#     # Δy = abs(maximum(y) - minimum(y))

#     Pcr = Tools.get_load_factor(model, eig)
#     mode_index = argmin(Pcr)
#     mode = model.shapes[mode_index][:, eig]

#     n = fill(5, length(t))

#     cross_section_coords, Δ, figure_max_dims = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)
   
#     Δx = figure_max_dims[1]
#     Δy = figure_max_dims[2]

#     # figure = Figure(resolution = (Δx*72, Δy*72) .* drawing_scale)
#     ax = Axis(figure[1, 1], aspect = Δx/Δy)
#     hidedecorations!(ax)  # hides ticks, grid and lables
#     hidespines!(ax)  # hide the frame
#     thickness_scale = maximum(t) * 72 * drawing_scale
#     linewidths = t ./ maximum(t) * thickness_scale
    
#     attributes = (color=:grey, linestyle=:solid, linewidth=linewidths, marker=:circle, markersize=0)
#     [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], deformation_scale, attributes, attributes.linewidth[i]) for i in eachindex(Δ)];
   

#     # cross_section_mode_shape!(ax, model.elem, model.node, mode, n, deformation_scale, attributes);
    
#     return ax, figure

# end







# function minimum_mode_shape(model, eig, t, deformation_scale, drawing_scale, figure = Figure())
    
#     # x = model.node[:, 2]
#     # y = model.node[:, 3]
#     # Δx = abs(maximum(x) - minimum(x))
#     # Δy = abs(maximum(y) - minimum(y))

#     Pcr = Tools.get_load_factor(model, eig)
#     mode_index = argmin(Pcr) 
#     mode = model.shapes[mode_index][:, eig]

#     n = fill(5, length(t))

#     cross_section_coords, Δ, figure_max_dims = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)
   
#     Δx = figure_max_dims[1]
#     Δy = figure_max_dims[2]

#     # figure = Figure(resolution = (Δx*72, Δy*72) .* drawing_scale)
#     ax = Axis(figure[1, 1], aspect = Δx/Δy)
#     hidedecorations!(ax)  # hides ticks, grid and lables
#     hidespines!(ax)  # hide the frame
#     thickness_scale = maximum(t) * 72 * drawing_scale
#     linewidths = t ./ maximum(t) * thickness_scale
    
#     attributes = (color=:grey, linestyle=:solid, linewidth=linewidths, marker=:circle, markersize=0)
#     [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], deformation_scale, attributes, attributes.linewidth[i]) for i in eachindex(Δ)];
   
#     figure
#     # cross_section_mode_shape!(ax, model.elem, model.node, mode, n, deformation_scale, attributes);
    
#     # return ax, figure

# end



function mode_shape(model, eig, mode_index, t, deformation_scale, drawing_scale)
    
    mode = model.shapes[mode_index][:, eig]

    n = fill(5, length(t))

    cross_section_coords, Δ, figure_max_dims = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)
   
    Δx = figure_max_dims[1]
    Δy = figure_max_dims[2]

    figure = Figure(size = (Δx*72, Δy*72) .* drawing_scale)
    ax = Axis(figure[1, 1], aspect = Δx/Δy)
    hidedecorations!(ax)  # hides ticks, grid and lables
    hidespines!(ax)  # hide the frame
    thickness_scale = maximum(t) * 72 * drawing_scale
    linewidths = t ./ maximum(t) * thickness_scale
    
    attributes = (color=:grey, linestyle=:solid, linewidth=linewidths, marker=:circle, markersize=0)
    [show_element_deformed_shape!(ax, cross_section_coords[i], Δ[i], deformation_scale, attributes, attributes.linewidth[i]) for i in eachindex(Δ)];
   

    # cross_section_mode_shape!(ax, model.elem, model.node, mode, n, deformation_scale, attributes);
    
    return ax, figure

end


function get_deformed_shape_extents(model, eig, deformation_scale)

        
    Pcr = Tools.get_load_factor(model, eig)
    mode_index = argmin(Pcr)
    mode = model.shapes[mode_index][:, eig]
    t = model.elem[:, 4]
    n = fill(5, length(t))
    cross_section_coords, Δ_nodes, figure_max_dim_range = cross_section_mode_shape_info(model.elem, model.node, mode, n, deformation_scale)

    return figure_max_dim_range

end

function define_drawing_size(Δ, max_pixels)

    max_index = argmax(Δ)
    min_index = argmin(Δ)

    if max_index == min_index
        min_index = 2
    end

    drawing_size = Vector{Float64}(undef, 2)
    
    drawing_size[max_index] = max_pixels
    drawing_size[min_index] = max_pixels * Δ[min_index] / Δ[max_index]

    thickness_scale = (max_pixels / Δ[max_index])

    return drawing_size, thickness_scale

end


end #module 