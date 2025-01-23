module Tools

using ..CUFSM 

# function view_closed_section_mode_shape(node, shapes, mode_index, scale_x, scale_y)

#     num_nodes = size(node)[1]

#     mode = shapes[mode_index]

#     mode_x = mode[1:2:2*num_nodes]
#     mode_y = mode[(2*num_nodes + 1):2:4*num_nodes]

#     defx = node[:, 2] .+ scale_x * mode_x
#     defy = node[:, 3] .+scale_y * mode_y

#     #For a closed cross-section, add first node to end for plotting.
#     defx = [defx; defx[1]]
#     defy = [defy; defy[1]]

#     #Define undeformed shape.
#     undefx = [node[:,2]; node[2,1]]
#     undefy = [node[:,3]; node[3,2]]

#     plot(undefx, undefy, size = (600,600), legend = false)
#     plot!(defx, defy, markershape = :o)

# end

function closed_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)

    #Define coord matrix.
    coord = [x_center y_center]

    #Define number of cross-section elements.
    num_elem = length(x_center)

    #Define number of nodes.
    num_nodes = length(x_center)

    start_nodes = 1:num_nodes
    end_nodes = [2:num_nodes; 1]
    ends = [start_nodes end_nodes t]

    # #Define element connectivity
    # ends = zeros(num_elem, 3)
    # ends[:, 1] .= 1:(num_nodes-1)
    # ends[:, 2] .= 2:num_nodes
    # ends[:, 3] .= t

    #Calculate section properties.
    section_properties = CUFSM.cutwp_prop2(coord,ends)

    #Map section properties to CUFSM.
    A = section_properties.A
    xcg = section_properties.xc
    zcg = section_properties.yc
    Ixx = section_properties.Ixx
    Izz = section_properties.Iyy
    Ixz = section_properties.Ixy
    thetap = section_properties.θ
    I11 = section_properties.I1
    I22 = section_properties.I2
    unsymm = 0;  #Sets Ixz=0 if unsymm = 0

    #Define the number of cross-section nodes.
    num_cross_section_nodes = size(coord)[1]

    #Initialize CUFSM node matrix.
    node = zeros(Float64, (num_cross_section_nodes, 8))

    #Add node numbers to node matrix.
    node[:, 1] .= 1:num_cross_section_nodes

    #Add nodal coordinates to node matrix.
    node[:, 2:3] .= coord

    #Add nodal restraints to node matrix.
    node[:, 4:7] .= ones(num_cross_section_nodes,4)

    #Initialize CUFSM elem matrix.
    elem = zeros(Float64, (num_elem, 5))

    #Add element numbers to elem matrix.
    elem[:, 1] = 1:num_elem

    #Add element connectivity and thickness to elem matrix.
    elem[:, 2:4] .= ends

    #Add element material reference to elem matrix.
    elem[:, 5] .= ones(num_elem) * 100
                            
    #There are no springs or constraints.
    # springs = []
    # constraints = 0

    #Define material properties.
    G = E / (2 *(1 + ν))
    prop = [100 E E ν ν G]

    neigs = 1  #just need the first mode 

    #Add reference stresses to node matrix.
    node = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    #Run CUFSM.
    curve, shapes = CUFSM.strip(prop, node, elem, lengths, springs, constraints, neigs)

    model = CUFSM.Model(prop=prop, node=node, elem=elem, lengths=lengths, springs=springs, constraints=constraints, neigs=neigs, curve=curve, shapes=shapes)


    return model

end


# function view_multi_branch_section_mode_shape(node, elem, shapes, mode_index, scale_x, scale_y, xlims, ylims)

#     num_nodes = size(node)[1]

#     mode = shapes[mode_index]

#     undefx = node[:, 2]
#     undefy = node[:, 3]

#     mode_x = mode[1:2:2*num_nodes]
#     mode_y = mode[(2*num_nodes + 1):2:4*num_nodes]

#     defx = node[:, 2] .+ scale_x * mode_x
#     defy = node[:, 3] .+scale_y * mode_y

#     #For a multi-branch cross-section, plot each element individually.

#     num_elem = size(elem)[1]

#     for i = 1:num_elem

#         node_i = Int(elem[i, 2])
#         node_j = Int(elem[i, 3])

#         if i == 1

#             plot([undefx[node_i], undefx[node_j]], [undefy[node_i], undefy[node_j]], size = (600, 600), legend = false, linecolor = :black, xlims = xlims, ylims = ylims, aspect_ratio = :equal)
#             plot!([defx[node_i], defx[node_j]], [defy[node_i], defy[node_j]], legend = false, linecolor = :blue)
    
#         else
#             plot!([undefx[node_i], undefx[node_j]], [undefy[node_i], undefy[node_j]], legend = false, linecolor = :black)
#             plot!([defx[node_i], defx[node_j]], [defy[node_i], defy[node_j]], legend = false, linecolor = :blue)

#         end

#     end

#     return current()

# end



function open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    ####Calculate section properties.

    #Define coord matrix.
    coord = [x_center y_center]

    #Define number of cross-section elements.
    num_elem = length(x_center) - 1

    #Define number of nodes.
    num_nodes = length(x_center)

    #Define element connectivity
    ends = zeros(num_elem, 3)
    ends[:, 1] .= 1:(num_nodes-1)
    ends[:, 2] .= 2:num_nodes
    ends[:, 3] .= t

    section_properties = CUFSM.cutwp_prop2(coord,ends)

    #Map section properties to CUFSM.
    A = section_properties.A
    xcg = section_properties.xc
    zcg = section_properties.yc
    Ixx = section_properties.Ixx
    Izz = section_properties.Iyy
    Ixz = section_properties.Ixy
    thetap = section_properties.θ
    I11 = section_properties.I1
    I22 = section_properties.I2
    unsymm = 0  #Sets Ixz=0 if unsymm = 0

    #Define the number of cross-section nodes.
    num_cross_section_nodes = size(coord)[1]

    #Initialize CUFSM node matrix.
    node = zeros(Float64, (num_cross_section_nodes, 8))

    #Add node numbers to node matrix.
    node[:, 1] .= 1:num_cross_section_nodes

    #Add nodal coordinates to node matrix.
    node[:, 2:3] .= coord

    #Add nodal restraints to node matrix.
    node[:, 4:7] .= ones(num_cross_section_nodes,4)

    #Initialize CUFSM elem matrix.
    elem = zeros(Float64, (num_elem, 5))

    #Add element numbers to elem matrix.
    elem[:, 1] = 1:num_elem

    #Add element connectivity and thickness to elem matrix.
    elem[:, 2:4] .= ends

    #Add element material reference to elem matrix.
    elem[:, 5] .= ones(num_elem) * 100
                            
    #There are no springs or constraints.
    # springs = []
    # constraints = 0

    #Define material properties.
    G = E / (2 *(1 + ν))
    prop = [100 E E ν ν G]

    # neigs = 1  #just need the first mode 

    #Add reference stresses to node matrix.
    node = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    #Run CUFSM.
    curve, shapes = CUFSM.strip(prop, node, elem, lengths, springs, constraints, neigs)

    #collect up 
    model = CUFSM.Model(prop, node, elem, lengths, springs, constraints, neigs, curve, shapes)

    return model

end




function open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, supports, neigs)

    ####Calculate section properties.

    #Define coord matrix.
    coord = [x_center y_center]

    #Define number of cross-section elements.
    num_elem = length(x_center) - 1

    #Define number of nodes.
    num_nodes = length(x_center)

    #Define element connectivity
    ends = zeros(num_elem, 3)
    ends[:, 1] .= 1:(num_nodes-1)
    ends[:, 2] .= 2:num_nodes
    ends[:, 3] .= t

    section_properties = CUFSM.cutwp_prop2(coord,ends)

    #Map section properties to CUFSM.
    A = section_properties.A
    xcg = section_properties.xc
    zcg = section_properties.yc
    Ixx = section_properties.Ixx
    Izz = section_properties.Iyy
    Ixz = section_properties.Ixy
    thetap = section_properties.θ
    I11 = section_properties.I1
    I22 = section_properties.I2
    unsymm = 0  #Sets Ixz=0 if unsymm = 0

    #Define the number of cross-section nodes.
    num_cross_section_nodes = size(coord)[1]

    #Initialize CUFSM node matrix.
    node = zeros(Float64, (num_cross_section_nodes, 8))

    #Add node numbers to node matrix.
    node[:, 1] .= 1:num_cross_section_nodes

    #Add nodal coordinates to node matrix.
    node[:, 2:3] .= coord

    #Add nodal restraints to node matrix.
    node[:, 4:7] .= ones(num_cross_section_nodes,4)

    #Initialize CUFSM elem matrix.
    elem = zeros(Float64, (num_elem, 5))

    #Add element numbers to elem matrix.
    elem[:, 1] = 1:num_elem

    #Add element connectivity and thickness to elem matrix.
    elem[:, 2:4] .= ends

    #Add element material reference to elem matrix.
    elem[:, 5] .= ones(num_elem) * 100
                            
    #There are no springs or constraints.
    # springs = []
    # constraints = 0

    #Define material properties.
    G = E / (2 *(1 + ν))
    prop = [100 E E ν ν G]

    # neigs = 1  #just need the first mode 

    #Add reference stresses to node matrix.
    node = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    #Add supports:
    for i in eachindex(supports)

        node[supports[i][1], 4:7] = supports[i][2:5]

    end


    #Run CUFSM.
    curve, shapes = CUFSM.strip(prop, node, elem, lengths, springs, constraints, neigs)

    #collect up 
    model = CUFSM.Model(prop=prop, node=node, elem=elem, lengths=lengths, springs=springs, constraints=constraints, neigs=neigs, curve=curve, shapes=shapes)

    return model

end




function open_section_analysis(x_center, y_center, ends::Matrix{Float64}, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

    ####Calculate section properties.

    #Define coord matrix.
    coord = [x_center y_center]

    #Define number of cross-section elements.
    num_elem = length(x_center) - 1

    #Define number of nodes.
    num_nodes = length(x_center)

    # #Define element connectivity
    # ends = zeros(num_elem, 3)
    # ends[:, 1] .= 1:(num_nodes-1)
    # ends[:, 2] .= 2:num_nodes
    # ends[:, 3] .= t

    section_properties = CUFSM.cutwp_prop2(coord,ends)

    #Map section properties to CUFSM.
    A = section_properties.A
    xcg = section_properties.xc
    zcg = section_properties.yc
    Ixx = section_properties.Ixx
    Izz = section_properties.Iyy
    Ixz = section_properties.Ixy
    thetap = section_properties.θ
    I11 = section_properties.I1
    I22 = section_properties.I2
    unsymm = 0  #Sets Ixz=0 if unsymm = 0

    #Define the number of cross-section nodes.
    num_cross_section_nodes = size(coord)[1]

    #Initialize CUFSM node matrix.
    node = zeros(Float64, (num_cross_section_nodes, 8))

    #Add node numbers to node matrix.
    node[:, 1] .= 1:num_cross_section_nodes

    #Add nodal coordinates to node matrix.
    node[:, 2:3] .= coord

    #Add nodal restraints to node matrix.
    node[:, 4:7] .= ones(num_cross_section_nodes,4)

    #Initialize CUFSM elem matrix.
    elem = zeros(Float64, (num_elem, 5))

    #Add element numbers to elem matrix.
    elem[:, 1] = 1:num_elem

    #Add element connectivity and thickness to elem matrix.
    elem[:, 2:4] .= ends

    #Add element material reference to elem matrix.
    elem[:, 5] .= ones(num_elem) * 100
                            
    #There are no springs or constraints.
    # springs = []
    # constraints = 0

    #Define material properties.
    G = E / (2 *(1 + ν))
    prop = [100 E E ν ν G]

    # neigs = 1  #just need the first mode 

    #Add reference stresses to node matrix.
    node = CUFSM.stresgen(node,P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    #Run CUFSM.
    curve, shapes = CUFSM.strip(prop, node, elem, lengths, springs, constraints, neigs)

    #collect up 
    model = CUFSM.Model(prop=prop, node=node, elem=elem, lengths=lengths, springs=springs, constraints=constraints, neigs=neigs, curve=curve, shapes=shapes)

    return model

end







function get_load_factor(model, eig)

    load_factor = [model.curve[i][eig, 2] for i in eachindex(model.curve)]

end


end #module 
