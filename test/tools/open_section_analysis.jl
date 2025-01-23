using CUFSM, CrossSectionGeometry


#AISI D100-17E Cold-Formed Steel Design Manual Volume 1
#Example III-3


t = 0.0566
section_dimensions = [0.500, 1.625, 3.625, 1.625, 0.500]
r = [0.0849+t, 0.0849+t, 0.0849+t, 0.0849+t]
n = [3, 3, 3, 3, 3]
n_r = [3, 3, 3, 3];
θ = [π/2, π, -π/2, 0.0, π/2]

section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(section_dimensions, θ, n, r, n_r, t, centerline = "to left", offset = (section_dimensions[2], section_dimensions[3] - section_dimensions[1]))

E = 29500.0
ν = 0.30
P = 1.0
Mxx = 0.0
Mzz = 0.0
M11 = 0.0
M22 = 0.0
constraints = []
springs = []


lengths = 2.0:1.0:30
neigs = 1

x_center = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
y_center = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

model = CUFSM.Tools.open_section_analysis(x_center, y_center, t, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs, neigs)

# scale = (1.0, 1.0)
# eig = 1
# ax, figure = CUFSM.Show.signature_curve(model, eig, scale)
# figure

# drawing_scale = 1.0
# eig = 1
# deformation_scale = (1.0, 1.0)
# mode_index = 3
# ax, figure = CUFSM.Show.mode_shape(all_models[i], eig, mode_index, data.t[i]*ones(Float64, size(all_models[i].elem)[1]), deformation_scale, drawing_scale)
# figure



# CorZ = 0
# θ_top = 90.0
# #Calculate top flange + lip section properties.
# CorZ = 0

# b = 1.625 - t
# d = 0.500 - t/2
# θ = 90.0
# ho = 3.625
# μ = 0.30
# E = 29500.0
# G = E / (2 * (1 + μ))
# kϕ = 0.0
# Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,xhf,yhf,yof = AISIS100.v16S3.table_2_3_3__1(CorZ,t,b,d,θ)

# #Calculate the purlin distortional buckling half-wavelength.
# Lm = 999999999.0
# Lcrd = AISIS100.v16S3.appendix2_2_3_3_1__7(ho, μ, t, Ixf, xof, xhf, Cwf, Ixyf, Iyf)