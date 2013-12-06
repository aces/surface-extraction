bin_PROGRAMS = surface_fit check_self_intersect create_anchor_constraints skel \
               adapt_object_mesh inflate_to_sphere interpolate_sphere \
               sphere_mesh

surface_fit_SOURCES = 	evaluate_fit_3d.c \
		 	boundary_search.c \
			conjugate_grad.c \
			poly_poly_dist.c \
			surf_surf_dist.c \
			self_intersect.c \
			clip_search.c \
			line_movement.c \
			surface_fit.c \
			gw_gradient.c

adapt_object_mesh_SOURCES = adapt_object_mesh.c

inflate_to_sphere_SOURCES = inflate_to_sphere.c

interpolate_sphere_SOURCES = interpolate_sphere.c \
                             octree.c

sphere_mesh_SOURCES = sphere.cc

check_self_intersect_SOURCES = 	check_self_intersect.c \
				self_intersect.c \
				poly_poly_dist.c

create_anchor_constraints_SOURCES = create_anchor_constraints.c

skel_SOURCES = skeletonization.c

noinst_HEADERS = conjugate.h\
	deform.h \
	deform_prototypes.h \
	fit_3d.h \
	fit_3d_prototypes.h \
	interval.h \
	octree.h \
	special_geometry.h \
	sp_geom_prototypes.h \
	subdiv.h \
	subdiv_prototypes.h \
	tri_mesh.h \
	skeletonization.h
