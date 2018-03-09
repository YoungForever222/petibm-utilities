bin_PROGRAMS = \
	petibm-vorticity2d \
	petibm-vorticity3d

petibm_vorticity2d_SOURCES = \
	vorticity2d.cpp
petibm_vorticity2d_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_vorticity2d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)

petibm_vorticity3d_SOURCES = \
	vorticity3d.cpp
petibm_vorticity3d_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_vorticity3d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
