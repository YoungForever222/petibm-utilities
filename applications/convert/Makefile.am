bin_PROGRAMS = \
	petibm-convert2d \
	petibm-convert3d

petibm_convert2d_SOURCES = \
	convert.cpp
petibm_convert2d_CPPFLAGS = \
	-DDIMENSIONS=2 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_convert2d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)

petibm_convert3d_SOURCES = \
	convert.cpp
petibm_convert3d_CPPFLAGS = \
	-DDIMENSIONS=3 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_convert3d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
