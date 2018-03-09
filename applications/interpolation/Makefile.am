bin_PROGRAMS = \
	petibm-interpolation2d \
	petibm-interpolation3d \
	petibm-interpolation2d-test \
	petibm-interpolation3d-test

petibm_interpolation2d_SOURCES = \
	interpolation.cpp
petibm_interpolation2d_CPPFLAGS = \
	-DDIMENSIONS=2 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_interpolation2d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)

petibm_interpolation3d_SOURCES = \
	interpolation.cpp
petibm_interpolation3d_CPPFLAGS = \
	-DDIMENSIONS=3 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_interpolation3d_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)

petibm_interpolation2d_test_SOURCES = \
	test.cpp
petibm_interpolation2d_test_CPPFLAGS = \
	-DDIMENSIONS=2 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_interpolation2d_test_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)

petibm_interpolation3d_test_SOURCES = \
	test.cpp
petibm_interpolation3d_test_CPPFLAGS = \
	-DDIMENSIONS=3 \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_interpolation3d_test_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
