lib_LTLIBRARIES = libpetibm-utilities.la

libpetibm_utilities_la_SOURCES = \
	field.cpp \
	grid.cpp \
	misc.cpp \
	timestep.cpp \
	vorticity.cpp

libpetibm_utilities_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)

libpetibm_utilities_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
