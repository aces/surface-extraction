AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src

bin_SCRIPTS = scripts/asp \
   scripts/calibrate_white \
   scripts/clasp \
   scripts/classify_correct \
   scripts/cortical_thickness \
   scripts/expand_from_white \
   scripts/extract_white_surface \
   scripts/make_asp_grid \
   scripts/mask_cortical_white_matter \
   scripts/surface_qc


pkgdata_DATA = models/white_matter_mask.obj \
   models/white_model_320.obj \
   models/ellipsoid_1024.obj.gz \
   models/ellipsoid_1280.obj.gz \
   models/ellipsoid_128.obj.gz \
   models/ellipsoid_131072.obj.gz \
   models/ellipsoid_16384.obj.gz \
   models/ellipsoid_20480.obj.gz \
   models/ellipsoid_2048.obj.gz \
   models/ellipsoid_256.obj.gz \
   models/ellipsoid_320.obj.gz \
   models/ellipsoid_32768.obj.gz \
   models/ellipsoid_4096.obj.gz \
   models/ellipsoid_5120.obj.gz \
   models/ellipsoid_512.obj.gz \
   models/ellipsoid_81920.obj.gz \
   models/ellipsoid_327680.obj.gz \
   models/ellipsoid_8192.obj.gz


EXTRA_DIST = \
   $(pkgdata_DATA) \
   $(bin_SCRIPTS) \
   epm-header.in \
   m4/smr_WITH_BUILD_PATH.m4 \
   m4/mni_REQUIRE_LIB.m4 \
   m4/mni_REQUIRE_MNILIBS.m4


