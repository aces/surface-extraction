perllibdir = $(pkglibdir)

# slightly confusing ... the sources for the two scripts below have
# a .plib extension, which will get changed (after parsing) to .pl

unchecked_libs = deform_utils.plib utils.plib
checked_libs = $(unchecked_libs:.plib=.pl)
$(checked_libs): %.pl: %.plib
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

perllib_DATA = deform_utils.pl utils.pl

bin_SCRIPTS = asp expand_from_white extract_white_surface \
			  mask_cortical_white_matter make_asp_grid \
              clasp1 clasp2 calibrate_white classify_correct

EXTRA_DIST =  asp.pl expand_from_white.pl extract_white_surface.pl \
			  mask_cortical_white_matter.pl deform_utils.plib utils.plib \
              make_asp_grid.pl clasp1.pl clasp2.pl calibrate_white.pl \
              classify_correct.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = asp expand_from_white extract_white_surface \
			  mask_cortical_white_matter deform_utils.pl utils.pl \
			  make_asp_grid clasp1 clasp2 calibrate_white
