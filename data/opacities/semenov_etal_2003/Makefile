OPAC_DIR = $(TOP)/data/opacities/$(OPACITY)

VPATH := $(VPATH):$(OPAC_DIR)

opacity_mod.o           : scaling_mod.o opacity

opacity:
	ln -s $(OPAC_DIR) opacity
