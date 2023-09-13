TOOL  = ../OMtopogen/create_topog_refinedSampling.py
HGRID = ./dataset/ocean_hgrid.Merc.4deg.nc
OCEAN_TOPOG_SRC  = --source_file ./dataset/GEBCO_2014_2D.nc --source_lon lon --source_lat lat --source_elv elevation

TARGS = ocean_topog.Merc.4deg_r3b1h1.nc \
        ocean_topog.Merc.4deg_r3b2h1.nc \
        ocean_topog.Merc.4deg_r3b4h1.nc \
        ocean_topog.Merc.4deg_r3b4c4noargsh1.nc

all: $(TARGS)
	cat hash.md5
	md5sum -c hash.md5

ocean_topog.Merc.4deg_r3b1h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=1 --max_refine=3 $(OCEAN_TOPOG_SRC) --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b2h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=2 --max_refine=3 $(OCEAN_TOPOG_SRC) --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b4h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=4 --max_refine=3 $(OCEAN_TOPOG_SRC) --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b4c4noargsh1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=4 --ncores=4 --max_refine=3 $(OCEAN_TOPOG_SRC) --outputfilename $@ --no_changing_meta

#hash.md5.gfdl-pan105: | $(TARGS)
#	md5sum $(TARGS) > $@
#	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) ocean_topog.*.nc
