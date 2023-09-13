TOOL  = ../OMtopogen/create_topog_refinedSampling.py
HGRID = /home/Niki.Zadeh/safe/topogentestInput/ocean_hgrid.Merc.4deg.nc
TARGS = ocean_topog.Merc.4deg_r3b1h1.nc \
        ocean_topog.Merc.4deg_r3b2h1.nc \
        ocean_topog.Merc.4deg_r3b4h1.nc

all: $(TARGS)
	cat hash.md5
	md5sum -c hash.md5

ocean_topog.Merc.4deg_r3b1h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b2h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b4h1.nc:
	time $(TOOL) --hgridfilename $(HGRID) --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta

#hash.md5.gfdl-pan105: | $(TARGS)
#	md5sum $(TARGS) > $@
#	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) ocean_topog.*.nc
