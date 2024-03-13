TOOL  = ../OMtopogen/create_topog_refinedSampling.py \
	--source_file /home/Niki.Zadeh/safe/topogentestInput/GEBCO_2014_2D.nc \
        --source_lon lon \
        --source_lat lat \
        --source_elv elevation
OM4_025_preprocessing = /home/Niki.Zadeh/safe/topogentestInput
HGRID4deg = /home/Niki.Zadeh/safe/topogentestInput/ocean_hgrid.Merc.4deg.nc
TARGS = ocean_topog.Merc.4deg_r3b1h1.nc \
        ocean_topog.Merc.4deg_r3b2h1.nc \
        ocean_topog.Merc.4deg_r3b4h1.nc \
	ocean_topog.SCDP.OM4_r3b1h1.nc  \
	ocean_topog.BP.OM4_r3b4h1.nc    \
	ocean_topog.BP.OM4_r3b2h1.nc    \
	ocean_topog.SO.OM4_r3b2h1.nc    \
	ocean_topog.SO.OM4_r3b1h1.nc    \
	ocean_topog.SO.OM4_r5b1h1.nc    

all: $(TARGS)
	cat hash.md5
	md5sum -c hash.md5

ocean_topog.Merc.4deg_r3b1h1.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b2h1.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r3b4h1.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SCDP.OM4_r3b1h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/scap_supergrid.nc --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.BP.OM4_r3b4h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/ncap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.BP.OM4_r3b2h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/ncap_supergrid.nc --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SO.OM4_r3b2h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SO.OM4_r3b1h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

#default nxblocks and max_refine
ocean_topog.SO.OM4_r5b1h1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --outputfilename $@ --no_changing_meta

#hash.md5.gfdl-pan105: | $(TARGS)
#	md5sum $(TARGS) > $@
#	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) ocean_topog.*.nc
