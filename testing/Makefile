#ln -s /home/Niki.Zadeh/datasets .datasets
TOOL  = ../OMtopogen/create_topog_refinedSampling.py \
	--source_file .datasets/GEBCO_2020.nc \
        --source_lon lon \
        --source_lat lat \
        --source_elv elevation
OM4_025_preprocessing = .datasets/OM4p25
HGRID4deg = .datasets/ocean_hgrid.Merc.4deg.nc
TARGS = ocean_topog.Merc.4deg_r5xb1.nc \
        ocean_topog.Merc.4deg_r5xb2.nc \
        ocean_topog.Merc.4deg_r5xb3.nc \
        ocean_topog.Merc.4deg_r5xb30.nc \
	ocean_topog.SCDP.OM4_r3xb1.nc  \
	ocean_topog.BP.OM4_r3xb4.nc    \
	ocean_topog.BP.OM4_r3xb2.nc    \
	ocean_topog.SO.OM4_r3xb2.nc    \
	ocean_topog.SO.OM4_r3xb1.nc    \
	ocean_topog.SO.OM4_r3xb2yb2.nc

all: $(TARGS)
	cat hash.md5
	md5sum -c hash.md5

ocean_topog.Merc.4deg_r5xb1.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=1 --max_refine=5 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r5xb2.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=2 --max_refine=5 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r5xb3.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=3 --max_refine=5 --outputfilename $@ --no_changing_meta

ocean_topog.Merc.4deg_r5xb30.nc:
	time $(TOOL) --hgridfilename $(HGRID4deg) --nxblocks=30 --max_refine=5 --outputfilename $@ --no_changing_meta

ocean_topog.SCDP.OM4_r3xb1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/scap_supergrid.nc --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.BP.OM4_r3xb4.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/ncap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.BP.OM4_r3xb2.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/ncap_supergrid.nc --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SO.OM4_r3xb2.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --nxblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SO.OM4_r3xb1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --nxblocks=1 --max_refine=3 --outputfilename $@ --no_changing_meta

ocean_topog.SO.OM4_r3xb2yb2.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --nxblocks=2 --nyblocks=2 --max_refine=3 --outputfilename $@ --no_changing_meta

#default nxblocks and max_refine
ocean_topog.SO.OM4_r5xb1.nc:
	time $(TOOL) --hgridfilename $(OM4_025_preprocessing)/antarctic_spherical_supergrid.nc --outputfilename $@ --no_changing_meta

hash.md5.gfdl-lscamd50-d: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS) ocean_topog.*.nc
