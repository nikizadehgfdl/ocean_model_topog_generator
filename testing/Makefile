#ln -s /home/Niki.Zadeh/datasets .datasets
TOOL_NUMPY  = ../OMtopogen/create_topog_refinedSampling.py
TOOL_TORCH  = ../OMtopogen/create_topog_refinedSampling_torch.py
SOURCE_ARGS = --source_file .datasets/GEBCO_2023.nc --source_lon lon --source_lat lat --source_elv elevation
HGRID4deg = .datasets/ocean_hgrid.Merc.4deg.nc
HGRID_OM4p25_subgrids = .datasets/OM4p25
HGRID_OM4p25 = .datasets/OM4p25.v20170622/ocean_hgrid.nc
HGRID_OM5p25 = .datasets/OM5p25.v20240410/ocean_hgrid.nc
HGRID_OM4p0833 = .datasets/ocean_hgrid_res0.0833.nc
HGRID_OM4p0625 = .datasets/ocean_hgrid_res0.0625.nc
HGRID_OM4p125 = .datasets/ocean_hgrid_OM4p125v20210310.nc

TARGS = $(FAST_NUMPY) $(FAST_TOCPU) $(FAST_TOGPU) #$(QUICK_NUMPY) $(QUICK_TOCPU)  $(QUICK_TOGPU)

FAST_NUMPY = ocean_topog.Merc.4deg_r4xb10_NUMPY.nc 
FAST_TOCPU = ocean_topog.Merc.4deg_r4xb10_torchCPU.nc 
FAST_TOGPU = ocean_topog.Merc.4deg_r4xb10_torchGPU.nc 

QUICK_NUMPY = ocean_topog.BP.OM4_r3xb4_NUMPY.nc ocean_topog.MERC.OM4_r3xb4_NUMPY.nc \
              ocean_topog.SCAP.OM4_r3xb4_NUMPY.nc ocean_topog.SO.OM4_r3xb4_NUMPY.nc ocean_topog.BP.OM4_r3xb24yb4_NUMPY.nc
QUICK_TOCPU = ocean_topog.BP.OM4_r3xb4_torchCPU.nc ocean_topog.MERC.OM4_r3xb4_torchCPU.nc \
              ocean_topog.SCAP.OM4_r3xb4_torchCPU.nc ocean_topog.SO.OM4_r3xb4_torchCPU.nc ocean_topog.BP.OM4_r3xb24yb4_torchCPU.nc
QUICK_TOGPU = ocean_topog.BP.OM4_r3xb4_torchGPU.nc ocean_topog.MERC.OM4_r3xb4_torchGPU.nc \
              ocean_topog.SCAP.OM4_r3xb4_torchGPU.nc ocean_topog.SO.OM4_r3xb4_torchGPU.nc ocean_topog.BP.OM4_r3xb24yb4_torchGPU.nc

SLOW_NUMPY  = ocean_topog.SO.OM4_r7xb4_NUMPY.nc
SLOW_TOCPU  = ocean_topog.SO.OM4_r7xb4_torchCPU.nc
SLOW_TOGOU  = ocean_topog.SO.OM4_r7xb4_torchGPU.nc

all: $(TARGS)
	cat hash.md5
	md5sum -c hash.md5

#QUICK TESTS
ocean_topog.Merc.4deg_r4xb10_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID4deg) --nxblocks=10 --max_refine=4 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.Merc.4deg_r4xb10_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID4deg) --nxblocks=10 --max_refine=4 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
ocean_topog.Merc.4deg_r4xb10_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID4deg) --nxblocks=10 --max_refine=4 --outputfilename $@ --no_changing_meta > $@.stdout

ocean_topog.BP.OM4_r3xb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r3xb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r3xb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
ocean_topog.MERC.OM4_r3xb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.MERC.OM4_r3xb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.MERC.OM4_r3xb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
ocean_topog.SCAP.OM4_r3xb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/scap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SCAP.OM4_r3xb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/scap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SCAP.OM4_r3xb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/scap_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
ocean_topog.SO.OM4_r3xb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SO.OM4_r3xb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SO.OM4_r3xb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#test xb and yb
ocean_topog.BP.OM4_r3xb24yb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r3xb24yb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r3xb24yb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=3 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#test max_refine
ocean_topog.BP.OM4_r5xb24yb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=5 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r5xb24yb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=5 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r5xb24yb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=5 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#
#SLOW!
ocean_topog.BP.OM4_r6xb24yb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=6 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r6xb24yb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=6 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r6xb24yb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=6 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#
#SLOWER!
ocean_topog.SO.OM4_r7xb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SO.OM4_r7xb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.SO.OM4_r7xb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/antarctic_spherical_supergrid.nc --nxblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#Very SLOW
ocean_topog.MERC.OM4_r7xb24yb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.MERC.OM4_r7xb24yb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.MERC.OM4_r7xb24yb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/mercator_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout
#
#Increasing max_refine to 7 the answers diverge between NUMPY and TORCH for the Bipolar Northern Cap grid
#
ocean_topog.BP.OM4_r7xb24yb4_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r7xb24yb4_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.BP.OM4_r7xb24yb4_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p25_subgrids)/ncap_supergrid.nc --nxblocks=24 --nyblocks=4 --max_refine=7 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout

#
#PRODUCTION TESTS
#
#VERY SLOW!
#OM5p25 global
ocean_topog.OM5p25v20240410_r8xb24yb387_NUMPY.nc:
	time $(TOOL_NUMPY) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM5p25) --nxblocks=24 --nyblocks=387 --max_refine=8 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.OM5p25v20240410_r8xb24yb387_torchCPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM5p25) --nxblocks=24 --nyblocks=387 --max_refine=8 --outputfilename $@ --no_changing_meta > $@.stdout
ocean_topog.OM5p25v20240410_r8xb24yb387_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM5p25) --nxblocks=24 --nyblocks=387 --max_refine=8 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout

ocean_topog.OM4p125_r6xb20yb20_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p125) --nxblocks=20 --nyblocks=20 --max_refine=6 --outputfilename $@ --no_changing_meta --device cuda > $@.stdout

ocean_topog_res0.0833_r8xb480yb18_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p0833) --nxblocks=480 --nyblocks=18 --max_refine=8 --outputfilename $@  --device cuda > $@.stdout

ocean_topog_res0.0625_r8xb640yb16_torchGPU.nc:
	time $(TOOL_TORCH) $(SOURCE_ARGS) --hgridfilename $(HGRID_OM4p0625) --nxblocks=640 --nyblocks=16 --max_refine=8 --outputfilename $@  --device cuda > $@.stdout

#If the hash file is not present for you platform, create it
#then make a symlink to it called hash.md5
#The checksums are not the same for all platforms or versions of the code
hash.md5.fast.gfdl-lscamd50-d: | $(TARGS)
	md5sum $(TARGS) > $@
	cat $@

check:
	md5sum -c hash.md5

clean:
	rm -f $(TARGS)
