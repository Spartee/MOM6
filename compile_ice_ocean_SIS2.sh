mkdir -p build/gnu/shared/repro/
(cd build/gnu/shared/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ../../../../src/FMS; \
../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/ncrc-gnu.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names)
(cd build/gnu/shared/repro/; source ../../env; make NETCDF=3 REPRO=1 libfms.a -j)

mkdir -p build/gnu/ice_ocean_SIS2/repro/
(cd build/gnu/ice_ocean_SIS2/repro/; rm -f path_names; \
../../../../src/mkmf/bin/list_paths -l ./ ../../../../src/MOM6/config_src/{dynamic,coupled_driver,external} ../../../../src/MOM6/src/{*,*/*}/ ../../../../src/{atmos_null,coupler,land_null,ice_param,icebergs,SIS2,FMS/coupler,FMS/include}/)
(cd build/gnu/ice_ocean_SIS2/repro/; \
../../../../src/mkmf/bin/mkmf -t ../../../../src/mkmf/templates/ncrc-gnu.mk -o '-I../../shared/repro' -p MOM6 -l '-L../../shared/repro -lfms' -c '-Duse_AM3_physics -D_USE_LEGACY_LAND_' path_names )
(cd build/gnu/ice_ocean_SIS2/repro/; source ../../env; make NETCDF=3 REPRO=1 MOM6 -j)

