###################################################
### Scratch space to store cdo terminal scripts ###
###################################################

# 2024-06-26

################################################
### Processing historical CMIP6 sea ice data ###
################################################

# use cdo scripts to convert from native grid to regular 1x1 degree grid
# outputted netcdf can then be read by R

### ACCESS-CM2 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/SImon/siconc/gn/20200817
cdo remapbil,global_1 siconc_SImon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### ACCESS-ESM1-5 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/SImon/siconc/gn/20200817
cdo remapbil,global_1 siconc_SImon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### BCC-CSM2-MR ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/BCC/BCC-CSM2-MR/historical/r1i1p1f1/SImon/siconc/gn/20200218
cdo remapbil,global_1 siconc_SImon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_BCC-CSM2-MR_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### CanESM5 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p1f1/SImon/siconc/gn/20190429
cdo remapbil,global_1 siconc_SImon_CanESM5_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_CanESM5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### CESM2-WACCM ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/SImon/siconc/gn/20190227
cdo remapbil,global_1 siconc_SImon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### CMCC-CM2-SR5 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CMCC/CMCC-CM2-SR5/historical/r1i1p1f1/SImon/siconc/gn/20200616
cdo remapbil,global_1 siconc_SImon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### EC-Earth3 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/historical/r1i1p1f1/SImon/siconc/gn/20200918
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

### EC-Earth3-Veg ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3-Veg/historical/r1i1p1f1/SImon/siconc/gn/20211207
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

### FGOALS-g3 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CAS/FGOALS-g3/historical/r1i1p1f1/SImon/siconc/gn/20210108
cdo remapbil,global_1 siconc_SImon_FGOALS-g3_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_FGOALS-g3_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### IPSL-CM6A-LR ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/SImon/siconc/gn/20180803
cdo remapbil,global_1 siconc_SImon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### MIROC6 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/MIROC/MIROC6/historical/r1i1p1f1/SImon/siconc/gn/20181212
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

### MPI-ESM1-2-LR ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/SImon/siconc/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

### MRI-ESM2-0 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/SImon/siconc/gn/20190904
cdo remapbil,global_1 siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_MRI-ESM2-0_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### NESM3 ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/NUIST/NESM3/historical/r1i1p1f1/SImon/siconc/gn/20190704
cdo remapbil,global_1 siconc_SImon_NESM3_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_NESM3_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### NorESM2-LM ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/NCC/NorESM2-LM/historical/r1i1p1f1/SImon/siconc/gn/20190815
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done



#################
#################
#################


##################
### ACCESS-CM2 ###
##################

### historical ###
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/Omon/mlotst/gn/20191108
cdo remapbil,global_1 mlotst_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc mlotst_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/Omon/sos/gn/20191108
cdo remapbil,global_1 sos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc sos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/Omon/tos/gn/20191108
cdo remapbil,global_1 tos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc tos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/Omon/zos/gn/20191108
cdo remapbil,global_1 zos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc zos_Omon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/CMIP/CSIRO-ARCCSS/ACCESS-CM2/historical/r1i1p1f1/SImon/siconc/gn/20200817
cdo remapbil,global_1 siconc_SImon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc siconc_SImon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412_bil_1x1.nc

### SSP245 ###
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245/r1i1p1f1/Omon/mlotst/gn/20191108
cdo remapbil,global_1 mlotst_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012.nc mlotst_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245/r1i1p1f1/Omon/sos/gn/20191108
cdo remapbil,global_1 sos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012.nc sos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245/r1i1p1f1/Omon/tos/gn/20191108
cdo remapbil,global_1 tos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012.nc tos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245/r1i1p1f1/Omon/zos/gn/20191108
cdo remapbil,global_1 zos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012.nc zos_Omon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp245/r1i1p1f1/SImon/siconc/gn/20200817
cdo remapbil,global_1 siconc_SImon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012.nc siconc_SImon_ACCESS-CM2_ssp245_r1i1p1f1_gn_201501-210012_bil_1x1.nc

### ssp585 ###
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585/r1i1p1f1/Omon/mlotst/gn/20210317
cdo remapbil,global_1 mlotst_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012.nc mlotst_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585/r1i1p1f1/Omon/sos/gn/20210317
cdo remapbil,global_1 sos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012.nc sos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585/r1i1p1f1/Omon/tos/gn/20210317
cdo remapbil,global_1 tos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012.nc tos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585/r1i1p1f1/Omon/zos/gn/20210317
cdo remapbil,global_1 zos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012.nc zos_Omon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012_bil_1x1.nc
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/CSIRO-ARCCSS/ACCESS-CM2/ssp585/r1i1p1f1/SImon/siconc/gn/20210317
cdo remapbil,global_1 siconc_SImon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012.nc siconc_SImon_ACCESS-CM2_ssp585_r1i1p1f1_gn_201501-210012_bil_1x1.nc

####################
### IPSL-CM6A-LR ###
####################

# there are major issues with the IPSL data sets - why?

# process netcdf files
cd /Users/home/ANTSIE\ data/CMIP6/ScenarioMIP/IPSL/IPSL-CM6A-LR/ssp245/r1i1p1f1/Omon/tos/gn/20190119

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapbil,r360x180  "mlotst_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012.nc" "mlotst_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012_bil.nc"
done

cdo remapbil,r360x180 "mlotst_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012.nc"


cdo remapnn,global_1 "tos_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012.nc" "test.nc"



cdo remapbil,global_1 "mlotst_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_185001-201412.nc" "test.nc"


cdo remapbil,global_1 tos_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012.nc test2.nc

cdo ncatted -O -a units,lat,c,c,"degrees north" -a units,lon,c,c,"degrees east" tos_Omon_IPSL-CM6A-LR_ssp245_r1i1p1f1_gn_201501-210012.nc


                
                
                