#-----------------------------------------------------------------------------------------#
###         Identifying hidden Zika hotspots in Pernambuco, Brazil:                     ###
###                               A spatial analysis                                    ###
#-----------------------------------------------------------------------------------------#

## Laís Picinini Freitas 

## Last update: 8 March 2021

## Preprint: 

## Objective: 
# To identify hotspots of Zika in Pernambuco state, Northeast Brazil, using Aedes-borne diseases (dengue, chikungunya and Zika) 
# and microcephaly data. Kulldorff’s Poisson purely spatial scan statistic was used to detect low- and high-risk clusters 
# and then we combined the results.

rm(list=ls())

# Packages ----------------------------------------------------------------

library(sf) # to read the shapefile
library(tidyverse) # to organize the data
library(colorspace) # palettes for the figures
library(rsatscan) # SaTScan (must be installed in the machine. To download SaTScan: www.satscan.org)


# Data --------------------------------------------------------------------

## Reading Aedes-borne diseases confirmed cases and microcephaly cases by Municipality and year, 2013-2017
# Dengue, Zika and Chikungunya data originally obtained from SINAN at <ftp://ftp.datasus.gov.br/dissemin/publicos/SINAN/DADOS>. 
# Live births and microcephaly data originally obtained from SINASC at <ftp://ftp.datasus.gov.br/dissemin/publicos/SINASC/>.
# Population projections by municipality estimated by Freire et al available at <https://demografiaufrn.net/laboratorios/lepp/>.
arboPE <- read_csv('data_by_year.csv')

## Reading shapefile for PE state (without Fernando de Noronha island) municipalities  
# Original source: IGBE at <https://www.ibge.gov.br/geociencias/organizacao-do-territorio/malhas-territoriais>.
shape <- read_sf('PE_semilha.shp')  %>% 
  mutate(ID_MUNICIP = str_sub(CD_GEOCMU,1,6)) %>% st_set_crs(4326)


## Organizing the data for the cluster and risk analysis
dzcPE <- arboPE %>% 
  ungroup() %>% 
  select(ID_MUNICIP,year,Population,dengue,zika,chikungunya) %>% 
  # Selecting data 2014-2017
  filter(year>2013) %>% 
  # Aggregating all years
  group_by(ID_MUNICIP) %>% 
  summarise(dengue=sum(dengue),
            zika=sum(zika),
            chikungunya=sum(chikungunya),
            population=mean(Population)) %>% 
  # Aggregating the three diseases and caculating the expected cases (E)
  mutate(dzc=dengue+zika+chikungunya,
         E=(sum(dzc)/sum(population))*population,
         # For the satscan analysis we need this variable:
         time='unspecified')

microPE <- arboPE %>% 
  ungroup() %>% 
  select(ID_MUNICIP,year,microcephaly,live_births) %>% 
  # Selecting data 2015-2017
  filter(year>2014) %>% 
  # Aggregating all years
  group_by(ID_MUNICIP) %>% 
  summarise(microcephaly=sum(microcephaly),
            live_births=sum(live_births)) %>% 
  # Caculating the expected cases (E)
  mutate(E=(sum(microcephaly)/sum(live_births))*live_births,
         # For the satscan analysis we need this variable:
         time='unspecified')


# Scan statistics (satscan) -----------------------------------------------

## Creating the centroids for each municipality

centroids <- st_coordinates(st_centroid(shape)[,1:2])
centroids <- as.data.frame(cbind(shape,centroids)) %>% select(ID_MUNICIP,X,Y,-geometry)


# Loading SaTScan and parameters for poisson spatiotemporal analysis

ss.local <- "/home/Lais/SaTScan/"  # this need to be changed to the path SaTScan is installed
ssenv$.ss.params.OLD <-  ssenv$.ss.params
ssenv$.ss.params <- readLines('template_satscan.txt')


# Dengue, Zika and Chikungunya --------------------------------------------

ARQ <- 'dataPE.cas'

write.table(dzcPE[,c(1,6)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(dzcPE[,c(1,8,5)],file='PE.pop' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(centroids[,c(1,3,2)],file='PE.geo' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

ss.options(list(CaseFile=ARQ,
                StartDate="2014/01/01",
                EndDate="2017/12/31", 
                PrecisionCaseTimes=0, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic     
                PopulationFile="PE.pop",
                CoordinatesFile="PE.geo", 
                CoordinatesType=1, # 0=Cartesian, 1=latitude/longitude
                AnalysisType=1, # 1=Purely Spatial, 3=Retrospective Space-Time
                ModelType=0, # 0=Discrete Poisson
                ScanAreas=3, # 1=High, 2=Low, 3=High and Low
                TimeAggregationUnits=1, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
                TimeAggregationLength=1,
                MaxSpatialSizeInPopulationAtRisk=50,
                MinimumTemporalClusterSize=1,
                MaxTemporalSize=1,
                MinimumCasesInHighRateClusters=5,
                MaxSpatialSizeInPopulationAtRisk_Reported=50,
                UseDistanceFromCenterOption_Reported='y',
                MaxSpatialSizeInDistanceFromCenter_Reported=50 # max size radius cluster in km                
))


ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))

modelo <- paste0(ARQ,'_modelo')
write.ss.prm(tmp,modelo)
result_dzcPE <-  satscan(tmp,modelo, sslocation="/home/lais/SaTScan/")

summary(result_dzcPE)
result_dzcPE


### Results

rgis_dzc <- result_dzcPE$gis %>% 
  filter(P_VALUE < 0.05) %>% 
  mutate(type = case_when(
    CLU_RR > 1 ~ 'high',
    CLU_RR < 1 ~ 'low'
  ))

rgis_dzc$CL2 <- factor(rgis_dzc$CLUSTER)

clus_dzc <- shape %>% 
  left_join(rgis_dzc,by=c('ID_MUNICIP'='LOC_ID')) 

ggplot() + 
  geom_sf(data=shape, color='grey20', size=0.1, fill=NA) +
  geom_sf(data = subset(clus_dzc, type=='high'),aes(fill=CL2),size=0.1,color='black')  +
  scale_fill_discrete_sequential(palette='Burg',name="High-risk",rev=FALSE,
                                 guide = guide_legend(ncol=2)) +
  theme_void() 

ggplot() + 
  geom_sf(data=shape, color='grey20', size=0.1, fill=NA) +
  geom_sf(data = subset(clus_dzc, type=='low'),aes(fill=CL2), size=0.1, color='black') +
  scale_fill_discrete_sequential(palette='Blues',name="Low-risk",rev=FALSE,
                                 guide = guide_legend(ncol=2)) +
  theme_void() 


# Microcephaly ------------------------------------------------------------

ARQ <- 'microPE.cas'

write.table(microPE[,c(1,2)],file=ARQ , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(microPE[,c(1,5,3)],file='PE.pop' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")
write.table(centroids[,c(1,3,2)],file='PE.geo' , 
            row.names = FALSE,col.names = FALSE,qmethod = "double",fileEncoding = "latin1")

ss.options(list(CaseFile=ARQ,
                StartDate="2015/01/01",
                EndDate="2017/12/31", 
                PrecisionCaseTimes=0, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic     
                PopulationFile="PE.pop",
                CoordinatesFile="PE.geo", 
                CoordinatesType=1, # 0=Cartesian, 1=latitude/longitude
                AnalysisType=1, # 1=Purely Spatial, 3=Retrospective Space-Time
                ModelType=0, # 0=Discrete Poisson
                ScanAreas=3, # 1=High, 2=Low, 3=High and Low
                TimeAggregationUnits=1, # 0=None, 1=Year, 2=Month, 3=Day, 4=Generic
                TimeAggregationLength=1,
                MaxSpatialSizeInPopulationAtRisk=50,
                MinimumTemporalClusterSize=1,
                MaxTemporalSize=1,
                MinimumCasesInHighRateClusters=5,
                MaxSpatialSizeInPopulationAtRisk_Reported=50,
                UseDistanceFromCenterOption_Reported='y',
                MaxSpatialSizeInDistanceFromCenter_Reported=50 # max size radius cluster in km
                
))


ss.options(c("NonCompactnessPenalty=0", "ReportGiniClusters=n", "LogRunToHistoryFile=n"))

modelo <- paste0(ARQ,'_modelo')
write.ss.prm(tmp,modelo)
result_microPE <-  satscan(tmp,modelo, sslocation="/home/lais/SaTScan/")

summary(result_microPE)
result_microPE


### Results

rgis_micro <- result_microPE$gis %>% 
  filter(P_VALUE < 0.05) %>% 
  mutate(type = case_when(
    CLU_RR > 1 ~ 'high',
    CLU_RR < 1 ~ 'low'
  ))

rgis_micro$CL2 <- factor(rgis_micro$CLUSTER)

clus_micro <- shape %>% 
  left_join(rgis_micro,by=c('ID_MUNICIP'='LOC_ID')) 

cores.clu.hr <- sequential_hcl(11, palette = "Burg")
cores.clu.lr <- sequential_hcl(15, palette = 'Blues')

ggplot() + 
  geom_sf(data=shape, color='grey20', size=0.1, fill=NA) +
  geom_sf(data = subset(clus_micro, type=='low'),aes(fill=CL2), lwd=0.1, color='black') +
  scale_fill_manual(values=cores.clu.lr[c(1,4)],name="Low-risk") +
  theme_void() 

ggplot() + 
  geom_sf(data=shape, color='grey20', size=0.1, fill=NA) +
  geom_sf(data = subset(clus_micro, type=='high'),aes(fill=CL2),lwd=0.1, color='black')  +
  scale_fill_manual(values=cores.clu.hr[c(1,4,6)],name="High-risk") +
  theme_void() 


# Zika burden classification ----------------------------------------------

## Organizing cluster analysis results to merge with the spatial models results

rgis_micro <- result_microPE$gis %>% 
  filter(P_VALUE < 0.05) %>% 
  mutate(micro_cluster_type = ifelse(CLU_RR>1, 'high','low'))
micro_cluster <- rgis_micro %>% 
  select(LOC_ID,micro_cluster_likelihood = CLUSTER,micro_cluster_type)

rgis_dzc <- result_dzcPE$gis %>% 
  filter(P_VALUE < 0.05) %>% 
  mutate(dzc_cluster_type = ifelse(CLU_RR>1, 'high','low'))
dzc_cluster <- rgis_dzc %>% 
  select(LOC_ID,dzc_cluster_likelihood = CLUSTER,dzc_cluster_type)


PE <- dzc_cluster %>% 
  full_join(micro_cluster)

PE2 <- shape %>% 
  left_join(as.data.frame(PE), by=c('ID_MUNICIP'='LOC_ID')) %>% 
  mutate(
      micro_dzc_cluster = case_when(
      micro_cluster_type == 'high' & dzc_cluster_type == 'high' ~ '5',
      micro_cluster_type == 'high' & dzc_cluster_type == 'low' ~ '3',
      micro_cluster_type == 'low' & dzc_cluster_type == 'high' ~ '1',
      micro_cluster_type == 'low' & dzc_cluster_type == 'low' ~ '-3',
      is.na(micro_cluster_type) & dzc_cluster_type == 'high' ~ '2',
      micro_cluster_type == 'high' & is.na(dzc_cluster_type) ~ '4',
      is.na(micro_cluster_type) & dzc_cluster_type == 'low' ~ '-1',
      micro_cluster_type == 'low' & is.na(dzc_cluster_type) ~ '-2',
      TRUE  ~ '0')
  )

PE2$micro_dzc_cluster <- factor(PE2$micro_dzc_cluster, levels= c(5,4,3,2,1,0,-1,-2,-3))


cores <- sequential_hcl(11, palette = "Heat")
cores2 <- sequential_hcl(10, palette = "Blues")
posi <- c(1,3,5,7,9)
coresc <- c(cores[posi],'white',cores2[c(9,7,4)]) 

ggplot() +
  geom_sf(data=PE2, aes(fill=micro_dzc_cluster),size=0.1) +
  scale_fill_manual(values=coresc,name="Zika burden\ngradient") +
  theme_void() +
  theme(text = element_text(size=12))

