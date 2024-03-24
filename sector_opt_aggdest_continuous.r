library(tidyverse)
library(ompr)
library(ompr.roi)
library(ompr.highs)
library(leaflet)
library(ggforce)
library(sf)


rm(list=ls())
gc(T,T,T)

options(scipen=9999)

array_2d_multiplication_fcn <- function(static_array, row_variable, column_variable){ 
  vapply(seq_along(static_array), function(k) static_array[
    row_variable[k], column_variable[k]], numeric(1L))  }

highs_solver_parameters <- list(mip_rel_gap=0.1,log_to_console=TRUE)
avg_sector_metric_pct_tolerance <- 0.10  #set tighter than the actual tolerance. Need some wiggle room
#for zips that are assigned to multiple sectors


sector_count <- 48 #48 prod recruiters + 12 SNCOICs in RS TWC


#import data from MCRC structure workbook
all_zips <- readxl::read_excel("C:/Users/r_j_a/Documents/USMC/MCRC Structure/FY20 - FY35 MCRC Structure Data for Distribution (Feb 2024).xlsx", 
                                  sheet = "MCRC Zip DATA")

sector_zips <- all_zips %>%
  janitor::clean_names() %>%
  dplyr::select(rs,rss,state,county,fips,zip,lat,long,metric=qep_2025) %>%
  dplyr::filter(rs == 'TWIN CITIES') %>%
  dplyr::mutate(zip4=substr(zip,1,4),
                  zip3=substr(zip,1,3)
               ) %>%
  dplyr::filter(!is.na(lat),!is.na(long),!is.nan(lat),!is.nan(long)) %>%
  dplyr::filter(metric>0) %>%
  dplyr::mutate(zip = as.character(zip))



#sector_zips <- zipcodeR::zip_code_db %>% dplyr::filter(state=='MN')%>%
#  dplyr::select(zipcode,state,lat,lng,population) %>%
#  dplyr::mutate(zip4=substr(zipcode,1,4),
#                zip3=substr(zipcode,1,3)
 #               ) %>%
#  dplyr::filter(!is.na(lat),!is.na(lng),!is.nan(lat),!is.nan(lng)) %>%
#  dplyr::filter(population>0)

#summarize destinations at zip4 level
dest_summarized <- sector_zips %>% 
  dplyr::group_by(zip4) %>%
  dplyr::summarise(lat=sum(lat*metric)/sum(metric),
                   long=sum(long*metric)/sum(metric)) %>%
  dplyr::rename(zip=zip4) 
  
zipcount <- nrow(sector_zips)
destcount <- nrow(dest_summarized)


bigM <- sum(sector_zips$metric) * 1.1
avg_sector_metric <- sum(sector_zips$metric)/sector_count


#distance matrix, weighted by metric, for objective function
zip_dest_dist <- geosphere::distm(
  x=cbind(sector_zips$long,sector_zips$lat),
  y=cbind(dest_summarized$long,dest_summarized$lat)) %>%
  measurements::conv_unit('m','mi') * 0.0001 #scale down for numeric performance in solver



mymodel <- ompr::MILPModel() %>%
  add_variable(zip_dest_assign[orig,dest],orig=1:zipcount,dest=1:destcount,
               type='continuous',lb=0) %>%
  add_variable(zipdest_popn[dest], dest=1:destcount,type='continuous',lb=0) %>%
  add_variable(zipdest_open[dest],dest=1:destcount,type='binary') %>%
  
#constraint: for each orig zip, sum of metric equals metric
  add_constraint(
    sum_expr(zip_dest_assign[orig,dest],dest=1:destcount) == sector_zips$metric[orig],
    orig = 1:zipcount) %>%

#constraint: for each dest zip, relate zipdest_popn to sum of assigned popn
  add_constraint(
    sum_expr(zip_dest_assign[orig,dest],orig=1:zipcount) == zipdest_popn[dest],
    dest = 1:destcount) %>%
    
  #constraint: activate zipdest_open if zipdest_popn greater than zero
  add_constraint(zipdest_popn[dest] <= bigM*zipdest_open[dest], dest=1:destcount) %>%
  
  #constraint: must open exactly sector_count destinations
  add_constraint(sum_expr(zipdest_open[dest],dest=1:destcount) == sector_count) %>%

#constraint: if opening sector, constrain popn to within X% of average

add_constraint(zipdest_popn[dest] >= (1-avg_sector_metric_pct_tolerance) * 
                 avg_sector_metric * zipdest_open[dest], dest=1:destcount) %>%
add_constraint(zipdest_popn[dest] <= (1+avg_sector_metric_pct_tolerance) * 
                 avg_sector_metric * zipdest_open[dest], dest=1:destcount) %>%
  
  
#objective function: minimize zipzipdist * assigned volume

  set_objective(sum_expr( ompr::colwise(
    array_2d_multiplication_fcn(static_array=zip_dest_dist,
row_variable=orig,column_variable=dest)) *
  zip_dest_assign[orig,dest], orig=1:zipcount,dest=1:destcount),sense='min')

solve_out <-  mymodel %>% solve_model(highs_optimizer(control=highs_solver_parameters))


 #retrieve solutions
zip_dest_assign_soln <- get_solution(solve_out , zip_dest_assign[orig,dest]) %>%
#assign each origin to only a single destination. winner takes all.
dplyr::filter(value>0) %>%
dplyr::group_by(orig) %>%
 dplyr::mutate(pct_popn = value/sum(value),
               count_dest = length(unique(dest)),
               total_value = sum(value) ) # %>%

# test not using winner takes all logic
  #winner takes all logic for zips assigned to multiple sectors
  #dplyr::mutate(max_pct = pct_popn == max(pct_popn)) %>%
  #dplyr::filter(max_pct == TRUE) %>%
  #dplyr::mutate(value = total_value) %>%
  #dplyr::select(-total_value)
  

#check by origin to ensure all assigned
origin_soln <- zip_dest_assign_soln %>% dplyr::group_by(orig )%>%
  dplyr::summarise(solver_metric=sum(value)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(metric = sector_zips$metric[orig]) %>%
  dplyr::mutate(absdiff = abs(solver_metric-metric))

summary(origin_soln$absdiff)

dest_soln  <- zip_dest_assign_soln %>% dplyr::group_by(dest )%>%
  dplyr::summarise(solver_metric=sum(value)) %>%
  dplyr::filter(solver_metric>0.01) %>%
  dplyr::mutate(pct = solver_metric / avg_sector_metric)

#check if any dest_soln are outside the allotted bounds. May happen when
#reassigning ZIPs to a single destination with winner-takes-all
summary(dest_soln$pct)


dest_open_soln <- get_solution(solve_out, zipdest_open[dest]) %>%
  dplyr::mutate(value=round(value)) %>%  dplyr::filter(value==1)

#mapping

#colors <- RColorBrewer::brewer.pal(n=sector_count,"Set1")
colors <- c('orange','blue','green','red','yellow','purple')
colors <- rep_len(colors,sector_count)
colordf <- cbind.data.frame(dest = unique(zip_dest_assign_soln$dest),colors,stringsAsFactors=FALSE)

zip_df <- zip_dest_assign_soln %>% 
  #since we are allowing a single zip to be aligned to multiple sectors for now,
  #limit to just a single row in zip_df and signify multi-alignment by black marker
  dplyr::group_by(orig) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::left_join(colordf) %>%
  dplyr::mutate(colors = dplyr::if_else(count_dest>1,"black",colors)) %>%
  dplyr::left_join(sector_zips %>% 
                     dplyr::ungroup() %>%
                     dplyr::mutate(orig = dplyr::row_number())) %>%
  dplyr::select(orig, zip,lat,long,metric,dest,count_dest,colors)






#choropleth map of all zips

#shapefile from:
#https://catalog.data.gov/dataset/tiger-line-shapefile-2019-2010-nation-u-s-2010-census-5-digit-zip-code-tabulation-area-zcta5-na
# Load the shapes

postal.codes=sf::st_read(
  "C:/Users/r_j_a/Documents/USMC/MCRC Structure/tl_2019_us_zcta510/tl_2019_us_zcta510.shp") %>%
  dplyr::filter(ZCTA5CE10 %in% sector_zips$zip)

postal.codes2 <- postal.codes %>% 
  dplyr::left_join(zip_df %>% 
  #dplyr::mutate(ZCTA5CE10 = as.character(zip)) %>%
    dplyr::ungroup() %>%
    dplyr::select(ZCTA5CE10 = zip, dest, colors)) 

leaflet_all_sectors <- leaflet(postal.codes2) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(color=~colors,weight = 1, smoothFactor = 0.5,
              opacity = 1.0, fillOpacity = 0.5)


#identify any sectors that consist of non-contiguous components

#use sf::st_union to check for disconnectedness
dest_soln$union_class <- NA

for(destindex in 1:sector_count){
  print(destindex)
  sector_postal_codes <- postal.codes %>% 
    filter(ZCTA5CE10 %in% zip_df$zip[zip_df$dest == dest_soln$dest[destindex]])
  
  sector_union <- sf::st_union(sector_postal_codes)
  
  dest_soln$union_class[destindex] <- class(sector_union)[1]
}  

#make leaflet maps for each sector
sector_maps <- list()

for(i in 1:sector_count){
  print(i)
  
  dest_ <- dest_soln$dest[i]
  union_class_ <- dest_soln$union_class[i]
  
  data_ <- postal.codes2 %>% dplyr::filter(dest == dest_)
  
  leaflet_sector <- leaflet(data_) %>% 
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(color=~colors,weight = 1, smoothFactor = 0.5,
                opacity = 1.0, fillOpacity = 0.5) %>%
    addControl(paste(
      "Dest ID: ", dest_,"<br>",
      "Union Class: ",union_class_,"<br>",
      "ZCTA Count: ", nrow(data_),
      sep= " "),position = "topleft") %>%
    #minimap
    leaflet::addMiniMap()
  
  sector_maps[[i]] <- leaflet_sector
  
}

leafsync::latticeview(sector_maps)
