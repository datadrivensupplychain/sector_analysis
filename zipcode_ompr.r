#template of network model for recruiting
library(tidyverse)
library(magrittr)
library(ompr)
library(ompr.roi)
library(ompr.highs)
library(leaflet)
library(leaflet.minicharts)
library(zipcodeR)
library(geosphere)
library(measurements)
library(RColorBrewer)
library(scdesign)

rm(list=ls())
gc(T,T,T)
options(scipen=99999)

array_2d_multiplication_fcn <- function(static_array, row_variable, column_variable){
  vapply(seq_along(static_array), function(k) static_array[
    row_variable[k], column_variable[k]], numeric(1L))  }


rs_zipcodes <- zipcodeR::zip_code_db %>% dplyr::filter(state %in% c('ND','SD','MN')) %>%
  dplyr::filter(!is.na(lat), !is.na(lng)) %>%
  dplyr::select(zipcode,lat,lng,population,population_density,land_area_in_sqmi) %>%
  dplyr::mutate(zip4 = substr(zipcode,1,4)) 

#clustering_count <- length(unique(rs_zipcodes$zip4))
sector_count <- 40


#clustering maybe unnecesesary
clustering_count <- 300

#geocluster zipcodes to generate potential sector centers
zip_clustering <- scdesign::geocluster_kmeans_kmedoids(locationdf=rs_zipcodes %>% dplyr::select(lng,lat),
                                                       cluster_count = clustering_count,cluster_method='kmeans')


rs_zipcodes$cluster <- zip_clustering$resultsdf$ClusterNumber

cluster_centers  <- zip_clustering$resultsdf %>% dplyr::distinct(Cluster_Lng,Cluster_Lat) %>%
  dplyr::rename(lat=Cluster_Lat,lng=Cluster_Lng)

rs_zipcodes_agg <- rs_zipcodes %>%
  
  #test: don't actually grup
  #dplyr::group_by(cluster) %>%
  dplyr::group_by(zip4) %>%
  #dplyr::group_by(zipcode) %>%
  dplyr::summarise(lat=sum(lat*population)/sum(population),lng=sum(lng*population)/sum(population),population=sum(population))

#leaflet map
map <- leaflet(rs_zipcodes_agg) %>% addTiles() %>%
  addCircleMarkers(lng=~lng,lat=~lat,radius=2)
map

zipcount <- nrow(rs_zipcodes_agg)

zipclusterdist <- geosphere::distm(x=cbind(rs_zipcodes_agg$lng,rs_zipcodes_agg$lat),
                                   y=cbind(cluster_centers$lng,cluster_centers$lat)) %>%
  measurements::conv_unit('m','mi')

#zipzipdist <- geosphere::distm(cbind(rs_zipcodes_agg$lng,rs_zipcodes_agg$lat)) %>%
#  measurements::conv_unit('m','mi')

popn_matrix <- matrix(nrow= zipcount, ncol=clustering_count, data= rs_zipcodes_agg$population,byrow=FALSE)

wtddistance_matrix <- zipclusterdist * popn_matrix

avg_popn <- sum(rs_zipcodes_agg$population)/sector_count

#generate potential sector centers

milp_model <- ompr::MILPModel() %>%
  
  #add_variable(zip_zipcenter[ziporigin,zipdest],
  #             ziporigin = 1 :zipcount, zipdest=1:zipcount, type='binary') %>%

  add_variable(zip_cluster[ziporigin,clusterdest],
               ziporigin = 1 :zipcount, clusterdest=1:clustering_count, type='binary') %>%
  
  #open/close for sector centers
  add_variable(center_open[clusterdest], clusterdest=1:clustering_count,type='binary') %>%
  #population within each sector
  
  add_variable(sector_popn[clusterdest], clusterdest=1:clustering_count,type='continuous',lb=0) %>%
  
  #require each zip to be aligned to exactly one center
  
  add_constraint(sum_expr(zip_cluster[ziporigin,clusterdest],clusterdest=1:clustering_count)==1, ziporigin=1:zipcount) %>%
  
  #force opening with bigM
  
  add_constraint(
    sum_expr(zip_cluster[ziporigin,clusterdest],ziporigin = 1:zipcount) <= 
      9999*center_open[clusterdest], clusterdest=1:clustering_count) %>%
  
  #required number of opens: corresponss to number of recruitesr
  add_constraint(sum_expr(center_open[clusterdest], clusterdest=1:clustering_count) == sector_count)

#use constraints to define population for each destination zip
for(i in 1:clustering_count){
milp_model <- milp_model %>%
add_constraint(
sum_expr(ompr::colwise(popn_matrix[ziporigin,i]) * zip_cluster[ziporigin,i], ziporigin=1:zipcount) == sector_popn[i] ) %>%
 # use capacity constraints to keep sector population within +/- x% of average
   add_constraint( sector_popn[i] <= avg_popn * 1.5 * center_open[i]) %>%
   add_constraint( sector_popn[i] >= avg_popn * 0.5 * center_open[i])

}


#add constraints to relate sector to RSS
#add constraints to relate sector to RSS
#add constraints to relate sector to RSS
#add constraints to relate sector to RSS
#add constraints to relate sector to RSS



milp_model <- milp_model %>%
  
  #set objective: minimize distance to center
  set_objective(sum_expr( ompr::colwise(
    array_2d_multiplication_fcn(static_array = wtddistance_matrix, #zipclusterdist,
    row_variable=ziporigin,column_variable=clusterdest)) *
      zip_cluster[ziporigin,clusterdest], 
    ziporigin=1:zipcount,clusterdest=1:clustering_count), sense='min')

highs_solver_parameters <- list(mip_rel_gap=.05,log_to_console=TRUE)

milp_model_out <-  milp_model %>% solve_model(highs_optimizer(control=highs_solver_parameters))
   
  # ompr::solve_model(with_ROI(solver='symphony',verbosity=1))

soln_out <- get_solution(milp_model_out,zip_zipcenter[ziporigin,zipdest]) %>%
  dplyr::mutate(value = round(value)) %>%
  dplyr::filter(value==1) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(zip4origin = rs_zipcodes_agg$zip4[ziporigin],
                zip4dest = rs_zipcodes_agg$zip4[zipdest])

opensoln_out <- get_solution(milp_model_out,zipcenter_open[zipdest]) %>%
  dplyr::mutate(value = round(value)) %>%
  dplyr::filter(value==1) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(zip4dest = rs_zipcodes_agg$zip4[zipdest])

sectorpopn_out <- get_solution(milp_model_out, sector_popn[zipdest]) %>%
  dplyr::mutate(value = plyr::round_any(value,1)) %>%
  dplyr::filter(value>0)

summary(sectorpopn_out$value/avg_popn)

#manually add up population
popn_df <- soln_out %>% dplyr::left_join(rs_zipcodes_agg %>% dplyr::select(zip4origin=zip4,population)) %>%
  dplyr::group_by(zip4dest) %>% dplyr::summarise(popn=sum(population))


summary(popn_df$popn/avg_popn)


sectorzip4 <- unique(soln_out$zip4dest)




colors <- brewer.pal(n=sector_count,"Set1")
#if not enough colors in the palette, repeat
colors <- rep_len(colors,sector_count)
colordf <- cbind.data.frame(zip4center=sectorzip4,colors,stringsAsFactors=FALSE)

rs_zipcodes <- rs_zipcodes %>% dplyr::left_join(
  soln_out %>% dplyr::select(zip4 = zip4origin,zip4center=zip4dest) ) %>%
  dplyr::left_join(colordf)


p2 <- leaflet(rs_zipcodes) %>% 
  addTiles() %>%
  addCircles(lng=~lng,lat=~lat,color=~colors,popup=~as.factor(zip4))
p2
