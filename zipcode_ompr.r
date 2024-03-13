library(tidyverse)
library(zipcodeR)
library(ompr)
library(ompr.roi)
library(ompr.highs)
#library(ROI.plugin.glpk)
#library(ROI.plugin.symphony)
library(leaflet)

rm(list=ls())
gc(T,T,T)

options(scipen=9999)

highs_solver_parameters <- list(mip_rel_gap=.05,log_to_console=TRUE)



array_2d_multiplication_fcn <- function(static_array, row_variable, column_variable){
  vapply(seq_along(static_array), function(k) static_array[
    row_variable[k], column_variable[k]], numeric(1L))  }

sector_zips <- zipcodeR::zip_code_db %>% dplyr::filter(state=='MN')%>%
  dplyr::select(zipcode,state,lat,lng,population) %>%
  dplyr::mutate(zip4=substr(zipcode,1,4)
                ) %>%
  dplyr::filter(!is.na(lat),!is.na(lng),!is.nan(lat),!is.nan(lng)) %>%
  dplyr::filter(population>0)

#choose level at which to summarize zips: zip5 or zip4
zips_summarized <- sector_zips %>% 
  #dplyr::group_by(zip4) %>%
  dplyr::group_by(zipcode) %>%
  dplyr::summarise(lat=sum(lat*population)/sum(population),
                   lng=sum(lng*population)/sum(population),
                   population=sum(population)) %>%
  #dplyr::rename(zip=zip4) 
  dplyr::rename(zip=zipcode)

zipcount <- nrow(zips_summarized)
sector_count <- 30
bigM <- sum(zips_summarized$population) * 1.1
avg_sector_popn <- sum(zips_summarized$population)/sector_count
avg_sector_popn_pct_tolerance <- 0.1

#distance matrix, weighted by population, for objective function
zip_zip_dist <- geosphere::distm(cbind(zips_summarized$lng,
                                         zips_summarized$lat)) %>%
  measurements::conv_unit('m','mi') * 0.0001 #scale down for numeric performance in solver



mymodel <- ompr::MILPModel() %>%
  add_variable(zip_zip_assign[orig,dest],orig=1:zipcount,dest=1:zipcount,
               type='continuous',lb=0) %>%
  add_variable(zipdest_popn[dest], dest=1:zipcount,type='continuous',lb=0) %>%
  add_variable(zipdest_open[dest],dest=1:zipcount,type='binary') %>%
  
#constraint: for each orig zip, sum of population equals population
  add_constraint(
    sum_expr(zip_zip_assign[orig,dest],dest=1:zipcount) == zips_summarized$population[orig],
    orig = 1:zipcount) %>%

#constraint: for each dest zip, relate zipdest_popn to sum of assigned popn
  add_constraint(
    sum_expr(zip_zip_assign[orig,dest],orig=1:zipcount) == zipdest_popn[dest],
    dest = 1:zipcount) %>%
    
  #constraint: activate zipdest_open if zipdest_popn greater than zero
  add_constraint(zipdest_popn[dest] <= bigM*zipdest_open[dest], dest=1:zipcount) %>%
  
  #constraint: must open exactly sector_count destinations
  add_constraint(sum_expr(zipdest_open[dest],dest=1:zipcount) == sector_count) %>%

#constraint: if opening sector, constrain popn to between 75% and 125% average

add_constraint(zipdest_popn[dest] >= (1-avg_sector_popn_pct_tolerance) * 
                 avg_sector_popn * zipdest_open[dest], dest=1:zipcount) %>%
add_constraint(zipdest_popn[dest] <= (1+avg_sector_popn_pct_tolerance) * 
                 avg_sector_popn * zipdest_open[dest], dest=1:zipcount) %>%
  
  
#objective function: minimize zipzipdist * assigned volume

  set_objective(sum_expr( ompr::colwise(
    array_2d_multiplication_fcn(static_array=zip_zip_dist,
row_variable=orig,column_variable=dest)) *
  zip_zip_assign[orig,dest], orig=1:zipcount,dest=1:zipcount),sense='min')

solve_out <-  mymodel %>% solve_model(highs_optimizer(control=highs_solver_parameters))


 #retrieve solutions
zip_zip_assign_soln <- get_solution(solve_out , zip_zip_assign[orig,dest]) %>%
#assign each origin to only a single destination. winner takes all.
dplyr::filter(value>0) %>%
dplyr::group_by(orig) %>%
 dplyr::mutate(pct_popn = value/sum(value),
               count_dest = length(unique(dest)),
               total_value = sum(value)
                ) #%>%
  #for zips aligned to multiple sectors, randomly select (equal weighting)
  ##dplyr::arrange(orig,desc(pct_popn)) %>%
  ##dplyr::mutate(cumpop = cumsum(pct_popn)) %>%
  #dplyr::mutate(runif = runif(n=1)) %>%
  #dplyr::mutate(cumpop_gte_runif = cumpop>=runif) %>%
  #dplyr::filter(cumpop_gte_runif == TRUE) %>%
  #dplyr::group_by(orig) %>%
  #dplyr::mutate(orig_rownum = dplyr::row_number()) %>%
  #dplyr::filter(orig_rownum == 1) %>%
  #dplyr::select(variable,orig,dest,value=total_value)
  #dplyr::mutate(highest_pct = pct_popn ==max(pct_popn)) %>%
  #dplyr::mutate(value_new = total_value * highest_pct) %>%
  #dplyr::select(variable,orig,dest,value=value_new)

#check by origin to ensure all assigned
origin_soln <- zip_zip_assign_soln %>% dplyr::group_by(orig )%>%
  dplyr::summarise(solver_popn=sum(value)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(popn = zips_summarized$population[orig]) %>%
  dplyr::mutate(absdiff = abs(solver_popn-popn))
summary(origin_soln$absdiff)

dest_soln  <- zip_zip_assign_soln %>% dplyr::group_by(dest )%>%
  dplyr::summarise(solver_popn=sum(value)) %>%
  dplyr::filter(solver_popn>0.1) %>%
  dplyr::mutate(pct = solver_popn / avg_sector_popn)

#check if any dest_soln are outside the allotted bounds. May happen when
#reassigning ZIPs to a single destination with winner-takes-all
summary(dest_soln$pct)


dest_open_soln <- get_solution(solve_out, zipdest_open[dest]) %>%
  dplyr::mutate(value=round(value)) %>%  dplyr::filter(value==1)

#mapping

colors <- RColorBrewer::brewer.pal(n=sector_count,"Set1")
colors <- rep_len(colors,sector_count)
colordf <- cbind.data.frame(dest = 1:sector_count,colors,stringsAsFactors=FALSE)

zip_df <- zip_zip_assign_soln %>% dplyr::left_join(colordf) %>%
  dplyr::left_join(zips_summarized %>% 
                     dplyr::ungroup() %>%
                     dplyr::mutate(orig = dplyr::row_number())) %>%
  dplyr::select(orig, zip,lat,lng,population,dest,colors)
  
map_ <- leaflet(zip_df) %>% 
  addTiles() %>%
  addCircles(lng=~lng,lat=~lat,color=~colors,radius=5)
map_

p1 <- ggplot2::ggplot(data=zip_df,aes(x=lng,y=lat,color=as.factor(dest)))+geom_point() + facet_grid(~dest)
p1
