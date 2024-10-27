
library(tidyverse)

################################
## Note, the first row of the csv file
##   is the variable names
## The second row is the units!
##   The data starts on line 3
## 
## So first read in the variable names only
##   then read in the data.
##
## Also note, the North Atlantic is coded as BASIN "NA"
##   which R reads as an NA - not available...
##
storm_col_names <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.since1980.list.v04r00.csv",
                            n_max = 0)

storms_raw <- read_csv("https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.since1980.list.v04r00.csv",
                       skip = 2, col_names = names(storm_col_names), na="" )

save(storms_raw, file="../data/rawIBTACS.RData")
load("../data/rawIBTACS.RData")
## The current year is still happening. So all data before this year!

current_year <- year(Sys.time())

storms_fct <- storms_raw %>%
  dplyr::select(SID, SEASON, BASIN, NAME, NATURE, WMO_WIND, USA_WIND, TOKYO_WIND, CMA_WIND, HKO_WIND, NEWDELHI_WIND, REUNION_WIND, BOM_WIND, NADI_WIND, WELLINGTON_WIND, DS824_WIND, TD9636_WIND, NEUMANN_WIND, MLC_WIND) %>%
  mutate(SEASON=as.numeric(SEASON)) %>% 
  filter(SEASON < current_year) %>% 
  group_by(SID) %>% 
  summarize(N=n(),                  ## Number of observations
            SEASON=first(SEASON),   ## Save the Season/year
            BASIN=first(BASIN),     ## Save the BASIN
            NAME=first(NAME),       ## Storm name (might be useful for annotation)
            NATURE=first(NATURE),   ## The storm classification (we create our own later, this is needed for one storm, see below)
                                    ## Then find the max wind speed -- its recorded across differnet variables for different basics
            Wind = max(c(WMO_WIND, USA_WIND, TOKYO_WIND, CMA_WIND, HKO_WIND, NEWDELHI_WIND, REUNION_WIND, BOM_WIND, NADI_WIND, WELLINGTON_WIND, DS824_WIND, TD9636_WIND, NEUMANN_WIND, MLC_WIND), na.rm=TRUE) ) %>%
  mutate(BASIN = factor(BASIN,     ## Label the BASIN with something informative
                        levels=c("NA", "EP", "WP", "NI", "SI", "SP", "SA"),
                        labels=c("North Atlantic",
                                 "Eastern North Pacific",
                                 "Western North Pacific",
                                 "Northern Indian",
                                 "Southern Indian",
                                 "Southern Pacific",
                                 "Southern Atlantic")),
         StormCat = case_when(Wind >= 137 ~ "Cat-5",
                              Wind >= 113 ~ "Cat-4",
                              Wind >= 96  ~ "Cat-3",
                              Wind >= 83  ~ "Cat-2",
                              Wind >= 64  ~ "Cat-1",
                              Wind >= 34  ~ "Cat-0",
                              is.infinite(Wind) & NATURE=="TS" ~ "Cat-0"),  ## A single storm, SID="1980085S11135" is marked as a tropical storm but without wind speeds
         StormCat = factor(StormCat, levels=c("Cat-0", "Cat-1", "Cat-2", "Cat-3", "Cat-4", "Cat-5") ) ) %>%
  dplyr::filter(!is.na(StormCat))

storms_duration <- storms_raw %>%
  dplyr::select(SID, SEASON, BASIN, ISO_TIME, NAME, NATURE, WMO_WIND, USA_WIND, TOKYO_WIND, CMA_WIND, HKO_WIND, NEWDELHI_WIND, REUNION_WIND, BOM_WIND, NADI_WIND, WELLINGTON_WIND, DS824_WIND, TD9636_WIND, NEUMANN_WIND, MLC_WIND) %>%
  mutate(SEASON=as.numeric(SEASON)) %>% 
  filter(SEASON < current_year) %>% 
  mutate(Wind = pmax(WMO_WIND, USA_WIND, TOKYO_WIND, CMA_WIND, HKO_WIND, NEWDELHI_WIND, REUNION_WIND, BOM_WIND, NADI_WIND, WELLINGTON_WIND, DS824_WIND, TD9636_WIND, NEUMANN_WIND, MLC_WIND, na.rm=TRUE) ) %>%
  mutate(Wind = ifelse(is.infinite(Wind), 35, Wind) ) %>%
  filter(Wind > 34) %>%
  group_by(SID) %>% 
  summarize(SID = first(SID),
            SEASON = first(SEASON),
            BASIN = first(BASIN),
            Duration = max(ISO_TIME)-min(ISO_TIME)) %>%
  mutate(Duration_Class = ifelse(Duration <= 48*60*60, "Shorty", "Non_Shorty")) %>%
  mutate(BASIN = factor(BASIN,     ## Label the BASIN with something informative
                        levels=c("NA", "EP", "WP", "NI", "SI", "SP", "SA"),
                        labels=c("North Atlantic",
                                 "Eastern North Pacific",
                                 "Western North Pacific",
                                 "Northern Indian",
                                 "Southern Indian",
                                 "Southern Pacific",
                                 "Southern Atlantic")) )

#########################
## In storm_counts_fct
##
## Each row is a storm from a specific basin of origin (first basin with record)
##   along with its classification

storm_counts <- storms_fct %>%
  group_by(SEASON, BASIN) %>%
  summarize(Counts=n() ) %>%
  filter(BASIN != "Southern Atlantic")


storm_counts_wide <- storm_counts %>%
  pivot_wider(id_cols=SEASON, 
              names_from=BASIN,
              values_from=Counts)

save(storm_counts, storm_counts_wide, file="./data/tropicalStormCounts.RData")



