library(jsonlite)
library(reshape2)

getTractAgeSex <- function(state, county) {
  varString <-
    "B01001_003E,B01001_004E,B01001_005E,B01001_006E,B01001_007E,B01001_008E,B01001_009E,B01001_010E,B01001_011E,B01001_012E,B01001_013E,B01001_014E,B01001_015E,B01001_016E,B01001_017E,B01001_018E,B01001_019E,B01001_020E,B01001_021E,B01001_022E,B01001_023E,B01001_024E,B01001_025E,B01001_027E,B01001_028E,B01001_029E,B01001_030E,B01001_031E,B01001_032E,B01001_033E,B01001_034E,B01001_035E,B01001_036E,B01001_037E,B01001_038E,B01001_039E,B01001_040E,B01001_041E,B01001_042E,B01001_043E,B01001_044E,B01001_045E,B01001_046E,B01001_047E,B01001_048E,B01001_049E"
  
  ACSpop <-
    as.data.frame(fromJSON(
      paste(
        "http://api.census.gov/data/2015/acs5?get=NAME,",
        varString,
        "&for=tract:*&in=state:",
        state,
        "+county:",
        county,
        "&key=f78d6b6c18608edc379b5a06c55407ceb45e7038",
        sep = ""
      )
    ))
  ACSpop <- ACSpop[-1,]
  colnames(ACSpop) <-
    c(
      "name",
      "M_Under 5 years",
      "M_5 to 9 years",
      "M_10 to 14 years",
      "M_15 to 17 years",
      "M_18 and 19 years",
      "M_20 years",
      "M_21 years",
      "M_22 to 24 years",
      "M_25 to 29 years",
      "M_30 to 34 years",
      "M_35 to 39 years",
      "M_40 to 44 years",
      "M_45 to 49 years",
      "M_50 to 54 years",
      "M_55 to 59 years",
      "M_60 and 61 years",
      "M_62 to 64 years",
      "M_65 and 66 years",
      "M_67 to 69 years",
      "M_70 to 74 years",
      "M_75 to 79 years",
      "M_80 to 84 years",
      "M_85 years and over",
      "F_Under 5 years",
      "F_5 to 9 years",
      "F_10 to 14 years",
      "F_15 to 17 years",
      "F_18 and 19 years",
      "F_20 years",
      "F_21 years",
      "F_22 to 24 years",
      "F_25 to 29 years",
      "F_30 to 34 years",
      "F_35 to 39 years",
      "F_40 to 44 years",
      "F_45 to 49 years",
      "F_50 to 54 years",
      "F_55 to 59 years",
      "F_60 and 61 years",
      "F_62 to 64 years",
      "F_65 and 66 years",
      "F_67 to 69 years",
      "F_70 to 74 years",
      "F_75 to 79 years",
      "F_80 to 84 years",
      "F_85 years and over",
      "state",
      "county",
      "tract"
    )
  
  ACSpop <- melt(ACSpop, id = c("name", "state", "county", "tract"))
  
  ACSpop$gender <-
    matrix(unlist(strsplit(as.character(ACSpop$variable), "_")), ncol = 2, byrow =
             T)[, 1]
  ACSpop$acsAge <-
    matrix(unlist(strsplit(as.character(ACSpop$variable), "_")), ncol = 2, byrow =
             T)[, 2]
  
  ITHIMageKey <-
    c(
      "ageClass1",
      "ageClass2",
      "ageClass2",
      "ageClass3",
      "ageClass3",
      "ageClass3",
      "ageClass3",
      "ageClass3",
      "ageClass3",
      "ageClass4",
      "ageClass4",
      "ageClass4",
      "ageClass5",
      "ageClass5",
      "ageClass5",
      "ageClass6",
      "ageClass6",
      "ageClass6",
      "ageClass6",
      "ageClass7",
      "ageClass7",
      "ageClass8",
      "ageClass8"
    )
  names(ITHIMageKey) <- unique(ACSpop$acsAge)
  ACSpop$ITHIMage <- ITHIMageKey[as.character(ACSpop$acsAge)]
  
  WONDERageKey <-
    c(
      "1-4 years",
      "5-14 years",
      "5-14 years",
      "15-24 years",
      "15-24 years",
      "15-24 years",
      "15-24 years",
      "15-24 years",
      "25-34 years",
      "25-34 years",
      "35-44 years",
      "35-44 years",
      "45-54 years",
      "45-54 years",
      "55-64 years",
      "55-64 years",
      "55-64 years",
      "65-74 years",
      "65-74 years",
      "65-74 years",
      "75-84 years",
      "75-84 years",
      "85+"
    )
  names(WONDERageKey) <- unique(ACSpop$acsAge)
  ACSpop$WONDERage <- WONDERageKey[as.character(ACSpop$acsAge)]
  
  return(ACSpop)
}

DallasCnty <- getTractAgeSex(state = 48,county = 113) 
DaneCnty <- getTractAgeSex(state = 55,county = 025) 

library(ggplot2)
library(plyr)

ggplot(DaneCnty, aes(x=tract,y=value, fill=gender)) + geom_bar(stat="identity") + facet_grid(. ~ acsAge) + coord_flip()


DallasWONDERpop <- ddply(DallasCnty, .(state, county, tract, gender, WONDERage), summarise,
                         Population = sum(as.numeric(as.character(value)), na.rm=T)
                                      
                         )

#write.csv(DallasWONDERpop, "./Box Sync/work/Louisville/Dallas/dallasWonderPop.csv", row.names=F)

