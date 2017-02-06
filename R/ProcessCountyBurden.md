---
title: "CDC WONDER county estimates for ITHIM"
output: html_document
---
The script creates county-level estimates of disesae burden, which can be used to update an ITHIM object.  It does so by estimating annual number of age-sex-cause-specific deaths and combining them with the Global Burden of Disease numbers for the US.  This script reads in text files output by CDC's WONDER database for deatiled mortality [https://wonder.cdc.gov/ucd-icd10.html]. We use the following causes of death for ITHIM from the ICD-10 113 Cause List. 

_Malignant neoplasms of colon, rectum and anus (C18-C21), Malignant neoplasms of trachea, bronchus and  lung (C33-C34), Malignant neoplasm of breast (C50), #Diabetes mellitus (E10-E14), #Alzheimer's disease (G30), Ischemic heart   diseases (I20-I25), Acute and subacute endocarditis (I33), Diseases of pericardium and acute myocarditis (I30-I31,I40),	Essential hypertension and hypertensive renal disease (I10,I12,I15), #Influenza and pneumonia (J09-J18), Other acute lower	respiratory infections (J20-J22,U04), #Chronic lower respiratory diseases (J40-J47), Other diseases of respiratory system	(J00-J06,J30- J39,J67,J70-J98), Motor vehicle accidents	(V02-V04,V09.0,V09.2,V12-V14,V19.0-V19.2,V19.4-V19.6,V20-V79,V80.3-V80.5,V81.0-V81.1,V82.0-V82.1,V83-V86,V87.0-V87.8,V88.0-V88.8, V89.0,V89.2),	Other land transport accidents (V01,V05-V06,V09.1,V09.3-V09.9,V10-V11,V15-V18,V19.3,V19.8-V19.9,V80.0-V80.2,V80.6-V80.9,V81.2-V81.9,V82.2-V82.9,V87.9,V88.9,	V89.1,V89.3,V89.9)_

Obtaining the full detailed mortality needed for ITHIM requires querying CDC WONDER no less than 4 times. You can run the queries and create these files for your self at [this website](https://wonder.cdc.gov/controller/saved/D76/D8F562).  Be sure to select your state and modify the filename and query properties appropriately between each one. 

1. A query by County (all counties in a state), Age (5 yr groups <1 - 84), Gender (both M and F), and ICD-10 113 Cause (see above) for years 99-2015 - named "ST_counties_99-14_5yr" (ie WI_counties_99-15_5yr.txt)
2. A query by County (all counties in a state), Age (5 yr groups <1 - 84), Gender (both M and F), and ICD-10 113 Cause (see above) for 2014 - named "ST_counties_2015_5yr" (ie WI_counties_2015_5yr.txt)
3. A query by County (all counties in a state), Age (10 yr groups 85+), Gender (both M and F), and ICD-10 113 Cause (see above) for years 99-2015 - named "ST_counties_99-15_10yr" (ie WI_counties_99-15_10yr.txt)
4. A query by County (all counties in a state), Age (10 yr groups 85+), Gender (both M and F), and ICD-10 113 Cause (see above) for 2014 - named "ST_counties_2015_10yr" (ie WI_counties_2015_10yr.txt)

The WONDER query must be divided into the 5yr and 10 year groupings because WONDER does not provide any population counts for age groups over 84 years unless they are all grouped into the '85+' group which is only available using WONDER's 10yr grouping. We use the 5 year groups so that we are able to reconstruct the same age groups used in previous versions of ITHIM. 

The WONDER queries vary in length and contain metadate at the end of the text file. This function will parse the text files downloaded directly from WONDER and prepare them for use in the rest of this script. 
```{r}
#function to read in and clean up each WONDER file
readWONDER <- function(file){
  
  foo<-readLines(file)
  foo<-foo[-1]
  
  if(length(grep("---",foo))>0){
    foo<-foo[-c(grep("---",foo):length(foo))]
  }
  n <- as.numeric(length(foo))
    
  wonder <- read.table(file, sep="\t", as.is = T, nrows = n, skip=1)
  header <- read.table(file, sep="\t", as.is = T, nrows = 1)
  colnames(wonder) <- unlist(header)
  
  wonder$Deaths <- ifelse(wonder$Population== "Suppressed" | wonder$Population== "Missing" , 0, wonder$Deaths)
  wonder$POPclean <- ifelse(wonder$Population== "Suppressed", sample(1:9), as.integer(as.character(wonder$Population))) 
  wonder$POPclean <- ifelse(wonder$Population== "Missing", sample(1:9), as.integer(as.character(wonder$POPclean))) 
  
  return(wonder)
  
}

```

This function will be used at the end of the script to prepare the data for import into an ITHIM object later. It converts the dataframe to a list and aggeregates seveal of the disease outcomes into one 'CVD' category. It is very similar to the readGBD function in the ITHIM package, but accepts a DF as its input instead of an RDS file.

```{R}

#function to isolate single county for import to ITHIM object
#reads in file of all counties GBD and creates a list of GBD tables for each County
readGBDcounties <- function(stateGBD){
  
  stateList <- split(stateGBD, stateGBD$county)
  
  STcounties <- lapply(stateList, function(x) {
    
        gbdList <-split(x, x$disease)
        gbdList[["CVD"]] <- data.frame(disease = "CVD", gbdList$IHD[, 
                                   c("sex", "ageClass")], gbdList$IHD[, c("dproj", "yll", "yld", "daly")] + 
                                   gbdList$InflammatoryHD[, c("dproj", "yll", "yld", "daly")] + 
                                   gbdList$HHD[, c("dproj", "yll", "yld", "daly")])
                                   gbdList2 <- lapply(gbdList, function(x) split(x, as.factor(x$sex)))
                                   gbdList2 <- lapply(gbdList2, function(x) list(M = x$M, F = x$F))
        return(gbdList2)
  })
  
  return(STcounties)
  
}
```
```{r}
library(ggplot2)
library(plyr)


#Identify the state as 2-letter abbreviation (all caps) for which you are running the script here. 
state <- "OR"

setwd(paste("~/Box Sync/work/CoBenefits/MaggieShare/VitalStats/CDCwonder/Counties/",state,sep=""))
  
fiveYearAllYears <- paste("./",state,"_counties_99-15_5yr.txt", sep="")
fiveYearOneYear <-  paste("./",state,"_counties_2015_5yr.txt", sep="")
tenYearAllYears <-  paste("./",state,"_counties_99-15_10yr.txt", sep="")
tenYearOneYear <-  paste("./",state,"_counties_2015_10yr.txt", sep="")
 
#read in data for ages 0-84 (AKA 5yr age groups)
wonder <- readWONDER(fiveYearAllYears)
wonder2 <- readWONDER(fiveYearOneYear)
wonder2$POP2015 <- as.integer(as.character(wonder2$POPclean))
wonder <- cbind(wonder, wonder2["POP2015"])

#read in data for ages 85plus (AKA 10yr age groups)
wonder10 <- readWONDER(tenYearAllYears)
wonder210 <- readWONDER(tenYearOneYear)
wonder210$POP2015 <- as.integer(as.character(wonder210$POPclean))
wonder10 <- cbind(wonder10, wonder210["POP2015"])
names(wonder10) <- names(wonder)

wonderall <- rbind(wonder, wonder10)

```

The code binds the data from the 5yr and 10yr age groups to give you all years and includes the age-gender-cause-specific death counts and rate for all years. A column of single-year population is also binded to the table to give you an estimate of the number of people in each age-gender group for that county. Next, the code imputes 'Surpessed' values (counts less than 10) for the all year deaths totals. County-age-gender-cause cells with less than 10 deaths in ~15years  will remian nearly absent from the final data. Though rare, the single year population may also need to be imputed if less than 10 (see above).

Using the imputed death tallies and all year popualtions we calculate a crude rate for each county-age-gener-cause combo, and apply that rate to the single year population to estimate the single year death count. 

```{r}

wonderall$EstImpDeaths <- as.integer(ifelse(wonderall$Deaths == "Suppressed", sample(0:9), as.integer(as.character(wonderall$Deaths))))
wonderall$EstRate <- wonderall$EstImpDeaths/as.integer(as.character(wonderall$Population))
##### if pop is small set rate to zero or a default
wonderall$EstDeaths <- wonderall$EstRate*wonderall$POP2015

wonderall <- wonderall[,c(2,3,4,6,8,14:17)]
names(wonderall) <- c("county", "countyFIPS", "wonderAge", "wonderGender","wonderCause","POP2015","EstImpDeaths","EstRate","EstDeaths")   

```
One cause in ITHIM is not included in the ICD-10, depression. Here a frame of county-age-gender-depression cells are created, mereged with the appropriate single year populations and merged with the data from WONDER before final application of the burden estimates. 

```{r}

# remove 'Not stated' Age entries
wonderall <- subset(wonderall, wonderAge !="NS")

#depression was not in the WONDER data, since there are no deaths, have to add it manually
depress <- expand.grid(county=unique(wonderall$county), countyFIPS=0, wonderAge=unique(wonderall$wonderAge), wonderGender=unique(wonderall$wonderGender), 
                       wonderCause = "Depression", POP2015= 0, EstImpDeaths= 0, EstRate=0, EstDeaths = 0, stringsAsFactors=T)

popSummary <- ddply(wonderall, .(county, wonderGender, wonderAge), summarise, 
                        depressPOP= as.integer(mean(as.numeric(as.character(POP2015)), na.rm=T)))

depress <- merge(depress, popSummary)

depress$POP2015 <- depress$depressPOP

depress <- depress[,c(1:9)]

wonderall <- rbind(wonderall, depress)

```

In order to use this WONDER data with the Global Burden of Disease data from ITHIM, we reclassify the age, gender, and causes to match. This is dependent on the structure of the WONDER query matching ours. 

```{r}

#recode cause of death for ITHIM
ITHIMcauseKey <- c("ColonCancer","LungCancer","BreastCancer","Diabetes","Dementia","IHD","InflammatoryHD","InflammatoryHD","HHD",
                   "Stroke","AcuteRespInfections","AcuteRespInfections","RespiratoryDisease","RespiratoryDisease","RTIs","RTIs", "Depression")
names(ITHIMcauseKey) <- unique(wonderall$wonderCause)
wonderall$ITHIMcause <- ITHIMcauseKey[as.character(wonderall$wonderCause)]

#recode aage for ITHIM
ITHIMageKey <- c("ageClass1", "ageClass1", "ageClass2", "ageClass2", "ageClass3", "ageClass3", "ageClass3", "ageClass4", "ageClass4", "ageClass4",
                 "ageClass5", "ageClass5", "ageClass5", "ageClass6", "ageClass6", "ageClass7", "ageClass7", "ageClass8", "ageClass8")
names(ITHIMageKey) <- unique(wonderall$wonderAge)
wonderall$ITHIMage <- ITHIMageKey[as.character(wonderall$wonderAge)]

sum(wonderall$POP2015)/16 # this can be used to check the waorking state population 

```
Finally, this table can be paired with the US GBD numbers to estimate a scaled county burden for each age-gender-cause combination, and this info can be used in ITHIM. 
```{r}

########## this 'wonderall' is a table you can crunch for anything ##################
#summarize by state>Urbanization>Gender<age<cause
wonderSummary <- ddply(wonderall, .(county, wonderGender, ITHIMage, ITHIMcause, wonderCause ), summarise, 
                       Deaths = sum(as.integer(as.character(EstDeaths)), na.rm=T),
                       POP = sum(POP2015, na.rm=T)
)

wonderSummary <- ddply(wonderSummary, .(county, wonderGender, ITHIMage, ITHIMcause), summarise, 
                       Deaths = sum(as.integer(as.character(Deaths)), na.rm=T),
                       POP = mean(POP, na.rm=T)
)

names(wonderSummary) <- c("county","gender","ITHIMage","ITHIMcause","deaths","POP")

USGBD <- read.csv(url("https://uwmadison.box.com/shared/static/3k579xkxoyjgdbn6odx15wkxmbskanq9.csv"), header=T)
USGBD$ITHIMage <- as.character(USGBD$ITHIMage)
USGBD$ITHIMcause <- as.character(USGBD$ITHIMcause)

burden <- merge(wonderSummary, USGBD, all.x=T)

#magic number here for setting RRatio to 1
burden$ratio <- ifelse(burden$GBDdeaths < 10, 1, ((burden$deaths/burden$POP)/(burden$GBDdeaths/burden$GBDpop)))
burden$YLL <- burden$ratio * (burden$GBDyll/burden$GBDpop) * burden$POP
burden$YLD <- burden$ratio * (burden$GBDyld/burden$GBDpop) * burden$POP
burden$DALY <- burden$ratio * (burden$GBDdaly/burden$GBDpop) * burden$POP

burden$sex <- factor(ifelse(burden$gender == "Male" | burden$gender == "Male ","M","F"), levels=c("M","F"))
    
burdenOut <- burden[,c("county","ITHIMcause","sex","ITHIMage","deaths","YLL","YLD","DALY")]
names(burdenOut) <- c("county","disease","sex","ageClass","dproj","yll","yld","daly")

burdenOut <- burdenOut[order(burdenOut$county, burdenOut$disease, burdenOut$sex, burdenOut$ageClass),]


#output <- readGBDcounties(burdenOut)

write.csv(output, paste("./",state,"CountiesBurdens2015.csv", sep=""), row.names=F)

```

We can check the distribution of burden by county-age-gender for particular cause (say, Breast Cancer) in all counties by plotting it. 

```{r}

gah <- subset(burdenOut, disease == "BreastCancer")

ggplot(gah, aes(x=ageClass, y=yll, fill=factor(sex))) + geom_bar(stat="identity", position="dodge") + facet_wrap(~county)

```
