

rm(list=ls())

# Top ---------------------------------------------------------------------

set.seed(42)

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/UNH Research/Data/CSV_Files/Movement")

library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(adehabitatLT)
library(ggrepel)
library(gganimate)
library(move)
library(moveVis)
library(viridis)
library(RColorBrewer)
library(lubridate)
library(argosfilter)
library(pvclust)
library(TraMineR)
library(lubridate)
library(ggmosaic)
library(glmm)
library(modelr)
library(geosphere)
library(adehabitatHR)
library(ggmap)
library(ggsn)
library(ggnewscale)
library(ggstar)
library(ggspatial)
library(maptools)
library(lme4)
library(moments)
library(magick)
library(remotes)
library(readr)
library(arm)
#Installing rtools40 here done externally
#remotes::install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F) #Vignette not working
library(RSP)
library(actel)
library(latticeDensity)



theme_set(theme_bw(base_size=25))

'%notin%'<-Negate('%in%')



# Handling time for surgeries ---------------------------------------------
surgery2017<-read.csv("../2017_Just_Sculpin.csv")%>%
  dplyr::select(Date,Anesthesia=Time.under,Surgery.start,Surgery.stop=Surgery.done,Recovery=Recovered,Release=Released)%>%
  filter(Anesthesia!="")%>%
  mutate(Anesthesia=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Anesthesia))),
         Surgery.start=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.start))),
         Surgery.stop=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.stop))),
         Recovery=dmy_hm(paste(Date,Recovery)),
         Release=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Release))))%>%
  filter(!is.na(Anesthesia))%>%dplyr::select(-Date)

surgery2018<-read.csv("../2018_Fish_Data.csv")%>%
  dplyr::select(Date,Anesthesia=Anesthesia.Start,Surgery.start=Tag.Start,Surgery.stop=Tag.Stop,Recovery=Recovery.time,Release=Release.Time)%>%
  filter(Anesthesia!="-")%>%
  mutate(Anesthesia=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Anesthesia))),
         Surgery.start=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.start))),
         Surgery.stop=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.stop))),
         Recovery=dmy_hm(paste(Date,Recovery)),
         Release=dmy_hm(paste(Date,gsub("([a,p])"," \\1",Release))))%>%
  filter(!is.na(Anesthesia))%>%dplyr::select(-Date)

surgery2019<-read.csv("../2019_Sculpin_Tagged.csv")%>%
  dplyr::select(Date,Anesthesia=Anestheisa,Surgery.start=Tag_Start,Surgery.stop=Tag_Stop,Recovery,Release)%>%
  mutate(Anesthesia=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Anesthesia))),
         Surgery.start=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.start))),
         Surgery.stop=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.stop))),
         Recovery=mdy_hm(paste(Date,Recovery)),
         Release=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Release))))%>%
  filter(!is.na(Anesthesia))%>%dplyr::select(-Date)

surgeries<-bind_rows(surgery2017,surgery2018,surgery2019)%>%
  mutate(Surgery.start=if_else(as.numeric(difftime(Surgery.start,Anesthesia,units="mins"))<0,
                               Surgery.start+days(1),Surgery.start),
         Surgery.stop=if_else(as.numeric(difftime(Surgery.stop,Surgery.start,units="mins"))<0,
                              Surgery.stop+days(1),Surgery.stop),
         Recovery=if_else(as.numeric(difftime(Recovery,Surgery.stop,units="mins"))<0,
                          Recovery+days(1),Recovery),
         Release=if_else(as.numeric(difftime(Release,Recovery,units="mins"))<0,
                         Release+days(1),Release))

#Time in anesthesia before surgery starts
mean(as.numeric(difftime(surgeries$Surgery.start,surgeries$Anesthesia,units="mins")),na.rm=T)
sd(as.numeric(difftime(surgeries$Surgery.start,surgeries$Anesthesia,units="mins")),na.rm=T)/sqrt(nrow(filter(surgeries,!is.na(Anesthesia))))
#Time under surgery
mean(as.numeric(difftime(surgeries$Surgery.stop,surgeries$Surgery.start,units="mins")),na.rm=T)
sd(as.numeric(difftime(surgeries$Surgery.stop,surgeries$Surgery.start,units="mins")),na.rm=T)/sqrt(nrow(filter(surgeries,!is.na(Surgery.stop))))
#Time for recovery
mean(as.numeric(difftime(surgeries$Release,surgeries$Recovery,units="mins")),na.rm=T)



# Detection Data ----------------------------------------------------------


tags<-read_csv("tagDetails.csv", 
               col_types = cols(Release_date = col_date(format = "%Y-%m-%d"), 
                                Tag_ID = col_character(), Year = col_character()))%>%
  mutate(Fish_ID=ifelse(Tag_Type=="V9A",as.character(as.numeric(Tag_ID)+0.5),
                        ifelse(Tag_Type=="V9P",as.character(as.numeric(Tag_ID)-0.5),Tag_ID)))

#Everything above this takes a VERY long time, so best to just read in the detections from here
V9detections<-read_csv("final.SculpinDetections.csv",
                       col_types = cols(Date.and.Time = col_datetime(format = "%Y-%m-%d %H:%M:%S"), 
                                        Fish_ID = col_character(), Receiver_ID = col_character(), 
                                        Release_date = col_date(format = "%Y-%m-%d"), 
                                        Sensor.Measure = col_double(), Sensor.Value = col_number(), 
                                        Tag_ID = col_character(), Year = col_character()))


V9s<-filter(tags,Tag_Type=="V9")
nrow(filter(V9s,Year=="2017")) #2017
nrow(filter(V9s,Year=="2018")) #2018
n_distinct(filter(V9s,Year=="2019")$Fish_ID) #2019
#Detected
table(V9s$Recovered)



# Mapping with detection counts -------------------------------------------
spatial<-read.csv("spatial.csv")


#Doing it by distances
tremblay<-st_read("C:/Users/nh1087/OneDrive - USNH/Documents/UNH Research/Data/Mapping Data/tremblay_and_Islands.shp")
tremblay<-fortify(tremblay)

tremblay<-filter(tremblay,group=="5.1")

#Overlaying the detection heatmap on bathymetry
library(sp)
library(sf)
spatialBuffer<-st_as_sf(filter(spatial,Station.name!="Camp"),coords=c("Longitude","Latitude"))
st_crs(spatialBuffer)<-CRS("+proj=longlat +datum=WGS84")  
spatialBuffer <- st_transform(spatialBuffer, CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84"))
spatialBuffer <- st_buffer(spatialBuffer, 400)
spatialBuffer <- as_Spatial(spatialBuffer)
spatialBuffer <- spTransform(spatialBuffer,CRS("+proj=longlat +datum=WGS84"))
plotBuffer<-fortify(spatialBuffer,region="Station.name")
plotBuffer<-left_join(plotBuffer,V9detections%>%group_by(Station.name)%>%summarise(N=n()),by=c("id"="Station.name"))
plotBuffer[is.na(plotBuffer)]<-0

rangeBuffer<-st_as_sf(filter(spatial,Station.name!="Camp"),coords=c("Longitude","Latitude"))
st_crs(rangeBuffer)<-CRS("+proj=longlat +datum=WGS84")  
rangeBuffer <- st_transform(rangeBuffer, CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84"))
rangeBuffer <- st_buffer(rangeBuffer, 260)
rangeBuffer <- as_Spatial(rangeBuffer)
rangeBuffer <- spTransform(rangeBuffer,CRS("+proj=longlat +datum=WGS84"))
rangeBuffer<-fortify(rangeBuffer,region="Station.name")
rangeBuffer<-left_join(rangeBuffer,V9detections%>%group_by(Station.name)%>%summarise(N=n()),by=c("id"="Station.name"))
rangeBuffer[is.na(rangeBuffer)]<-0


bathy<-raster("C:/Users/nh1087/OneDrive - USNH/Documents/UNH Research/Data/Mapping Data/GEBCO_2020_15_Apr_2021_ed9a233c5ed5/gebco_2020_n72.59971618652344_s72.22480773925781_w-81.42654418945312_e-80.47622680664062.nc",varname="elevation")
bathyDF<-as.data.frame(bathy[[1]],xy=T)
colnames(bathyDF)<-c("Long","Lat","z")
bathyDF$z<-ifelse(bathyDF$z>1,NA,bathyDF$z)
bathyDF$absZ<-abs(bathyDF$z)
bathyDF<-filter(bathyDF,z>-270)

bathyPoly<-rasterToPolygons(bathy)
bathyPoly<-fortify(bathyPoly,region="Elevation.relative.to.sea.level")
bathyPoly$id<-as.numeric(bathyPoly$id)
bathyPolyT<-filter(bathyPoly,id<=0 & !(long>-80.95&lat<72.35) & 
                     !(long>-80.7 &lat<72.45) & !long>-80.64) #Keep only the Tremblay area (below 0 and not over in the next area)
bathyPolyW<-filter(bathyPoly,id<=0 & id>-270)
bathyPolyL<-filter(bathyPoly,id>0)


load("~/UNH Research/Data/Mapping Data/googleTR.rda")


water.col <- colorRampPalette(c("purple4", "navy", "blue", "skyblue", "lightblue", "white"))


#Complete detections Map
bathyMap<-image_graph(width=1330,height=1250,res=96)
ggmap(googleTR)+
  geom_contour_filled(data=filter(bathyPolyT,lat<75.502),aes(long,lat,z=id),color="black",breaks=c(-300,-250,-100,-50,-10,0))+
  scale_fill_manual(values=water.col(5),name="Depth (m)",labels=c(">250 m","250-100 m","100-50 m","50-10 m","<10 m"))+
  geom_star(data=filter(spatial,Type=="Release"),aes(Longitude,Latitude),fill="firebrick2",size=11)+ #Campsite
  geom_point(aes(x=-81.10518, y=72.35267),size=8,shape=22,fill="seagreen3")+ #Fyke Net
  new_scale_fill()+
  #geom_polygon(data=plotBuffer,aes(long,lat,group=group),fill="white",alpha=0.25,color="black",size=0.75)+
  geom_polygon(data=rangeBuffer,aes(long,lat,group=group),fill="transparent",color="black",size=1.1)+
  geom_polygon(data=filter(rangeBuffer,N>0),aes(long,lat,group=group,fill=N),color="black",size=1.25)+
  scale_fill_distiller(palette = "YlOrRd",trans="log",na.value="transparent",direction=1,
                       name="Number of  \nDetections",breaks=c(1,20,400,8000,160000),
                       labels=c("","20","400","8,000","160,000"))+
  #Gate labels
  geom_label(aes(x=-80.945,y=72.471,label="A"),color="white",fill=viridis(7,option="B",end=0.9)[1],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.025,y=72.429,label="B"),color="white",fill=viridis(7,option="B",end=0.9)[2],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.096,y=72.391,label="C"),color="white",fill=viridis(7,option="B",end=0.9)[3],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.154,y=72.364,label="D"),color="white",fill=viridis(7,option="B",end=0.9)[4],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.078,y=72.305,label="E"),color="black",fill=viridis(7,option="B",end=0.9)[5],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.119,y=72.277,label="F"),color="black",fill=viridis(7,option="B",end=0.9)[6],size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.188,y=72.246,label="G"),color="black",fill=viridis(7,option="B",end=0.9)[7],size=18,vjust=0.25,hjust=0.5)+
  #"Legend" of shape types
  annotate("rect",xmin=-80.85,xmax=-80.675,ymin=72.379,ymax=72.425,fill="white",alpha=0.85,color="black")+
  geom_text(aes(x=-80.84,y=72.42,label="Location"),hjust=0,size=10)+
  geom_text(aes(x=-80.815,y=72.406,label="Campsite"),hjust=0,size=8)+
  geom_text(aes(x=-80.815,y=72.396,label="Fyke Net"),hjust=0,size=8)+
  geom_text(aes(x=-80.815,y=72.386,label="Receivers"),hjust=0,size=8)+
  geom_star(aes(x=-80.83,y=72.405),size=8,fill="firebrick2")+ #Camp
  geom_point(aes(x=-80.83, y=72.395),size=8,shape=22,fill="seagreen3")+ #Fyke Net
  geom_point(aes(x=-80.83, y=72.385),size=8,shape=21,fill="white")+ #All the receiver points
  geom_point(aes(x=-80.83, y=72.385),size=2)+ #All the receiver points
  theme(legend.position=c(0.025,0.51),legend.justification=0,legend.title=element_text(size=30),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.88),color="black",size=1),
        legend.key.height=unit(26,"points"),legend.text=element_text(size=25),legend.margin=margin(10,20,20,15),
        axis.title=element_text(size=32),axis.text=element_text(size=25),
        axis.ticks=element_line(size=2),axis.ticks.length=unit(0.2,"cm"),
        plot.margin=margin(60,100,20,20), panel.border = element_rect(fill=NA,color=NA,size=1.25))+
  scale_y_continuous(expand=expansion(0),breaks=c(72.25,72.3,72.35,72.4,72.45,72.5),
                     labels=c("72.25°","72.30°","72.35°","72.40°","72.45°","72.50°"),name="Latitude")+
  scale_x_continuous(expand=expansion(0),breaks=c(-81.4,-81.2,-81,-80.8),
                     labels=c("-81.4°","-81.2°","-81.0°","-80.8°"),name="Longitude")+
  scalebar(location="bottomleft",anchor=c(x=-81.33,y=72.4935),dist=5,dist_unit="km",transform=T,height=0.025,
           x.min=-81.49,x.max=-80.61,y.min=72.24,y.max=72.501,st.dist=0.04,st.size=13,st.bottom=F)+
  north(x.min=-81.49,x.max=-80.61,y.min=72.24,y.max=72.501,
        anchor=c(x=-81.4475,y=72.5),location="topleft",symbol=12)+
  annotate("text",x=-81.4025,y=72.47,label="N",color="black",size=18,fontface="bold")+
  coord_cartesian(ylim=c(72.236,72.502),clip="off")
#x=-81.039,y=72.329 #Alternate D, below the line

dev.off()


####Inset####
load("~/UNH Research/Data/CSV_Files/Movement/insetMap_googleMap.R")
countries<-map_data("world")%>%filter(region%in%c("Canada","USA"))
countriesLabels<-countries%>%
  group_by(region,subregion)%>%summarise(long=mean(long),lat=mean(lat))
canada<-readOGR("~/UNH Research/Data/Mapping Data/Canada/lpr_000b16a_e.shp")
canada <- spTransform(canada, CRSobj = CRS("+proj=longlat +ellps=WGS84"))

canadaDF<-fortify(canada)
canadaDF<-canadaDF%>%
  mutate(id=factor(id),
         idName=ifelse(id=="0",canada@data$PRENAME[1],ifelse(id=="1",canada@data$PRENAME[2],ifelse(id=="2",canada@data$PRENAME[3],
                                                                                                   ifelse(id=="3",canada@data$PRENAME[4],ifelse(id=="4",canada@data$PRENAME[5],ifelse(id=="5",canada@data$PRENAME[6],
                                                                                                                                                                                      ifelse(id=="6",canada@data$PRENAME[7],ifelse(id=="7",canada@data$PRENAME[8],ifelse(id=="8",canada@data$PRENAME[9],
                                                                                                                                                                                                                                                                         ifelse(id=="9",canada@data$PRENAME[10],ifelse(id=="10",canada@data$PRENAME[11],
                                                                                                                                                                                                                                                                                                                       ifelse(id=="11",canada@data$PRENAME[12],canada@data$PRENAME[13])))))))))))))
#filter(canadaDF,id%in%c("0","1","2","3","4","5","6","12"))
canadaMap<-canadaDF%>%
  group_by(idName)%>%
  mutate(labLong=mean(long),labLat=ifelse(idName=="Ontario",mean(lat)+4,mean(lat)))%>%
  dplyr::select(labLong,labLat,idName)%>%unique()

inset<-image_graph(width=525,height=525,res=72,bg="transparent")
ggmap(insetMap)+
  geom_polygon(data=canadaDF,aes(long,lat,group=group),fill="white",color="black",size=0.5,alpha=0.2)+
  geom_text(data=canadaMap%>%filter(idName%in%c("Nunavut","Manitoba","Ontario","Quebec")),
            aes(labLong,labLat,label=idName),color="white")+
  #scale_x_continuous(limits=c(-115,-45))+
  #scale_y_continuous(limits=c(40,74))+
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
        plot.background = element_blank(), panel.border=element_rect(fill=NA,color="black",size=2))+
  scalebar(location="bottomleft",dist=1000,dist_unit="km",transform=T,st.color="white",st.bottom=F,st.dist=0.05,
           x.min=-115,x.max=-45,y.min=40,y.max=74)+
  #geom_point(aes(-70.934167,43.135),size=4,color="red") #Added to get UNH in right place
  geom_rect(aes(xmin=-81.89,xmax=-79.6,ymin=71.7,ymax=72.7),fill="transparent",color="firebrick1",size=2)+ #Added to surround Tremblay
  coord_map(xlim=c(-118,-45),ylim=c(38,74))

dev.off()

fullMap<-image_composite(bathyMap,inset,offset="+800+710") #Fit outside the map
fullMap





# Multivariate stats ------------------------------------------------------
# Mirroring methods from Lowe et al. 2020 (Animal Biotelemetry) for this procedure

detections<-ungroup(V9detections)
#Spatial df

#Need to be in the same wd as spatial.csv
receiver_distMat<-read.csv("distances.csv")
rownames(receiver_distMat)<-receiver_distMat$X
receiver_distMat<-receiver_distMat[,-1]

#### Creating the sequences ####
dayDetections<-dplyr::select(detections,Date.and.Time,Fish_ID,Station.name,Array,Longitude,Latitude,Receiver_Type)%>%
  mutate(Receiver_Frequency=ifelse(grepl("180",Receiver_Type),180,69),
         Date=floor_date(Date.and.Time,unit="day"),
         timeIn=as.numeric(difftime(Date.and.Time,Date,units="hours")))
dayDetections2<-dayDetections%>% 
  dplyr::group_by(Date,Fish_ID,Array)%>%
  dplyr::summarise(stateTime=max(timeIn)) #difference between last time in state and what it rounded down to
#Important to set when it "left" that state within the time step


detailedDays<-left_join(dayDetections2,tags)%>%
  dplyr::select(Date,Fish_ID,Array,stateTime,Release_date,Species,Tag_Type,Year)%>%
  mutate(days_since_722=ifelse(Year=="2017", difftime(Date,as.Date("2017-07-22",tz="America/Iqaluit"),units="days"),
                               ifelse(Year=="2018",difftime(Date,as.Date("2018-07-22",tz="America/Iqaluit"),units="days"),
                                      difftime(Date,as.Date("2019-07-22",tz="America/Iqaluit"),units="days"))))%>% #Again, change to the time steps (but from the earliest release time)
  ungroup()%>%
  dplyr::group_by(Fish_ID)%>%
  mutate(start_date=min(days_since_722),
         end_date=max(days_since_722))



# Defining seasons plotting sequence ---------------------------------------

#Defining the seasons of each detection--using the raw detection times and locations
rawICE<-dayDetections%>%
  mutate(ice_season=ifelse(Date.and.Time>="2017-07-03"&Date.and.Time<"2017-07-30","Ice Melt", #Ice melt of 2017
                           ifelse(Date.and.Time>="2017-07-30"&Date.and.Time<="2017-10-09","Ice Free", #Ice free of 2017
                                  ifelse(Date.and.Time>"2017-10-09"&Date.and.Time<="2017-10-17","Ice Freeze", #Ice freeze 2017
                                         ifelse(Date.and.Time>="2018-07-04"&Date.and.Time<"2018-08-13","Ice Melt",
                                                ifelse(Date.and.Time>="2018-08-13"&Date.and.Time<="2018-09-22","Ice Free",
                                                       ifelse(Date.and.Time>"2018-09-22"&Date.and.Time<="2018-10-22","Ice Freeze",
                                                              ifelse(Date.and.Time>="2019-06-28"&Date.and.Time<"2019-07-23","Ice Melt",
                                                                     ifelse(Date.and.Time>="2019-07-23"&Date.and.Time<="2019-10-19","Ice Free",
                                                                            ifelse(Date.and.Time>"2019-10-19"&Date.and.Time<="2019-11-09","Ice Freeze",
                                                                                   ifelse(Date.and.Time>"2019-11-09"&Date.and.Time<="2020-07-02","Frozen",
                                                                                          ifelse(Date.and.Time>"2020-07-02"&Date.and.Time<="2020-07-24","Ice Melt",
                                                                                                 ifelse(Date.and.Time>"2020-07-24","Ice Free","Frozen")))))))))))),
         ice_season_length=ifelse(Date.and.Time>="2017-07-03"&Date.and.Time<"2017-07-30",
                                  as.numeric(difftime("2017-07-30","2017-07-03")),
                                  ifelse(Date.and.Time>="2017-07-30"&Date.and.Time<="2017-10-09",
                                         as.numeric(difftime("2017-10-09","2017-07-30")),
                                         ifelse(Date.and.Time>"2017-10-09"&Date.and.Time<="2017-10-17",
                                                as.numeric(difftime("2017-10-17","2017-10-09")),
                                                ifelse(Date.and.Time>"2017-10-17"&Date.and.Time<"2018-07-04",
                                                       as.numeric(difftime("2018-07-04","2017-10-17")),
                                                       ifelse(Date.and.Time>="2018-07-04"&Date.and.Time<"2018-08-13",
                                                              as.numeric(difftime("2018-08-13","2018-07-04")),
                                                              ifelse(Date.and.Time>="2018-08-13"&Date.and.Time<="2018-09-22",
                                                                     as.numeric(difftime("2018-09-22","2018-08-13")),
                                                                     ifelse(Date.and.Time>"2018-09-22"&Date.and.Time<="2018-10-22",
                                                                            as.numeric(difftime("2018-10-22","2018-09-22")),
                                                                            ifelse(Date.and.Time>"2018-10-22"&Date.and.Time<"2019-06-28",
                                                                                   as.numeric(difftime("2019-06-28","2018-10-22")),
                                                                                   ifelse(Date.and.Time>="2019-06-28"&Date.and.Time<"2019-07-23",
                                                                                          as.numeric(difftime("2019-07-23","2019-06-28")),
                                                                                          ifelse(Date.and.Time>="2019-07-23"&Date.and.Time<="2019-10-19",
                                                                                                 as.numeric(difftime("2019-10-19","2019-07-23")),
                                                                                                 ifelse(Date.and.Time>"2019-10-19"&Date.and.Time<="2019-11-09",
                                                                                                        as.numeric(difftime("2019-11-09","2019-10-19")),
                                                                                                        ifelse(Date.and.Time>"2019-11-09"&Date.and.Time<="2020-07-02",
                                                                                                               as.numeric(difftime("2020-07-02","2019-11-09")),
                                                                                                               ifelse(Date.and.Time>"2020-07-02"&Date.and.Time<="2020-07-23",
                                                                                                                      as.numeric(difftime("2020-07-23","2020-07-02")),
                                                                                                                      ifelse(Date.and.Time>"2020-07-23",
                                                                                                                             floor(as.numeric(difftime(max(Date.and.Time),"2020-07-23"))),NA)))))))))))))),
         season=ifelse(ice_season=="Frozen","Winter","Summer"),
         season_year=ifelse(season=="Summer",year(Date.and.Time),
                            ifelse(yday(Date.and.Time)>214,year(Date.and.Time),year(Date.and.Time)-1)))%>%
  left_join(tags)%>%
  mutate(season_length=ifelse(season_year==2017 & Year==2017 & season=="Summer", #Summer 2017 (2017 tags)
                              as.numeric(difftime("2017-10-17",Release_date)),
                              ifelse(Date.and.Time>"2017-10-17"&Date.and.Time<"2018-07-04" & Tag_Type=="V9", #Winter 2017 (2017 V9s)
                                     as.numeric(difftime("2018-07-04","2017-10-17")),
                                     ifelse(Date.and.Time>"2017-10-17"&Date.and.Time<"2018-07-04" & Tag_Type=="V5", #Winter 2017 (2017 V5s)
                                            as.numeric(difftime(Release_date+150,"2017-10-17")),
                                            ifelse(season_year==2018 & Year==2017, #Summer 2018 (2017 tags)
                                                   as.numeric(difftime(Release_date+405,"2018-07-04")),
                                                   ifelse(season_year==2018 & Year==2018 & season=="Summer", #Summer 2018 (2018 tags)
                                                          as.numeric(difftime("2018-10-22",Release_date)),
                                                          ifelse(Date.and.Time>"2018-10-22"&Date.and.Time<"2019-06-28", #Winter 2018 (2018 tags)
                                                                 as.numeric(difftime("2019-06-28","2018-10-22")),
                                                                 ifelse(season_year==2019 & Year==2018, #Summer 2019 (2018 tags)
                                                                        as.numeric(difftime(Release_date+405,"2019-06-28")),
                                                                        ifelse(season_year==2019 & Year==2019 & season=="Summer", #Summer 2019 (2019 tags)
                                                                               as.numeric(difftime("2019-11-09",Release_date)),
                                                                               ifelse(Date.and.Time>"2019-11-09"&Date.and.Time<"2020-07-02" & Tag_Type=="V9", #Winter 2019 (2019 V9s)
                                                                                      as.numeric(difftime("2020-07-02","2019-11-09")),
                                                                                      ifelse(Date.and.Time>"2019-11-09"&Date.and.Time<"2020-07-02" & Tag_Type%in%c("V9A","V9P"), #Winter 2019 (2019 V9APs)
                                                                                             as.numeric(difftime(Release_date+255,"2019-11-09")),
                                                                                             ifelse(Date.and.Time>="2020-07-02", #Summer 2020 (2019 tags)
                                                                                                    as.numeric(difftime(Release_date+405,"2020-07-02")),NA))))))))))))



rawICE%>%
  dplyr::select(season_year,season,ice_season_length)%>%distinct()%>%
  group_by(season_year,season)%>%
  summarise(season_length=sum(ice_season_length))




#### Impute across sequences ####



seqDays<-as.data.frame(cbind(0:max(detailedDays$days_since_722),sort(rep(unique(detailedDays$Fish_ID),434))))
colnames(seqDays)<-c("days_since_722","Fish_ID")
seqDays$days_since_722<-as.numeric(seqDays$days_since_722)

test<-full_join(detailedDays,seqDays)%>%
  arrange(Fish_ID,days_since_722)%>%
  dplyr::group_by(Fish_ID)%>%
  dplyr::mutate(Array=ifelse(days_since_722==0 & is.na(Array),"D",Array))%>%
  fill(Array,Release_date:end_date,.direction="downup")%>%
  mutate(Date=as.Date(paste0(Year,"-07-22"))+days_since_722)

#When is a fish at more than one line in the same day
testDups<-test[duplicated(test[,c("Fish_ID","days_since_722")])|
                 duplicated(test[,c("Fish_ID","days_since_722")],fromLast=T),]
#Make a rule to only keep the time from the majority of that time (where it was detected over 12 hours at)
testKeep<-testDups%>%
  dplyr::group_by(Fish_ID,days_since_722)%>%
  mutate(keep=ifelse(min(stateTime)>12,min(stateTime),max(stateTime)))%>% #This is the rule for a day, but would need to adjust to other time steps to whatever the majority is (like 30 min for an hour time step)
  subset(keep==stateTime)%>%
  dplyr::select(-keep)%>%
  unique()
#All the non-duplicated ones
testOrig<-test[(!duplicated(test[,c("Fish_ID","days_since_722")]) &
                  !duplicated(test[,c("Fish_ID","days_since_722")],fromLast=T)),]
#Put the filter duplicates and non-duplicates back together
filledDetections<-bind_rows(testKeep,testOrig)

#Percent Days with double detections
nrow(testKeep)/nrow(filledDetections)*100

#Imputation, following Lowe et al. style reporting
imputation<-filledDetections%>%
  filter(days_since_722>=start_date & days_since_722<=end_date)%>%
  dplyr::group_by(Fish_ID)%>%
  dplyr::summarise(impute_num=sum(is.na(stateTime)),
                   impute_tot=n(),
                   impute=impute_num/impute_tot)
#MEAN
mean(imputation$impute) #Right in the range, Lowe et al. had 50-63%
#SEM
sd(imputation$impute)/sqrt(nrow(imputation))


filledDetections<-filledDetections%>%
  mutate(ice_season=case_when(Date>="2017-07-03"&Date<"2017-07-30"~"Ice Melt", #Ice melt of 2017
                              Date>="2017-07-30"&Date<="2017-10-09"~"Ice Free", #Ice free of 2017
                              Date>"2017-10-09"&Date<="2017-10-17"~"Ice Freeze", #Ice freeze 2017
                              Date>="2018-07-04"&Date<"2018-08-13"~"Ice Melt",
                              Date>="2018-08-13"&Date<="2018-09-22"~"Ice Free",
                              Date>"2018-09-22"&Date<="2018-10-22"~"Ice Freeze",
                              Date>="2019-06-28"&Date<"2019-07-23"~"Ice Melt",
                              Date>="2019-07-23"&Date<="2019-10-19"~"Ice Free",
                              Date>"2019-10-19"&Date<="2019-11-09"~"Ice Freeze",
                              Date>"2019-11-09"&Date<="2020-07-02"~"Frozen",
                              Date>"2020-07-02"&Date<="2020-07-24"~"Ice Melt",
                              Date>"2020-07-24"~"Ice Free",
                              TRUE~"Frozen"),
         season=ifelse(ice_season=="Frozen","Winter","Summer"),
         Array=gsub("0","",Array))


#### Dissimilarity ####

#Right way
locationMat<-filledDetections%>%
  ungroup()%>%
  filter(days_since_722>=start_date & days_since_722<=end_date)%>%
  dplyr::select(Fish_ID,Tag_Type,Date=days_since_722,Array)%>%
  arrange(Date,Fish_ID)%>%
  mutate(numLine=match(Array,LETTERS),
         numLine=ifelse(is.na(numLine),0,numLine))%>% #Convert to numbers to do an operation on them for distance of change
  group_by(Fish_ID)%>%
  mutate(dist=abs(lead(numLine)-numLine)) #Just to look at the movement order

#Percent of movements that are 2+ order movements after imputing 
nrow(subset(locationMat,dist>1))/nrow(subset(locationMat,!is.na(dist)))*100 #0.03% of movements (5/15064)

#Create the sequence matrix of all the different individuals locations at all times
sequence.Mat<-pivot_wider(locationMat,id_cols="Fish_ID",values_from="Array",names_from="Date",
                          values_fn=max) #Rotate to rows of sequences, fn because the AP have an A and P, but they're the same location so it's just keeping 1 of the 2 and it doesn't matter which
sequence.Mat<-arrange(sequence.Mat,Fish_ID)
sequence.Mat2<-as.matrix(sequence.Mat[,2:ncol(sequence.Mat)]) #Crop to just the sequences
rownames(sequence.Mat2)<-sequence.Mat$Fish_ID

#Make all missing data the same--give it a state of it's own
#Can also accomplish the same thing with "right=NA" in the seqdef function
#If I do it with the right=NA then I have to remember to do that down below with rep seqs too
sequence.Mat2[is.na(sequence.Mat2)]<-"NA"


#Making it a state sequence object
seqForm<-seqdef(sequence.Mat2)
rownames(seqForm)<-sequence.Mat$Fish_ID

customCosts<-matrix(data=c(0,1,2,2,2,2,2,1,
                           1,0,1,2,2,2,2,1,
                           2,1,0,1,2,2,2,1,
                           2,2,1,0,1,2,2,1,
                           2,2,2,1,0,1,2,1,
                           2,2,2,2,1,0,1,1,
                           2,2,2,2,2,1,0,1,
                           1,1,1,1,1,1,1,0),nrow=8,byrow=T)

rownames(seqForm)<-sequence.Mat$Fish_ID
seqDist<-seqdist(seqForm,method="OM",sm=customCosts,indel=0.95,with.missing=F)
max(seqDist[row(seqDist)!=col(seqDist)])
#Converting to a matrix
plotDist<-seqDist
rownames(plotDist)<-sequence.Mat$Fish_ID
colnames(plotDist)<-sequence.Mat$Fish_ID
for (i in 1:59) {
  for (j in 1:59) {
    plotDist[i,j]<-ifelse(j>i,NA,plotDist[i,j])
  }
}
plotDist_gg<-as.data.frame(plotDist)%>%
  mutate(trans1=rownames(plotDist))%>%
  pivot_longer(cols=c("1613":"6458"),names_to = "trans2", values_to = "Distance")

#Plotting that matrix
ggplot(plotDist_gg)+
  geom_tile(aes(fct_reorder(trans1,trans1,.desc=T),trans2,fill=Distance))+
  xlab("Transmitter")+
  scale_y_discrete(expand=expansion(0),name="Transmitter")+
  scale_fill_viridis_c(option="A",na.value = "transparent",breaks=c(0,100,200,300,400,500,600,700))+
  theme_classic(base_size=25)+
  theme(panel.border = element_blank(),panel.grid=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,size=20),
        legend.position=c(0.85,0.7),legend.key.height = unit(60,"pt"),legend.key.width = unit(20,"pt"),
        legend.title=element_text(size=30),axis.title=element_text(size=31),legend.text=element_text(size=25))

#### Clustering ####

set.seed(42)
par(mar=c(7,5,5,5))

rownames(seqDist)<-sequence.Mat$Fish_ID #Need the names in the matrix
colnames(seqDist)<-sequence.Mat$Fish_ID
seqClust<-pvclust(seqDist,method.hclust="ward.D2")
#Plotting the clusters
plot(seqClust,print.pv=c("au"),col.pv=c(au=2))
pvrect(seqClust,alpha=0.94,type="gt",border="dodgerblue3",lwd=2,col=rgb(0,0.4,0.9,alpha=0.1))

#Assigning out the groups
bootGroups<-pvpick(seqClust,alpha=0.94,type="gt")
trans<-as.character(sequence.Mat$Fish_ID)
groups<-numeric(length(trans))
for (t in 1:length(trans)) {
  for (i in 1:length(bootGroups$clusters)) { #Need to make the set equal to the number of groups identified from bootstrapping
    groups[t]=ifelse(sum(bootGroups$clusters[[i]]==trans[t])==1,i,groups[t])
  }
}
groups<-data.frame(Fish_ID=as.character(sequence.Mat$Fish_ID),Group=groups)%>%
  left_join(dplyr::select(tags,Fish_ID,Year))



#### Representative Sequences ####
par(mfrow=c(1,1),mar=c(5,5,2,2))


#Start with everybody
seqdplot(seqForm,group=numeric(59),with.legend="none",main="All Tags",cpal(seqForm)<-c(pnew,"grey50"),border=NA)

##### Group 1 #####
seq1<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[1]]) #Crop the sequence matrix down to the first group
seq1<-seqdef(seq1) #make it a sequence object again
cost1<-matrix(data=c(0,1,2,2,1,
                     1,0,1,2,1,
                     2,1,0,1,1,
                     2,2,1,0,1,
                     1,1,1,1,0),nrow=5,byrow=T)
plot(seqrep(seq1,criterion="dist",coverage=0.5,
            diss=seqdist(seq1,method="OM",sm=cost1,indel=0.95)),dmax=110,
     cpal(seq1)<-c(viridis(7,option="B",end=0.9)[4:7],"grey50"),main=list("Group 1, n=2",cex=2),border=NA,ylab=NA)


##### Group 2 #####
seq2<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[2]]) #Crop the sequence matrix down to the first group
seq2<-seqdef(seq2) #make it a sequence object again
cost2<-matrix(data=c(0,1,2,1,
                     1,0,1,1,
                     2,1,0,1,
                     1,1,1,0),nrow=4,byrow=T)


plot(seqrep(seq2,criterion="dist",coverage=0.5,
            diss=seqdist(seq2,method="OM",sm=cost2,indel=0.95)),dmax=110,
     cpal(seq2)<-c(viridis(7,option="B",end=0.9)[3:5],"grey50"),main=list("Group 2, n=2",cex=2),border=NA,ylab=NA)

##### Group 3 #####
seq3<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[3]]) #Crop the sequence matrix down to the first group
seq3<-seqdef(seq3) #make it a sequence object again
cost3<-customCosts

plot(seqrep(seq3,criterion="dist",coverage=0.5,
            diss=seqdist(seq3,method="OM",sm=cost3,indel=0.95)),dmax=110,
     cpal(seq3)<-c(viridis(7,option="B",end=0.9),"grey50"),main=list("Group 3, n=18",cex=2),border=NA,ylab=NA)

##### Group 4 #####
seq4<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[4]]) #Crop the sequence matrix down to the first group
seq4<-seqdef(seq4) #make it a sequence object again
cost4<-matrix(data=c(0,1,2,1,
                     1,0,1,1,
                     1,1,0,1,
                     1,1,1,0),nrow=4,byrow=T)

plot(seqrep(seq4,criterion="dist",coverage=0.5,
            diss=seqdist(seq4,method="OM",sm=cost4,indel=0.95)),dmax=110,
     cpal(seq4)<-c(viridis(7,option="B",end=0.9)[4:6],"grey50"),main=list("Group 4, n=2",cex=2),border=NA,ylab=NA)

##### Group 5 #####
seq5<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[5]]) #Crop the sequence matrix down to the first group
seq5<-seqdef(seq5) #make it a sequence object again
cost5<-matrix(data=c(0,1,2,2,2,2,1,
                     1,0,1,2,2,2,1,
                     2,1,0,1,2,2,1,
                     2,2,1,0,1,2,1,
                     2,2,2,1,0,1,1,
                     2,2,2,2,1,0,1,
                     1,1,1,1,1,1,0),nrow=7,byrow=T)


plot(seqrep(seq5,criterion="dist",coverage=0.5,
            diss=seqdist(seq5,method="OM",sm=cost5,indel=0.95)),dmax=110,
     cpal(seq5)<-c(viridis(7,option="B",end=0.9)[2:7],"grey50"),main=list("Group 5, n=4",cex=2),border=NA,ylab=NA)

##### Group 6 #####
seq6<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[6]]) #Crop the sequence matrix down to the first group
seq6<-seqdef(seq6) #make it a sequence object again
cost6<-matrix(data=c(0,1,1,
                     1,0,1,
                     1,1,0),nrow=3,byrow=T)


plot(seqrep(seq6,criterion="dist",coverage=0.5,
            diss=seqdist(seq6,method="OM",sm=cost6,indel=0.95)),dmax=110,
     cpal(seq6)<-c(viridis(7,option="B",end=0.9)[4:5],"grey50"),main=list("Group 6, n=3",cex=2),border=NA, ylab=NA)


##### Group 7 #####
seq7<-subset(sequence.Mat2,rownames(sequence.Mat2)%in%bootGroups$clusters[[7]]) #Crop the sequence matrix down to the first group
seq7<-seqdef(seq7) #make it a sequence object again
cost7<-matrix(data=c(0,1,2,2,2,2,1,
                     1,0,1,2,2,2,1,
                     2,1,0,1,2,2,1,
                     2,2,1,0,1,2,1,
                     2,2,2,1,0,1,1,
                     2,2,2,2,1,0,1,
                     1,1,1,1,1,1,0),nrow=7,byrow=T)

plot(seqrep(seq7,criterion="dist",coverage=0.5,
            diss=seqdist(seq7,method="OM",sm=cost7,indel=0.95)),dmax=110,
     cpal(seq7)<-c(viridis(7,option="B",end=0.9)[1:6],"grey50"),main=list("Group 7, n=28",cex=2),border=NA,ylab=NA)





#### Plotting all sequences in one figure (with group IDs) ####


realtransitions<-filledDetections%>%
  arrange(Fish_ID,days_since_722)%>%
  filter(stateTime>0)%>% #Lets me crop them to the battery life
  dplyr::group_by(Fish_ID)%>%
  filter(days_since_722==min(days_since_722)| #Need the 0
           days_since_722==max(days_since_722)| #And the end to stop the last bar
           Array!=lag(Array)| #Any starts of a transition
           Array!=lead(Array)) #And ends of a transition



plotT<-realtransitions%>%
  dplyr::group_by(Fish_ID)%>%
  mutate(start=ifelse(lag(Array)!=Array | days_since_722==start_date,days_since_722,NA),
         end=ifelse(lead(Array)==Array,lead(days_since_722),NA),
         end=ifelse(days_since_722==end_date&lag(Array)!=Array,days_since_722,end),
         end=ifelse(start>0 & is.na(end),start,end))
plotT<-unique(plotT[!is.na(plotT$start),c("Date","Fish_ID","Year","Array","start","end")])


#For dates on the geom_rects:
#Using 7/22 as 0, so the day is DATE-7/22, so any dates before that in the year are -Inf
#7/22 in the next year is 365 so 365-(7/22-DATE) is the latter summer
#First rect is the start of ice-off and end of ice-on, so the whole "summer" period
#Inner rect, second, is the end and start of ice-off and ice-on, respectively, so is the ice-free season
#Best plot
ggplot(plotT)+
  geom_rect(data=data.frame(Year=2017), 
            aes(xmin=-Inf,xmax=87,ymin=-Inf,ymax=Inf), ############# Summer of 2017 end (10-17)
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_rect(data=data.frame(Year=2017),
            aes(xmin=347,xmax=Inf,ymin=-Inf,ymax=Inf), ############# Summer of 2018 (for the 2018s) start (7-4) 
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_rect(data=data.frame(Year=2018),
            aes(xmin=-Inf,xmax=92,ymin=-Inf,ymax=Inf), ############# Summer of 2018 (for the 2019s) end (10-22)
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_rect(data=data.frame(Year=2018),
            aes(xmin=341,xmax=Inf,ymin=-Inf,ymax=Inf), ############# Summer of 2019 (for the 2018s) start (6-28) 
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_rect(data=data.frame(Year=2019),
            aes(xmin=-Inf,xmax=110,ymin=-Inf,ymax=Inf), ############ Summer of 2019 (for the 2019s) end (11-9)
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_rect(data=data.frame(Year=2019),
            aes(xmin=345,xmax=Inf,ymin=-Inf,ymax=Inf), ############# Summer of 2020 start (7-02) 
            fill=rgb(0.5,0.5,0.5,alpha=0.5))+
  geom_segment(aes(y=Fish_ID,yend=Fish_ID,
                   x=start,xend=end+1,color=Array),lwd=7)+
  geom_label(data=groups,aes(y=Fish_ID,x=440,label=Group),hjust=0.5,vjust=0.5,label.size=NA,size=6)+
  facet_grid(Year~.,scales="free_y",space="free")+
  scale_color_viridis_d(option="B",end=0.9,name="Gate")+
  scale_y_discrete(expand=c(0,0.75))+
  scale_x_continuous(expand=c(0,5),breaks=c(0,100,200,300,400),
                     label=c("7-22","10-30","2-7","5-18","8-26"),name="Date")+
  theme(axis.text.y=element_text(size=19),axis.text.x=element_text(size=25),axis.title=element_text(size=30),
        legend.margin=margin(0,0,0,-10),legend.title=element_text(size=30),
        legend.text=element_text(size=25),panel.grid.minor=element_blank(),
        legend.key.size=unit(1.15,"cm"),strip.text.y=element_text(angle=0,size=27.5),
        plot.title.position = "plot",plot.title = element_text(hjust=0.86,size=18,vjust=0))+
  ggtitle(label="Group")+
  ylab("Transmitter ID")




# Inter-gate Movements (by season) ---------------------------------------------------

# Number of total movements normal to length of season
seasonMoves<-rawICE%>%
  filter(Tag_Type!="V5")%>%
  arrange(Date.and.Time)%>%
  group_by(Fish_ID,season,season_year,season_length)%>%
  dplyr::summarise(Move=ifelse(Array==lag(Array),0,1),
                   nMoves=sum(Move,na.rm=T),
                   moves_norm=nMoves/season_length)%>%dplyr::select(-Move)%>%unique()

seasonMoves%>%group_by(season)%>%summarise(min=min(nMoves),max=max(nMoves),mean=mean(nMoves),sem=sd(nMoves)/sqrt(n()))


#Or maybe a paired t-test makes more sense
pairedMove<-pivot_wider(seasonMoves,id_cols="Fish_ID",names_from=c("season","season_year"),values_from="moves_norm",values_fn=mean)
t.test(Pair(Summer_2017,Winter_2017)~1,data=pairedMove) #2017 p=0.005715,
t.test(Pair(Summer_2018,Winter_2018)~1,data=pairedMove) #2018 p=0.008149,
t.test(Pair(Summer_2019,Winter_2019)~1,data=pairedMove) #2019 p=0.4851

#N
N<-rawICE%>%
  filter(Tag_Type!="V5")%>%
  group_by(season_year,season)%>%summarise(N=n_distinct(Fish_ID))%>%
  mutate(season=factor(season,levels=c("Summer","Winter")))

#Plot
ggplot(seasonMoves)+
  geom_boxplot(aes(season,moves_norm,fill=season),width=0.9,outlier.shape=NA,show.legend=F)+
  geom_path(aes(season,moves_norm,group=Fish_ID),alpha=0.5,
            position=position_jitter(width=0.25,seed=42))+
  geom_point(aes(season,moves_norm,fill=season),shape=21,size=4.1,show.legend=F,
             position=position_jitter(width=0.25,seed=42))+
  geom_label(data=N,aes(season,0.081,label=paste("N =",N)),size=8)+
  geom_segment(data=data.frame(season_year=c("2017","2018")),
               aes(x="Summer",xend="Winter",y=0.0865,yend=0.0865),lineend="round",size=2)+
  geom_text(data=data.frame(season_year=c("2017","2018")),aes(1.5,0.0875,label="***"),size=8)+
  scale_fill_viridis_d(begin=0.2,direction=-1)+
  scale_x_discrete(expand=expansion(mult=0.5),name="Season",labels=c("Ice-Free","Ice-Covered"))+
  scale_y_continuous(name="Inter-gate movements/day")+
  theme(panel.grid.major.x=element_blank())+
  facet_wrap(~season_year,nrow=1)




# Roaming Index  ----------------------------------------------




roamingICE<-rawICE%>%
  filter(Tag_Type!="V5")%>%
  group_by(Fish_ID,season_year,season,season_length)%>%
  summarise(locations=n_distinct(Station.name),
            roam=locations/37)%>%
  ungroup()%>%
  mutate(roam_norm=roam/season_length,
         arcsine_roam=asin(sqrt(roam_norm)),
         season=factor(season,levels=c("Summer","Winter")))

seasonroaming<-roamingICE%>%group_by(season)%>%summarise(mean=mean(roam_norm),sem=sd(roam_norm)/sqrt(n()))
seasonroaming$mean[seasonroaming$season=="Summer"]/seasonroaming$mean[seasonroaming$season=="Winter"]

#Or maybe a paired t-test makes more sense
pairedRoaming<-pivot_wider(roamingICE,id_cols="Fish_ID",names_from=c("season","season_year"),values_from="roam_norm",values_fn=mean)
t.test(Pair(Summer_2017,Winter_2017)~1,data=pairedRoaming) #2017 p<0.001, 
t.test(Pair(Summer_2018,Winter_2018)~1,data=pairedRoaming) #2018 p<0.001,
t.test(Pair(Summer_2019,Winter_2019)~1,data=pairedRoaming) #2019 p=0.038


#N
N<-rawICE%>%
  filter(Tag_Type!="V5")%>%
  group_by(season_year,season)%>%summarise(N=n_distinct(Fish_ID))%>%
  mutate(season=factor(season,levels=c("Summer","Winter")))

#Plotting
ggplot(roamingICE,aes(season,roam_norm,size=Burden))+
  geom_boxplot(aes(season,roam_norm,fill=season),color="black",lwd=0.9,width=0.9,outlier.shape=NA,show.legend = F)+
  geom_path(aes(group=Fish_ID),color="black",size=0.5,alpha=0.25,
            position=position_jitter(width=0.25,seed=42))+
  geom_point(aes(fill=season),shape=21,show.legend = F,
             position=position_jitter(width=0.25,seed=42))+
  geom_label(data=N,aes(season,y=0.0051,label=paste("N =",N)),size=8)+
  geom_segment(data=data.frame(season_year=c("2017","2018","2019")),
               aes(x="Summer",xend="Winter",y=0.0053,yend=0.0053),lwd=1.5,lineend="round")+
  geom_text(data=data.frame(season_year=c("2017","2018")),aes(x=1.5,y=0.00535,label="***"),size=8)+
  geom_text(data=data.frame(season_year="2019"),aes(x=1.5,y=0.00535,label="*"),size=8)+
  scale_fill_viridis_d(begin=0.2,direction=-1)+
  scale_y_continuous(limits=c(0,0.0056),expand=c(0,0),name="Normalized Roaming Index")+
  scale_x_discrete(expand=expansion(mult=0.5),name="Season",labels=c("Ice-Free","Ice-Covered"))+
  theme(panel.grid.major.x=element_blank())+
  facet_wrap(~season_year,nrow=1)



# Residency ---------------------------------------------------------------

allSeasons<-data.frame(unique(dplyr::select(rawICE,season_year,season)))%>%
  arrange(season_year,season)
allSeasons$start_date<-ymd(c("2017-07-22","2017-10-17","2018-07-04","2018-10-22","2019-06-28","2019-11-09","2020-07-02"))
allSeasons$end_date<-ymd(c("2017-10-17","2018-07-04","2018-10-22","2019-06-28","2019-11-09","2020-07-02","2020-10-10"))



#Residency length by the two seasons and normalized to the length of the season
seasonResid<-left_join(rawICE,allSeasons)%>%
  filter(Tag_Type!="V5")%>%
  dplyr::select(Fish_ID,Date.and.Time,season,season_year,season_length,Array)%>%
  arrange(Fish_ID,Date.and.Time)%>%unique()%>%
  group_by(Fish_ID,season,season_year,season_length)%>%
  mutate(residencyEnd=ifelse(Array!=lag(Array) | 
                               Date.and.Time==min(Date.and.Time) | 
                               as.numeric(difftime(Date.and.Time,lag(Date.and.Time),units="days"))>30,"start",
                             ifelse(Array!=lead(Array) |
                                      Date.and.Time==max(Date.and.Time) |
                                      as.numeric(difftime(lead(Date.and.Time),Date.and.Time,units="days"))>30,"end","no")),
         residency=1)%>%
  filter(residencyEnd!="no")

n=0
for (i in 2:nrow(seasonResid)) {
  n<-ifelse(seasonResid$residencyEnd[i]=="start",n+1,n) #Add to the residency number
  seasonResid$residency[i]<-seasonResid$residency[i]+n
}

seasonResid<-seasonResid%>%
  group_by(Fish_ID,season,season_year,season_length,Array,residency)%>%
  summarise(start=min(Date.and.Time),
            end=max(Date.and.Time),
            resid=ceiling(as.numeric(difftime(end,start,units="days"))))%>%
  ungroup()%>%group_by(Fish_ID,season,season_year)%>%
  mutate(propResid=resid/season_length,
         totalResid=sum(resid),
         totalPropResid=totalResid/season_length)%>%
  ungroup()%>%group_by(Fish_ID,season,season_year,Array)%>%
  mutate(totalResid.Array=sum(resid),
         totalPropResid.Array=totalResid.Array/season_length)

seasonResid%>%group_by(season)%>%summarise(mean=mean(propResid),sem=sd(propResid)/sqrt(n()))
seasonresid_total<-seasonResid%>%group_by(season)%>%summarise(mean=mean(resid),sem=sd(resid)/sqrt(n()))
seasonresid_total$mean[seasonresid_total$season=="Winter"]/seasonresid_total$mean[seasonresid_total$season=="Summer"]

#A paired t-test
pairedResid<-pivot_wider(seasonResid,id_cols="Fish_ID",names_from=c("season","season_year"),values_from="propResid",values_fn=max)
t.test(Pair(Summer_2017,Winter_2017)~1,data=pairedResid) #2017 p=0.7773,
t.test(Pair(Summer_2018,Winter_2018)~1,data=pairedResid) #2018 p=0.2325,
t.test(Pair(Summer_2019,Winter_2019)~1,data=pairedResid) #2019 p=0.7682

pseasonResid<-seasonResid%>%
  group_by(season,season_year,Fish_ID)%>%
  summarise(max=max(propResid),
            sd=sd(propResid))

#N
N<-seasonResid%>%
  group_by(season_year,season)%>%
  summarise(N=n_distinct(Fish_ID))

#Plot
ggplot(pseasonResid)+
  geom_boxplot(aes(season,max,fill=season),width=0.9,outlier.shape=NA,show.legend=F)+
  geom_path(aes(season,max,group=Fish_ID),alpha=0.5,
            position=position_jitter(width=0.25,seed=42))+
  geom_point(aes(season,max,fill=season),shape=21,size=4.1,show.legend=F,
             position=position_jitter(width=0.25,seed=42))+
  geom_label(data=N,aes(season,1.05,label=paste("N =",N)),size=7)+
  scale_fill_viridis_d(begin=0.2,direction=-1)+
  scale_x_discrete(expand=expansion(mult=0.5),name="Season",labels=c("Ice-Free","Ice-Covered"))+
  scale_y_continuous(name="Max Proportion Residency Length",expand=expansion(add=0.01,0.035),
                     breaks=c(0,0.25,0.5,0.75,1))+
  theme(panel.grid.major.x=element_blank())+
  facet_wrap(~season_year,nrow=1)




# Home Range Estimation and Modeling --------------------------------------------------------------


#Getting the right data formatting for ACTEL
bio<-read_csv("tagDetails.csv", 
              col_types = cols(Release_date = col_datetime(format = "%Y-%m-%d"), 
                               Tag_ID = col_character(), Year = col_character()))%>%
  mutate(Fish_ID=ifelse(Tag_Type=="V9A",as.character(as.numeric(Tag_ID)+0.5),
                        ifelse(Tag_Type=="V9P",as.character(as.numeric(Tag_ID)-0.5),Tag_ID)))%>%
  filter(Tag_Type=="V9")
#Getting the right data they want
bio<-dplyr::select(bio,Signal=Tag_ID,Release.date=Release_date,Species,Tag_Type,Year)%>%
  mutate(Release.date=ymd_hms(paste(Release.date,"07:00:00")),
         Signal=ifelse(Tag_Type=="V9P",as.integer(Signal-1),as.integer(Signal)),
         Tag_Type=ifelse(Tag_Type%in%c("V9A","V9P"),"V9AP",Tag_Type))%>%
  filter(!duplicated(Signal))%>%
  mutate(Release.site="Camp")


detect <-read_csv("all3downloads_sculpDetections+.csv", 
                  col_types = cols(Date.and.Time = col_datetime(format = "%Y-%m-%d %H:%M:%S"), 
                                   Fish_ID = col_character(), Receiver_ID = col_character(), 
                                   Release_date = col_date(format = "%Y-%m-%d"), 
                                   Sensor.Measure = col_double(), Sensor.Value = col_number(), 
                                   Tag_ID = col_character(), Year = col_character()))

detect<-filter(detect,Tag_Type=="V9")

# Don't include detections at non-adjacent receivers with just a single detection
detect<-detect%>%
  arrange(Tag_ID,Date.and.Time)%>%
  group_by(Tag_ID)%>%
  mutate(last=ifelse(is.na(lead(Array)),1,0),
         move=match(lag(Array),LETTERS)-match(Array,LETTERS),
         nadj=ifelse(abs(move)>1,1,0))%>%
  filter(!(nadj>0 & last==1))

# Setting up detections like they want
detect<-detect%>%
  ungroup()%>%
  mutate(Receiver_Frequency=69,
         CodeSpace=paste(Receiver_Type,Receiver_Frequency,sep="-"),
         CodeSpace=gsub("W","",CodeSpace),CodeSpace=gsub("AR","",CodeSpace),
         Sensor.Value=1,Sensor.Unit="D",Signal=as.integer(Fish_ID))%>%
  dplyr::select(Timestamp=Date.and.Time,Receiver=Receiver_ID,CodeSpace,Signal,Sensor.Value,Sensor.Unit)

# Deployment df
deployments <- read_csv("deployments.csv", 
                        col_types = cols(Receiver = col_character(), 
                                         Start = col_datetime(format = "%m/%d/%Y %H:%M"), 
                                         Stop = col_datetime(format = "%m/%d/%Y %H:%M")))
# For all 2019 (and unrecovered 2018) moorings, stop date listed at 9-1-2020

#Spatial df
spatial<-read.csv("spatial.csv")%>%
  mutate(Range=130)
#Location for B2 (which didn't get recovered) was adjusted to it's deployment point from 2017

#Need to be in the same wd as spatial.csv
receiver_distMat<-read.csv("distances.csv")
rownames(receiver_distMat)<-receiver_distMat$X
receiver_distMat<-receiver_distMat[,-1]

#Table of array connections
dot<-c("A0--A--B--C--D--E--F--G")

#Transition layer
t.layer<-transitionLayer(actel::loadShape(shape="tremblay_and_Islands.shp",size=0.001,buffer=0),directions=16)

#Preload from actel package
preload<-actel::preload(biometrics = bio, deployments = deployments, spatial = spatial,
                        detections = detect, dot=dot, distances = receiver_distMat, tz="America/Iqaluit")


#Explore first (quickly get a summary of your data)
explored<-actel::explore(tz="America/Iqaluit",datapack=preload,
                         max.interval=44444,speed.warning=10,jump.warning=1,inactive.warning=400,
                         report=F,auto.open=T,discard.orphans=T,discard.first=0,GUI="needed",print.releases=T)
#New movement events are created when a fish either: moves to a new gate (array) or is not detected for ~1 month (44,444 minutes)
n #Responses to questions from actel to save output
n


#End with runRSP
RSP<-runRSP(input=explored,t.layer=t.layer,coord.x="Longitude",coord.y="Latitude",
            time.step=1440,min.time=1440,max.time=24*400,distance=500) 
#Time step=the interval at which the RSP is calculated (when min time gets exceeded between detections) (in minutes) if not moving
#Distance=the distance between RSP positions when min time is exceeded between detections at different receivers
#Min time=Shortest gap in detections for an RSP position to be calculated (in minutes)
#Max time=Maximum time between detections where an RSP will be calculated (in hours)

RSPdetections<-bind_rows(RSP[["detections"]],.id="id")%>%
  mutate(idTrack=paste(id,Track),
         season=ifelse(Date>="2017-07-03"&Date<"2017-10-17","Summer",
                       ifelse(Date>="2018-07-04"&Date<"2018-10-22","Summer",
                              ifelse(Date>="2019-06-28"&Date<"2019-11-09","Summer",
                                     ifelse(Date>="2020-07-02","Summer","Winter")))))%>%
  group_by(id)%>%
  mutate(season_year=ifelse(season=="Summer",year(Date),
                            ifelse(yday(Date)>214,year(Date),year(Date)-1)))%>%
  ungroup()


#1 day COAs for all the original detection points
coaOut<-RSPdetections%>%
  filter(Position=="Receiver")%>%
  mutate(TimeStamp.coa=floor_date(Timestamp,unit="1 day"))%>%
  group_by(Signal, TimeStamp.coa) %>%
  summarise(Latitude.coa = mean(Latitude, na.rm = TRUE),
            Longitude.coa = mean(Longitude, na.rm = TRUE),
            Season=season,
            Year=season_year)%>%distinct()%>%ungroup()%>%
  mutate(source="COA")
#Median RSP position for each 1 day (same interval as the COAs)
rspOut<-RSPdetections%>%
  filter(Position=="RSP")%>%
  mutate(TimeStamp.coa=floor_date(Timestamp,unit="1 day"))%>%
  group_by(Signal, TimeStamp.coa)%>%
  mutate(n=row_number(),med=median(n))%>%
  filter(abs(n-med)<=0.5)%>%
  group_by(Signal,TimeStamp.coa)%>%
  mutate(Latitude.coa = mean(Latitude, na.rm = TRUE),
         Longitude.coa = mean(Longitude, na.rm = TRUE))%>%
  dplyr::select(Signal,TimeStamp.coa,Season=season,Year=season_year,Latitude.coa,Longitude.coa)%>%distinct()%>%
  mutate(source="RSP")

coarspOut<-bind_rows(coaOut,rspOut)%>%
  arrange(Signal,TimeStamp.coa)%>%
  group_by(Signal,TimeStamp.coa)%>%
  mutate(N=n())%>%
  filter(N==1 | (N>1 & source=="COA"))



# I needed to modify the code for the homerange function from the latticeDensity package so that I could save the output to an object

homerange.custom <- 
  function(densityOut, percent = 0.95, output=FALSE){
    #
    if(class(densityOut)!="densityOut"){
      stop("Should be the output from the function createDensity")}
    nodes <- densityOut$nodes
    poly <- densityOut$poly
    area <- densityOut$area
    z <- densityOut$probs
    cmstz <- cumsum(sort(z))
    count <- sum(cmstz <= (1-percent))
    ind <- (z>sort(z)[count])
    plot(nodes,cex=0.1,main="Homerange")
    points(nodes[ind,],pch=19,cex=0.6,col="firebrick2")
    lines(rbind(poly,poly[1,]))
    #
    #  Compute proportion of total area in homerange
    #
    proportion <- sum(ind)/length(nodes[,1])
    prop.area <- proportion*area
  }


#Tremblay outline
Polygon <- read.csv("Points.csv")
colnames(Polygon) <- c('x', 'y')
plot(Polygon)

## converting coord system of points to epsg:3995
PConv <- Polygon
colnames(PConv) <- c('lon','lat')
coordinates(PConv) <- c('lon', 'lat')
proj4string(PConv) <- CRS("+proj=longlat +ellps=WGS84")
CRS.new <- CRS("+init=epsg:3413")
Pconv2 <- spTransform(PConv, CRS.new)
plot(PConv)
TransPoints <- data.frame(Pconv2)

## changing name from lon lat to x y ##

colnames(TransPoints) <- c('x','y')
plot(TransPoints)

NewNodeFill <- nodeFilling(poly = TransPoints, node_spacing = 300)
NewFormLattice <- formLattice(NewNodeFill)
plot(NewFormLattice)

seasons <- c("Summer","Winter")
years <- c(2017,2018,2019,2020)
IDs <- sort(unique(coaOut$Signal))


coordinates(coarspOut) <- c('Longitude.coa', 'Latitude.coa')
proj4string(coarspOut) <- CRS("+proj=longlat +ellps=WGS84")
cRSPs <- spTransform(coarspOut, CRS.new)
cRSPs<-data.frame(cRSPs)

#All fish combined home range
pdatayear <- (cRSPs[,c("Longitude.coa","Latitude.coa")])
colnames(pdatayear) <- c("x", "y")
point <- as.matrix(distinct(pdatayear))

setwd("../../../Thesis Figs/Manuscripts/Broad Movements/Lattice Density HRs/")

density <- try(createDensity(NewFormLattice, PointPattern = point,
                             k = 15, intensity = FALSE, sparse = TRUE))
if(is(density, "try-error")){
  HR <- 0
  
} else {
  tiff(paste("allfish_HR_", "coaSculpin",".tiff", sep=""),
       compression = "lzw", height=8, width=12, res=600, units="in")
  
  HR <-  homerange.custom(density, percent = 0.95, output = FALSE)
  polygon(TransPoints)
  
  dev.off()
}
HR*0.000001
HR/areaRegion(NewFormLattice)


#doing all the sculpin individually
sculp.HR <- expand.grid(IDs, seasons, years)
sculp.HR <- sculp.HR[order(sculp.HR$Var1),]
colnames(sculp.HR)<-c("Fish_ID","Season","Year")
sculp.HR[,"homerange.50"] <- NA
sculp.HR[,"homerange.57.5"] <- NA
sculp.HR[,"homerange.65"] <- NA
sculp.HR[,"homerange.72.5"] <- NA
sculp.HR[,"homerange.80"] <- NA
sculp.HR[,"homerange.87.5"] <- NA
sculp.HR[,"homerange.95"] <- NA

#Only keep those ones that have detections in them (there's no need to calculate a homerange without detections)
actual<-dplyr::select(RSPdetections,Signal,season,season_year)%>%unique()%>%
  rename(Fish_ID=Signal,Season=season,Year=season_year)
sculp.HR<-right_join(sculp.HR,actual)


homerange.custom <- 
  function(densityOut, percent = 0.95, output=FALSE){
    #
    if(class(densityOut)!="densityOut"){
      stop("Should be the output from the function createDensity")}
    nodes <- densityOut$nodes
    poly <- densityOut$poly
    area <- densityOut$area
    z <- densityOut$probs
    cmstz <- cumsum(sort(z))
    count <- sum(cmstz <= (1-percent))
    ind <- (z>sort(z)[count])
    plot(nodes,cex=0.1,main=paste(100*pers[k],"% Homerange: ",year$Signal[1]," ",year$Season[1]," ",year$Year[1],sep=""),
         xlab="Longitude",ylab="Latitude") #Change this for the main title desired
    points(nodes[ind,],pch=19,cex=0.6,col="firebrick2")
    lines(rbind(poly,poly[1,]))
    #
    #  Compute proportion of total area in homerange
    #
    proportion <- sum(ind)/length(nodes[,1])
    prop.area <- proportion*area
  }


## this loop is for the homerange size
pers<-seq(0.5,0.95,by=0.075)


for(j in 1:nrow(sculp.HR)) {
  for (k in 1:length(pers)) {
    
    ID <- subset(cRSPs, Signal == sculp.HR$Fish_ID[j])
    season <- subset(ID, Season == sculp.HR$Season[j])
    year <- subset(season, Year == sculp.HR$Year[j])
    
    pdatayear <- year[,c("Longitude.coa","Latitude.coa")]
    colnames(pdatayear) <- c("x", "y")
    point <- as.matrix(pdatayear)
    
    density <- try(createDensity(NewFormLattice, PointPattern = point,
                                 k = 15, intensity = FALSE, sparse = TRUE))
    if(is(density, "try-error")){
      sculp.HR[j,4] <- 0
      
    } else {
      tiff(paste("HR", sculp.HR$Fish_ID[j], sculp.HR$Season[j], sculp.HR$Year[j],pers[k],".tiff", sep="_"),
           compression = "lzw", height=5, width=10, res=600, units="in")
      
      sculp.HR[j,3+k] <-  homerange.custom(density, percent = pers[k])
      polygon(TransPoints)
      
      dev.off()
    }
  }
  print(j/nrow(sculp.HR)*100)
}


sculp.HR_long<-pivot_longer(sculp.HR,cols=c("homerange.50":"homerange.95"),names_to="Percent",values_to="Area_km2")%>%
  na.omit()%>%
  mutate(Percent=as.numeric(substr(Percent,11,14)),
         Area_km2=Area_km2*0.000001)

max(sculp.HR_long$Area_km2) #Area of the largest home range for an individual in a season
max(sculp.HR_long$Area_km2)/(areaRegion(NewFormLattice)*0.000001)*100 #Percentage of Tremblay that the largest homerange is
areaRegion(NewFormLattice)*0.000001 #Total area of Tremblay

seasonHR<-sculp.HR_long%>%group_by(Season)%>%summarise(mean=mean(Area_km2),se=sd(Area_km2)/sqrt(n()))
seasonHR$mean[seasonHR$Season=="Summer"]/seasonHR$mean[seasonHR$Season=="Winter"]
seasonpercentHR<-sculp.HR_long%>%group_by(Season,Percent)%>%summarise(mean=mean(Area_km2),se=sd(Area_km2)/sqrt(n()))
relChange<-seasonpercentHR%>%
  pivot_wider(id_cols=c(Percent),names_from = Season,values_from=mean)%>%
  group_by(Percent)%>%
  summarise(relChange=(Summer-Winter)/Summer)
plot(relChange)
mean(relChange$relChange)*100
sd(relChange$relChange)/sqrt(nrow(relChange))*100
sculp.HR_long%>%group_by(Season)%>%summarise(max=max(Area_km2),min=min(Area_km2))%>%unique()

#Changing the mass for 1632 that's probably off by an order of magnitude
tags$Mass<-ifelse(tags$Fish_ID=="1632",tags$Mass*10,tags$Mass)
tags$Burden<-ifelse(tags$Tag_Type=="V5",0.65/tags$Mass,3.6/tags$Mass)
sculp.HR_long2<-sculp.HR_long%>%
  mutate(Fish_ID=as.character(Fish_ID))%>%
  left_join(tags,by="Fish_ID")

summary(latHR.GLMER<-glmer(Area_km2~Season*arm::rescale(Percent)+(1|Fish_ID),
                           data=sculp.HR_long2,family=Gamma(link="log"),na.action=na.fail))
summary(latHR.GLMER.b<-glmer(Area_km2~Season*arm::rescale(Percent)*arm::rescale(Burden)+(1|Fish_ID),
                             data=sculp.HR_long2,
                             family=Gamma(link="log"),na.action=na.fail))
summary(latHR.GLM<-glm(Area_km2~Season*arm::rescale(Percent)*arm::rescale(Burden),
                       data=sculp.HR_long2,family=Gamma(link="log"),na.action=na.fail))
anova(latHR.GLMER.b,latHR.GLMER,latHR.GLM)


library(MuMIn)
d<-dredge(latHR.GLMER.b)
d
glme.95ave <- model.avg(d, cumsum(weight) <= 0.95,fit=T)
glme.95ave
#Subset models to get coefficients
summary(latHR.GLMER.b<-glmer(Area_km2~Season+arm::rescale(Percent)+arm::rescale(Burden)+
                               Season:arm::rescale(Percent)+arm::rescale(Percent):arm::rescale(Burden)+arm::rescale(Burden):Season+(1|Fish_ID),
                             data=sculp.HR_long2,
                             family=Gamma(link="log"),na.action=na.fail))
#Top model
summary(latHR.GLMER.b<-glmer(Area_km2~Season+arm::rescale(Percent)+arm::rescale(Burden)+
                               arm::rescale(Percent):arm::rescale(Burden)+arm::rescale(Burden):Season+(1|Fish_ID),
                             data=sculp.HR_long2,
                             family=Gamma(link="log"),na.action=na.fail))

predOutSum.b<-cbind(sculp.HR_long2,pred=predict(latHR.GLMER.b,type="response"))%>%
  group_by(Season,Fish_ID,Burden)%>%
  summarise(mean=mean(pred),
            lower=mean-1.96*sd(pred),upper=mean+1.96*sd(pred))
predOutSum.p<-cbind(sculp.HR_long2,pred=predict(latHR.GLMER,type="response"))%>%
  group_by(Season,Percent)%>%
  summarise(mean=mean(pred),
            lower=mean-1.96*sd(pred),upper=mean+1.96*sd(pred))
predOut<-cbind(sculp.HR_long2,pred=predict(latHR.GLMER,type="response"))




#With mass/burden as the size
set.seed(42)
#Full lat density percentage plots
ggplot(predOut,aes(Percent,Area_km2,fill=Season,size=Burden))+
  scale_x_continuous(expand=expansion(mult=c(0.05,0.05)),name="Home Range Contour",
                     breaks=c(50,57.5,65,72.5,80,87.5,95),labels=c("50%","57.5%","65%","72.5%","80%","87.5%","95%"))+
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,27.5),expand=expansion(mult=0.001),
                     name=expression(paste("Lattice Density Homerange (km"^"2"*")")))+
  #theme(legend.position = c(0.5,0.2),legend.background = element_rect(color="black"))+
  guides(fill=guide_legend(override.aes=list(shape=21,size=10,stroke=1.1)),
         size=guide_legend(override.aes=list(stroke=c(0.5,1,1.5,2))))+
  theme(legend.position=c(0.29,0.7),legend.background=element_rect(color="black",size=1),
        legend.box="horizontal")+
  scale_fill_viridis_d(begin=0.2,name="Season",labels=c("Ice-Free","Ice-Covered"),direction=-1)+
  scale_color_viridis_d(begin=0.2,name="Season",direction=-1,guide="none")+
  scale_size_continuous(range=c(1,8))+
  geom_boxplot(aes(group=paste(Percent,Season)),outlier.shape=NA,show.legend = F,size=0.95,color="black")+
  geom_point(aes(stroke=Burden*20),position=position_jitterdodge(0.95,dodge.width=5),shape=21)+
  #geom_line(data=predOut,aes(Percent,pred,group=Fish_ID,color=Season),inherit.aes=F)+
  geom_ribbon(data=predOutSum.p,aes(Percent,ymin=lower,ymax=upper,color=Season,fill=Season),
              alpha=0.33,inherit.aes=F,show.legend = F)+
  geom_line(data=predOutSum.p,aes(Percent,mean,color=Season),inherit.aes=F,lwd=2.5)




