MAMBI.DJG <- function(BenthicData, EG_File_Name="Ref - EG Values 2018.csv", EG_Scheme="Hybrid") {
	require(tidyverse)
	require(reshape2)
	require(vegan)
	require(readxl)
	source("EQR.R")

	Saline_Standards <- read_xlsx("Pelletier2018_Standards.xlsx", sheet = "Saline Sites")
	TidalFresh_Standards <- read_xlsx("Pelletier2018_Standards.xlsx", sheet = "Tidal Fresh Sites")


	Input_File<- as_tibble(BenthicData) %>% 
	  mutate(Species_ended_in_sp=(str_detect(Species," sp$")), Taxon=(str_replace(Species, " sp$",""))) %>% 
	  mutate(Coast=(ifelse(Longitude<=-115,"West","Gulf-East")))%>% 
	  mutate(SalZone=case_when(Salinity>30&Salinity<=40&Coast=="Gulf-East"~"EH", Salinity>18&Salinity<=30&Coast=="Gulf-East"~"PH", Salinity>5&Salinity<=18~"MH", 
				   Salinity>0.2&Salinity<=5~"OH", Salinity<=0.2~"TF", Salinity>40~"HH",
				   Salinity>30&Salinity<=40&Coast=="West"~"WEH", Salinity>18&Salinity<=30&Coast=="West"~"WPH")) 

	  

	EG_Ref<-read.csv(EG_File_Name, stringsAsFactors = F, na.strings = "") %>% select(.,Taxon, Exclude, EG=EG_Scheme) %>% mutate(EG=(ifelse(Taxon=="Oligochaeta", "V", EG))) 




	total.abundance<-Input_File %>% group_by(StationID, Replicate, SampleDate) %>% summarise(Tot_abun=sum(Abundance))

	Sample.info<-Input_File %>% select(StationID, Replicate, SampleDate, Latitude, Longitude, Salinity, Coast, SalZone) %>% unique()

	Input_File2<-Input_File %>% filter(!is.na(Salinity))

	EG.Assignment<-Input_File %>% left_join(., EG_Ref, by="Taxon") %>% #filter(Exclude!="Yes") %>% 
	  left_join(.,total.abundance, by=c("StationID", "Replicate", "SampleDate")) %>% mutate(Rel_abun=((Abundance/Tot_abun)*100))

	AMBI.applicability<-EG.Assignment %>% mutate(EG_Test=ifelse(is.na(EG),"NoEG", "YesEG")) %>% dcast(.,StationID+Replicate+SampleDate~EG_Test, value.var = "Rel_abun", fun.aggregate = sum) %>% 
	   mutate(Use_AMBI=case_when(NoEG<=20~"Yes", NoEG>20&NoEG<=50~"With Care", NoEG>50~"Not Recommended", is.na(NoEG)~"Yes"))

	MAMBI.applicability<-Sample.info %>% mutate (Use_MAMBI=ifelse(is.na(SalZone),"No - No Salinity","Yes")) %>% 
	  select(StationID, Replicate, SampleDate, Use_MAMBI)

	Sal_range.dataset<-unique(Input_File2$SalZone)

	AMBI.Scores<-EG.Assignment %>% group_by(StationID, Replicate, SampleDate,Tot_abun,EG) %>% summarise(Sum_Rel=sum(Rel_abun)) %>% replace_na(list(EG="NoEG")) %>% 
		mutate(EG_Score= case_when(EG=="I"~Sum_Rel*0, EG=="II"~Sum_Rel*1.5, EG=="III"~Sum_Rel*3, EG=="IV"~Sum_Rel*4.5, EG=="V"~Sum_Rel*6, EG=="NoEG"~0)) %>% 
	 	mutate(EG_Score=ifelse(Tot_abun==0,700,EG_Score)) %>% 
	  	group_by(StationID, Replicate, SampleDate) %>% summarise(AMBI_Score=(sum(EG_Score)/100))

	Rich<-Input_File %>% group_by(StationID, Replicate, SampleDate) %>% summarise(S=length(Taxon))
	Rich$S<-as.numeric(Rich$S)
	Divy<-Input_File %>% dcast(StationID+Replicate+SampleDate~Taxon, value.var = "Abundance", fill=0) %>%
		mutate(H=diversity((select(.,4:(ncol(.)))),index = "shannon", base = 2)) %>% select(.,StationID, Replicate, SampleDate,H)


	 
	metrics<-AMBI.Scores %>% left_join(.,Rich, by=c("StationID", "Replicate", "SampleDate")) %>% left_join(.,Divy, by=c("StationID", "Replicate", "SampleDate")) %>% 
		mutate(S=(ifelse(AMBI_Score==7,0,S)),H=(ifelse(AMBI_Score==7,0,H)))
	metrics.1<-Sample.info %>% left_join(., metrics, by=c("StationID", "Replicate", "SampleDate")) %>% 
		select(StationID, Replicate, SampleDate,AMBI_Score, S, H, SalZone)


	metrics.2<-bind_rows(metrics.1,Saline_Standards)


	saline.mambi<-map_dfr(Sal_range.dataset, function(sal){
	sal.df<-filter(metrics.2, SalZone==sal)
	METRICS.tot<-sal.df[,c(4:6)]


	  options(warn = -1)
	    METRICS.fa2 <- princomp(METRICS.tot, cor = T, covmat = cov(METRICS.tot))
	   options(warn = 0)
	    METRICS.fa2.load <- loadings(METRICS.fa2) %*% diag(METRICS.fa2$sdev)
	    METRICS.fa2.load.varimax <- loadings(varimax(METRICS.fa2.load))
	    METRICS.scores2 <- scale(METRICS.tot) %*% METRICS.fa2.load.varimax
	    colnames(METRICS.scores2) <- c("x", "y", "z")
	  METRICS.tr <- METRICS.scores2
	  


	eqr <-EQR(METRICS.tr)
	colnames(eqr)<-c("MAMBI_Score")
	eqr<-data.frame(eqr)

	results<-sal.df %>% bind_cols(.,eqr) %>% left_join(.,Sample.info, by=c("StationID", "Replicate", "SampleDate", "SalZone")) %>% 
	  select(1,2,3,9,10,7,4:6,8) %>% filter(!StationID%in%Saline_Standards$StationID, SalZone!="TF") %>% 
	  mutate(Orig_MAMBI_Condition=case_when(MAMBI_Score<0.2~"Bad", MAMBI_Score>=0.2&MAMBI_Score<0.39~"Poor", MAMBI_Score>=0.39&MAMBI_Score<0.53~"Moderate", MAMBI_Score>=0.53&MAMBI_Score<0.77~"Good", MAMBI_Score>=0.77~"High"),
		 New_MAMBI_Condition=case_when(MAMBI_Score<=0.387~"High Disturbance", MAMBI_Score>0.387&MAMBI_Score<0.483~"Moderate Disturbance", MAMBI_Score>=0.483&MAMBI_Score<0.578~"Low Disturbance",MAMBI_Score>=0.578~"Reference"))
	saline.mambi<-results
	})

if(any(Sal_range.dataset=="TF"))
{
  
  TF.EG.Assignment <- EG.Assignment %>% filter(SalZone=="TF")
  TF.EG_Ref<-EG_Ref<-read.csv(EG_File_Name, stringsAsFactors = F, na.strings = "") %>% select(.,Taxon, Exclude, EG=EG_Scheme, Oligochaeta)
  
TF.AMBI.Scores<-TF.EG.Assignment %>% group_by(StationID, Replicate, SampleDate,Tot_abun,EG) %>% summarise(Sum_Rel=sum(Rel_abun)) %>% replace_na(list(EG="NoEG")) %>% 
  mutate(EG_Score= case_when(EG=="I"~Sum_Rel*0, EG=="II"~Sum_Rel*1.5, EG=="III"~Sum_Rel*3, EG=="IV"~Sum_Rel*4.5, EG=="V"~Sum_Rel*6, EG=="NoEG"~0)) %>% 
  mutate(EG_Score=ifelse(Tot_abun==0,700,EG_Score)) %>% 
  group_by(StationID, Replicate, SampleDate) %>% summarise(AMBI_Score=(sum(EG_Score)/100))

TF.Oligos <- Input_File %>% left_join(., total.abundance, by =c("StationID", "Replicate", "SampleDate")) %>% left_join(., TF.EG_Ref, by="Taxon" ) %>% 
  filter(Oligochaeta=="Yes", SalZone=="TF") %>% group_by(StationID, Replicate, SampleDate) %>% 
  summarise(Oligo_pct=(sum(Abundance/Tot_abun))*100)

TF.Divy<-Input_File %>% filter(SalZone=="TF") %>% dcast(StationID+Replicate+SampleDate~Taxon, value.var = "Abundance", fill=0) %>%
  mutate(H=diversity((select(.,4:(ncol(.)))),index = "shannon", base = 2)) %>% select(.,StationID, Replicate, SampleDate,H)


TF.metrics<-TF.AMBI.Scores %>% left_join(.,TF.Divy, by=c("StationID", "Replicate", "SampleDate")) %>% left_join(.,TF.Oligos, by=c("StationID", "Replicate", "SampleDate")) %>% 
  mutate(Oligo_pct=(ifelse(AMBI_Score==7,0,Oligo_pct)),H=(ifelse(AMBI_Score==7,0,H)), Oligo_pct=(ifelse(is.na(Oligo_pct),0,Oligo_pct)))

TF.metrics.1<-Sample.info %>% left_join(., TF.metrics, by=c("StationID", "Replicate", "SampleDate")) %>% 
  select(StationID, Replicate, SampleDate,AMBI_Score, H, Oligo_pct, SalZone) %>% filter(SalZone=="TF")

TF.metrics.2<-bind_rows(TF.metrics.1,TidalFresh_Standards)

TF.METRICS.tot<-TF.metrics.2[,c(4:6)]

options(warn = -1)
TF.METRICS.fa2 <- princomp (TF.METRICS.tot, cor = T, covmat = cov(TF.METRICS.tot))
options(warn = 0)
TF.METRICS.fa2.load <- loadings(TF.METRICS.fa2) %*% diag(TF.METRICS.fa2$sdev)
TF.METRICS.fa2.load.varimax <- loadings(varimax(TF.METRICS.fa2.load))
TF.METRICS.scores2 <- scale(TF.METRICS.tot) %*% TF.METRICS.fa2.load.varimax
colnames(TF.METRICS.scores2) <- c("x", "y", "z")
TF.METRICS.tr <- TF.METRICS.scores2

TF.eqr <-EQR(TF.METRICS.tr)
colnames(TF.eqr)<-c("MAMBI_Score")
TF.eqr<-data.frame(TF.eqr)

TF.mambi<-TF.metrics.2 %>% bind_cols(.,TF.eqr) %>% left_join(.,Sample.info, by=c("StationID", "Replicate", "SalZone", "SampleDate")) %>% 
  select(1,2,3,9,10,7,4:6,8) %>% filter(!StationID%in%TidalFresh_Standards$StationID) %>% 
  mutate(Orig_MAMBI_Condition=case_when(MAMBI_Score<0.2~"Bad", MAMBI_Score>=0.2&MAMBI_Score<0.39~"Poor", MAMBI_Score>=0.39&MAMBI_Score<0.53~"Moderate", MAMBI_Score>=0.53&MAMBI_Score<0.77~"Good", MAMBI_Score>=0.77~"High"),
         New_MAMBI_Condition=case_when(MAMBI_Score<=0.387~"High Disturbance", MAMBI_Score>0.387&MAMBI_Score<0.483~"Moderate Disturbance", MAMBI_Score>=0.483&MAMBI_Score<0.578~"Low Disturbance",MAMBI_Score>=0.578~"Reference"))

saline.mambi<- saline.mambi %>% mutate(Oligo_pct=NA) %>% select(1:9,13,10,11,12)
TF.mambi<-TF.mambi %>% mutate(S=NA) %>% select(1:7,13,8:12)

Overall.Results<-bind_rows(saline.mambi, TF.mambi) %>% left_join(.,AMBI.applicability[,c(1,2,3,5,6)], by=c("StationID", "Replicate", "SampleDate")) %>% 
  left_join(MAMBI.applicability,.,by=c("StationID", "Replicate", "SampleDate")) %>% 
  select(1:3,5:14,4,16,15)
}


else
{
   Overall.Results<-saline.mambi%>% left_join(.,AMBI.applicability[,c(1,2,3,5,6)], by=c("StationID", "Replicate", "SampleDate")) %>% 
     left_join(MAMBI.applicability,.,by=c("StationID", "Replicate", "SampleDate")) %>% 
     select(1:3,5:13,4,15,14)
}


}
