scope.data<-read.csv("extracted_data_CDA_Wilson_2022_final set.csv")

require(ggplot2)
require(ggpubr)

#How many articles either 1) describe the number of years that a healthcare worker worked in the healthcare industry before asthma onset,
#2) describe the number of years working in the healthcare industry and its relationship with asthma rates among participants, or
#3) described or mentioned the word "latency"?

num.articles.latency<-length(scope.data$Timestamp[scope.data$Does.the.article.include.the.term...latency..=="Yes" | 
                            scope.data$Does.the.article.address.the.number.of.years.working.in.healthcare.industry.before.asthma.symptom.onset.=="Yes"|
                            scope.data$Does.the.article.address.the.number.of.years.working.in.the.healthcare.industry.and.asthma.rates.outcome.among.participants.=="Yes" ])

# % of articles with latency info
num.articles.latency/length(scope.data$Timestamp)*100


View(scope.data[scope.data$Does.the.article.include.the.term...latency..=="Yes" | 
                              scope.data$Does.the.article.address.the.number.of.years.working.in.healthcare.industry.before.asthma.symptom.onset.=="Yes"|
                              scope.data$Does.the.article.address.the.number.of.years.working.in.the.healthcare.industry.and.asthma.rates.outcome.among.participants.=="Yes" ,])

#How many articles provided information on 1) frequency of cleaning/disifnection, concentrations of exposure (including in animal studies), information on relationships
#between frequency or concentration of exposure and an asthma-related outcome, and/or included the term dose-response?

num.articles.exposure.freq<-length(scope.data$Timestamp[scope.data$Does.it.address.measurements..frequency.and.or.concentration..of.cleaning.disinfection.exposures.during.cleaning.and.disinfection.in.healthcare.environments...or.animal.studies.=="Yes"|
                                                          scope.data$Does.the.article.address.frequency.of.cleaning.and.disinfection.activity.=="Yes"|
                                                          scope.data$Does.the.article.include.likelihood.or.asthma.rate.as.a.function.of.frequency.of.cleaning.and.disinfection.activity.=="Yes"|
                                                          scope.data$Does.the.article.include.the.phrase...dose.response..=="Yes"|
                                                          scope.data$Does.the.article.address.concentration.of.an.exposure.to.cleaning.and.disinfection.=="Yes"])

View(scope.data[scope.data$Does.it.address.measurements..frequency.and.or.concentration..of.cleaning.disinfection.exposures.during.cleaning.and.disinfection.in.healthcare.environments...or.animal.studies.=="Yes"|
                            scope.data$Does.the.article.address.frequency.of.cleaning.and.disinfection.activity.=="Yes"|
                            scope.data$Does.the.article.include.likelihood.or.asthma.rate.as.a.function.of.frequency.of.cleaning.and.disinfection.activity.=="Yes"|
                            scope.data$Does.the.article.include.the.phrase...dose.response..=="Yes"|
                            scope.data$Does.the.article.address.concentration.of.an.exposure.to.cleaning.and.disinfection.=="Yes",])

# % of articles with exposure/freq info
num.articles.exposure.freq/length(scope.data$Timestamp)*100

#---------------latency data fit-----------------------------------------------------------------------------------------------

#NEED TO DOUBLE-CHECK THIS FILE
latency<-read.csv("latency data_final.csv")
require(ggplot2)
require(ggpubr)
library(fitdistrplus)
library(MASS)

exp.fit<-fitdist(latency$Years[latency$Outcome=="Occupational asthma" & latency$Years!=0],"exp")
#weibull.fit<-fitdist(latency$Years[latency$Outcome=="Occupational asthma" & latency$Years!=0],"weibull",lower=c(0,0))
fit<-gofstat(exp.fit)
fit$cvmtest
fit$adtest
fit$kstest

#denscomp(list(exp.fit,weibull.fit),legendtext=c("Exponential","Weibull"),xlab="Latency (Years)",ylab="Density",fitcol=c("black","black"),
#         main="",fitlwd = 3,breaks=10)
#legend(10,0.25,legend=c("Exponential","Weibull"),col=c("black","black"),lty=c(1,2),lwd=3)

denscomp(list(exp.fit),legendtext=c("Exponential"),xlab="Latency (Years)",ylab="Density",fitcol=c("black","black"),
         main="",fitlwd = 3,breaks=10)
legend(7,0.18,legend=c("Exponential"),col=c("black"),lwd=3)


A<-recordPlot()

#only OA
latency<-latency[latency$Outcome=="Occupational asthma",]
#A<-ggplot(latency)+geom_histogram(aes(x=Years),color="black",fill="grey")+
#  scale_x_continuous(name="Latency (Years)")+
#  scale_y_continuous(name="Count")+
#  theme_bw()+
#  theme(axis.title = element_text(size=16),axis.text=element_text(size=16))

B<-ggplot(latency[latency$Outcome=="Occupational asthma" & latency$Years!=0,])+geom_histogram(aes(x=Years),color="black",fill="grey")+
  scale_x_continuous(name="Latency (Years)")+
  scale_y_continuous(name="Count")+
  facet_wrap(~Chemical,ncol=3)+
  theme_bw()+
  theme(axis.title = element_text(size=16),axis.text=element_text(size=16),strip.text = element_text(size=13))

ggarrange(A,B,ncol=2)

summary(latency$Years[latency$Outcome=="Occupational asthma" & latency$Years!=0])

#---------------Dose-response---------------------------------------------------------------------------------

measured_conc<-read.csv('Copy of measured_concentrations_v3.csv')

measured_conc$Concentration[measured_conc$Concentration=="Below LOD"]<-measured_conc$LOD[measured_conc$Concentration=="Below LOD"]/sqrt(2)
measured_conc$Concentration<-as.numeric(measured_conc$Concentration)

measured_conc$Concentration[measured_conc$Units=="ug/m^3"]<-measured_conc$Concentration[measured_conc$Units=="ug/m^3"]*(1/1000)
measured_conc$Units[measured_conc$Units=="ug/m^3"]<-"mg/m^3"

measured_conc$Concentration[measured_conc$Units=="ppb"]<-measured_conc$Concentration[measured_conc$Units=="ppb"]*(1/1000)
measured_conc$Units[measured_conc$Units=="ppb"]<-"ppm"

#----------convert mg/m^3 to ppm (assumes 25 deg C and 1 atm)-----------------------------------------------

#Entering molecular weights

#removing TVOC11 and TVOC14 since they are not comparable to other studies and didn't specify the mixture
measured_conc<-measured_conc[measured_conc$Chemical!="TVOC11" & measured_conc$Chemical!="TVOC14",]

#Assume glutaraldialdehyde = glutaraldehyde
measured_conc$Chemical[measured_conc$Chemical=="Glutaradialdehyde"]<-"Glutaraldehyde"

#Remove mixtures & TVOC (TVOC not specified)
measured_conc<-measured_conc[measured_conc$Chemical!="Hydrogen peroxcide, peracetic acid, acetic acid"&
                               measured_conc$Chemical!="Hydrogen peroxide and peraceteic acid",]

chemicalnames<-c("2-Propanol","Ethanol","Glutaraldehyde","didecyl-dimethylammonium chloride (DDAC)","Formaldehyde",
                 "Succinaldehyde","Acetic acid","Acetone","alpha-Pinene","Chloroform","d-Limonene",
                 "ortho-phthalaldehyde (OPA)","Peracetic acid","Benzene","Toluene","Ethylbenzene","m,p-xylene",
                 "o-xylene")

#Current uncertainty regarding mixture for TVOC and method of collection https://www.cdc.gov/niosh/nioshtic-2/20043503.html
molecular.weights<-c(60.1,46.07,100.11,362.08,30.031,86.09,60.052,58.08,136.23,119.38,136.23,134.13,76.0514,
                     78.11,92.14,106.167,106.16,106.16)

frame.molecular.weight<-data.frame(chemicalnames,molecular.weights)

limits<-read.csv('limits.csv')


measured_conc$limits<-NA
measured_conc$limittype<-NA

for (i in 1:length(chemicalnames)){
  measured_conc$Concentration[measured_conc$Chemical==chemicalnames[i] & measured_conc$Units!="ppm"]<-measured_conc$Concentration[measured_conc$Chemical==chemicalnames[i] & measured_conc$Units!="ppm"]*24.45/molecular.weights[i]
  measured_conc$Units[measured_conc$Chemical==chemicalnames[i] & measured_conc$Units!="ppm"]<-"ppm"
  #measured_conc$limits[measured_conc$Chemical==chemicalnames[i]]<-limits$Limit[limits$Chemical==chemicalnames[i]]
  #measured_conc$limittype[measured_conc$Chemical==chemicalnames[i]]<-limits$Limit.Type[limits$Chemical==chemicalnames[i]]
}

windows()
ggplot(measured_conc)+geom_histogram(aes(x=Concentration,y=..density..,group=Chemical),alpha=0.3)+
  geom_density(aes(x=Concentration,group=Chemical),alpha=0.3)+
  geom_vline(data=limits[!is.na(limits$Limit),],aes(xintercept=Limit,color=Limit.Type),size=1.5,linetype="dashed")+
  #geom_vline(data=measured_conc[!is.na(measured_conc$limits),],aes(xintercept=limits,color=limittype),size=1.5,linetype="dashed")+
  scale_color_discrete(name="")+
  theme(legend.position = "top")+
  scale_x_continuous(trans="log10",name="Concentration (ppm)")+
  scale_y_continuous(name="Density")+
  scale_fill_discrete(name="")+
  theme_pubr()+
  facet_wrap(~Chemical)+
  #ggtitle("B")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),legend.text = element_text(size=14),legend.position = "top")

#------------more dose-response checking-------------------------------------------------------------------

scope.data<-read.csv("extracted_data_CDA_Wilson_2022_final set.csv")

#Does it address frequency?
dose.response.sub.freq<-scope.data[scope.data$Does.the.article.address.frequency.of.cleaning.and.disinfection.activity.=="Yes" ,]
write.csv(dose.response.sub.freq,"frequency.csv")

#loading in csv file with OR data
data.frequency<-read.csv("freq_and_asthma_data_extract.csv")

#removing reference groups just for plotting

data.frequency<-data.frequency[data.frequency$Value.of.Outcome!=1 & data.frequency$X95..CI.Lower!=1 & data.frequency$X95..CI.Upper!=1,]

require(ggh4x)
require(ggplot2)
require(ggpubr)

A<-ggplot(data.frequency[data.frequency$Outcome.2=="ACT score",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(Outcome~Type.of.Outcome)+
  theme_pubr()+
  ggtitle("ACT Score")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

B.1<-ggplot(data.frequency[data.frequency$Outcome.2=="Asthma" & data.frequency$Outcomev2=="Asthma",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(~Type.of.Outcome,nrow=1)+
  theme_pubr()+
  ggtitle("Asthma")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

B.2<-ggplot(data.frequency[data.frequency$Outcome.2=="Asthma" & data.frequency$Outcomev2!="Asthma" & 
                             data.frequency$Outcomev2=="Work-related asthma" | 
                           data.frequency$Outcomev2=="Occupational asthma" | data.frequency$Outcomev2=="Post-hire asthma",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(~Type.of.Outcome,nrow=1)+
  theme_pubr()+
  ggtitle("Asthma")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

B.3<-ggplot(data.frequency[data.frequency$Outcome.2=="Asthma" & data.frequency$Outcomev2!="Asthma" & 
                             data.frequency$Outcomev2!="Work-related asthma" & 
                             data.frequency$Outcomev2!="Occupational asthma" & data.frequency$Outcomev2!="Post-hire asthma",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure,Outcomev2),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(~Type.of.Outcome,nrow=1)+
  theme_pubr()+
  ggtitle("Asthma")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

C<-ggplot(data.frequency[data.frequency$Outcome.2=="BHR",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(Outcome~Type.of.Outcome)+
  theme_pubr()+
  ggtitle("BHR")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

D<-ggplot(data.frequency[data.frequency$Outcome.2=="COPD",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(Outcome~Type.of.Outcome)+
  theme_pubr()+
  ggtitle("COPD")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

E<-ggplot(data.frequency[data.frequency$Outcome.2=="Exacerbation",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(Outcome~Type.of.Outcome)+
  theme_pubr()+
  ggtitle("Exacerbation")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))

G<-ggplot(data.frequency[data.frequency$Outcome.2=="Respiratory symptoms",])+
  geom_point(aes(x=interaction(Categories,Frequency.Measure,Outcome),color=X1st.Author,y=Value.of.Outcome),size=3)+
  geom_errorbar(aes(x=interaction(Categories,Frequency.Measure,Outcome),color=X1st.Author,ymin=X95..CI.Lower,ymax=X95..CI.Upper),size=1)+
  geom_hline(yintercept=1,linetype="dashed",color="black",size=1)+
  scale_x_discrete(guide="axis_nested",name="")+
  scale_y_continuous(name="")+
  scale_color_discrete(name="1st Author")+
  coord_flip() +
  facet_wrap(~Type.of.Outcome)+
  theme_pubr()+
  ggtitle("Respiratory symptoms")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),strip.text=element_text(size=13),
        legend.text = element_text(size=13),title=element_text(size=14))
