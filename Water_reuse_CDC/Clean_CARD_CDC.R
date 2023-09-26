#Make sure you know which CARD Database you are using!
#make sure you are in the right folder!
setwd ("~/Desktop/CDC RW/")


#loads your libraries
library(reshape2)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(stringr)
library(ggpubr)
library(vegan)


RawData<- read_excel("CARD080620CDC.xlsx",sheet = 1)
RawData$gene<-toupper(RawData$gene)

FRSheet<-read_excel("CARD - Manual Classification.xlsx")
FRSheet$Resistance_Genes<-toupper(FRSheet$Resistance_Genes)
#Adds Gene Categories 

Alldata_MergedResistances<-merge.data.frame(RawData,
                                       FRSheet,
                                       by.y = "Resistance_Genes",
                                       by.x="gene",
                                       sort = TRUE, all.x=TRUE)

remove_args<-read.csv("args.csv")
Alldata_MergedResistances_sd<-subset(Alldata_MergedResistances, !(gene %in% remove_args$Removed.ARGs ))
Alldata_MergedResistances_sd$type.y<-str_replace(Alldata_MergedResistances_sd$type.y,"mls","MLS")

#Adds in sample data 
Sampledata<-read_excel("Sample Data_CDCProject.xlsx",sheet=2)
Alldata_MergedResistances_sd<-merge.data.frame(Alldata_MergedResistances_sd,
                                            Sampledata,
                                            by = "SampleID",
                                            sort = TRUE)

Alldata_MergedResistances_sd$type.y <-factor(Alldata_MergedResistances_sd$type.y,levels=c("fosfomycin",      "triclosan",       "phenicol"  ,      "rifamycin" ,      "sulfonamide"   ,  "fluoroquinolone",
                                                                                      "peptide"  ,       "beta-lactam"  ,   "glycopeptide"   , "aminoglycoside" , "tetracycline"  ,  "MLS"   ,         
                                                                                      "other"   ,        "multidrug"  ))

Samples.genetotals<-  ddply(Alldata_MergedResistances_sd,.(SampleID,`Sampling Trip`,Location,Plant,`Treatment Stage`, gene,type.y), summarise,
                            gc_samp= mean(LogCopiesmL,na.rm=TRUE),
                            avClassfreqrel_samp = mean(`16S_Normalization`,na.rm=TRUE))



Other<-subset(Samples.genetotals, (Plant == "Virginia" & type.y == "other" & (Location ==10 | Location ==11)))

SamplesClass.totals<-  ddply(Samples.genetotals,.(SampleID,`Sampling Trip`,Location,Plant,`Treatment Stage`, type.y), summarise,
                             GCtot_samp=sum(gc_samp),
                             freq_class_samp = sum(avClassfreqrel_samp))




SamplesClass.totals$type.y <-factor(SamplesClass.totals$type.y,levels=c("other"   ,        "multidrug" ,"fosfomycin",      "triclosan",       "phenicol"  ,      "rifamycin" ,      "sulfonamide"   ,  "fluoroquinolone",
                                                                                        "peptide"  ,       "beta-lactam"  ,   "glycopeptide"   , "aminoglycoside" , "tetracycline"  ,  "MLS"           
                                                                                         ))
#Figure 2: Profiles
#NMDS Vegan
library(vegan)
braycurt<-na.omit(reshape2::dcast(Alldata_MergedResistances_sd, SampleID +Plant +`Treatment Stage` +Location ~ type.y , value.var = "count",sum))

#test<-na.omit(reshape2::dcast(Samples.Classtotals, SampleID +Plant +`Treatment Stage` +Location ~ type.y , value.var = "sumClassfreqrel",sum))
rownames(braycurt)<-braycurt$SampleID

braydist<-vegdist(braycurt[5:15],method="bray",binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE)

nmds_occurances<-metaMDS(na.omit(braydist),scale=TRUE,trymax = 100)

#plot(nmds_occurances,type="t")
NMDSSD<-cbind(data.frame(x=nmds_occurances$points[,1],y=nmds_occurances$points[,2]),braycurt[1:4])

colors = c(brewer.pal(name="Set2",n=8),brewer.pal(name="Paired",n=8) )


NMDSSD$`Treatment Stage`<-factor(NMDSSD$`Treatment Stage`,levels=c("Influent","Primary","Act. Sludge","Secondary",
                                                                    "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination","3 days","6 months", "Background","Denitrification",
                                                                   "Cl Storage", "Short Dist.","Long Dist."))

#NMDSSD<-subset(NMDSSD, Plant== "Virginia")
Fig2<-ggplot(NMDSSD,aes(x=x,y=y,fill =`Treatment Stage`,color = `Treatment Stage`)) +
  geom_point(aes(shape= Plant),size=9,alpha= 0.5)+
  scale_shape_manual(values=c(15,17))+
  geom_text(aes(label=Location),color = "Black",size=4,position = "identity") +
  #geom_text(aes(label=SampleID),color = "Black",size=2,position = "identity") +
  scale_color_manual(values =colors)+
 theme(axis.title.y=element_blank(), axis.title.x=element_blank())

Fig2

SamplesClass.totals<- subset(SamplesClass.totals,(type.y != "NA") )

Samples.Classtotals_alltrips<-  ddply(SamplesClass.totals,.(Location,Plant,`Treatment Stage`,type.y), summarise,
                                AvClass = mean(freq_class_samp,na.rm = TRUE))
Samples.Classtotals_alltrips<-Samples.Classtotals_alltrips[apply(Samples.Classtotals_alltrips!=0, 1, all),]

Samples.Classtotals_av<-  ddply(SamplesClass.totals,.(Location,Plant,`Sampling Trip`,`Treatment Stage`), summarise,
                                      AvtotClass = sum(freq_class_samp,na.rm = TRUE))
Samples.Classtotals_av<-Samples.Classtotals_av[apply(Samples.Classtotals_av!=0, 1, all),]

Samples.Classtotals_stdev1<-  ddply(Samples.Classtotals_av,.(Location,Plant,`Treatment Stage`), summarise,
                                   StdevtotClass = sd(AvtotClass,na.rm = TRUE),
                                   AvtotClass = mean(AvtotClass,na.rm = TRUE))
Samples.Classtotals_stdev1<-na.omit(Samples.Classtotals_stdev1)
Samples.Classtotals_stdev<-subset(Samples.Classtotals_stdev1,(Plant == "Virginia") )
Samples.Classtotals_stdev$`Treatment Stage`<-factor(Samples.Classtotals_stdev$`Treatment Stage`,levels<-c("Influent","Primary","Act. Sludge","Secondary",
                                                                                                          "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination",
                                                                                                          "3 days","6 months","Background"))



AllSampClassprephelp<- subset(Samples.Classtotals_alltrips,(Plant == "Virginia") )
AllSampClassprephelp<-na.omit(AllSampClassprephelp)
AllSampClassprephelp$`Treatment Stage`<-factor(AllSampClassprephelp$`Treatment Stage`,levels<-c("Influent","Primary", "Act. Sludge","Secondary",
                                                                                                "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination",
                                                                                                "3 days","6 months", "Background"))
AllSampClassprephelp$type.y <-factor(AllSampClassprephelp$type.y,levels=c("other"   ,        "multidrug" ,"fosfomycin",      "triclosan",       "phenicol"  ,      "rifamycin" ,      "sulfonamide"   ,  "fluoroquinolone",
                                                                                          "peptide"  ,       "beta-lactam"  ,   "glycopeptide"   , "aminoglycoside" , "tetracycline"  ,  "MLS"           
                                                                                           ))

mycolors = colors <- c("seagreen2",  "mediumblue",  "mediumvioletred", "gray17", "dodgerblue", 
                       "indianred2", "darkmagenta", "green", "pink", "skyblue", 
                       "dark red", "yellow", "gray47","magenta2", "blue",  "darkorange1", "purple", "darkturquoise", 
                       "gray67", "goldenrod2")

VAAllSampGenesplotrel<- ggplot(data= AllSampClassprephelp, aes(x = `Treatment Stage`))+
 geom_bar(aes(y =AvClass,fill = type.y),stat = "identity",color="black",width=0.9)+ #factor(Resistance_Genes,levels = lvls)
  geom_errorbar(data=Samples.Classtotals_stdev, colour="black",aes(ymin=AvtotClass,
                                                                   ymax = AvtotClass+ StdevtotClass)) +
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key.size = unit(0.5, 'lines'))+
  theme(axis.title.y=element_text(size=18))+
  theme(legend.title = NULL) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1, size=20,color="black"))+
  ylab("ARG copies / 16S rRNA \n")+
  theme(axis.title.y=element_text(vjust=-1.5))+
  theme(axis.title.x=element_blank())+
  #scale_y_continuous(breaks=seq(0,1.0,0.25),limits=c(0,1.0))+
  theme(legend.text = element_text(size=22,color="black")) +
  theme(axis.text.y=element_text( size=14,color="black"))

VAAllSampGenesplotrel

AllSampClassprephelp<- subset(Samples.Classtotals_alltrips,(Plant == "Florida") )
AllSampClassprephelp<-na.omit(AllSampClassprephelp)

Samples.Classtotals_stdev<-subset(Samples.Classtotals_stdev1,(Plant == "Florida") )
AllSampClassprephelp$`Treatment Stage`<-factor(AllSampClassprephelp$`Treatment Stage`,levels=c("Influent"  ,"Primary","Act. Sludge","Secondary",
                                                                                            "Denitrification","Chlorination","Cl Storage","Short Dist.", "Long Dist."))

Samples.Classtotals_stdev$`Treatment Stage`<-factor(Samples.Classtotals_stdev$`Treatment Stage`,levels=c("Influent","Primary","Act. Sludge","Secondary",
                                                                                                         "Denitrification","Chlorination","Cl Storage","Short Dist.", "Long Dist."))

FLAllSampGenesplotrel<- ggplot(data= AllSampClassprephelp, aes(x = `Treatment Stage`))+
 # geom_bar(aes(y = AvClass,fill = type.y),stat = "identity",color="black",width=0.9)+
  geom_bar(data = AllSampClassprephelp,aes(y = AvClass,fill = type.y),stat = "identity",color="black",width=0.9)+
  geom_errorbar(data=Samples.Classtotals_stdev, colour="black",aes(ymin=AvtotClass,
                                                               ymax = AvtotClass+ StdevtotClass)) +
scale_fill_manual(values = mycolors)+
  theme_classic()+
  #theme(legend.justification=c(0,1), legend.position="right")+
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key.size = unit(0.5, 'lines'))+
  #facet_grid(rows = vars(Plant),switch = 'x')+
  theme(axis.title.y=element_text(size=18))+
  theme(legend.title = NULL) + 
  ylab("ARG copies / 16S rRNA \n")+
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=16,color="black"))+
  theme(axis.title.y=element_text(vjust=-1.5))+
  theme(axis.title.x=element_blank())+
  #scale_y_continuous(breaks=seq(0,2,0.5),limits=c(0,2))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=20,color="black"))
  #theme(axis.text.x=element_blank())

FLAllSampGenesplotrel 

pl<-list((VAAllSampGenesplotrel+theme(legend.position="none")),(FLAllSampGenesplotrel+theme(legend.position="none")))
pl2<-grid.arrange(grobs = pl,nrow=2)

ggsave("Fig3.pdf", plot = pl2, scale = 1, width = 14, height = 9, units = "in", dpi = 600) #save pdf in working directory folder


#Fig 4
#Change by class!
data<-read_excel("classchange081120_2.xlsx",sheet = 3)
data_long <- melt(data, id = c("Relative","Plant","Location","Treatment Stage"))
data_long$value <- as.numeric(data_long$value)

#specify the order you want the treatment plants to appear in
#data_long$Replace = factor(data_long$Replace, levels=c("IND-1","IND-2","HKG-1","HKG-2","PHL-1","PHL-2","USA-1","USA-2","CHE-1","CHE-2","SWE-1","SWE-2"))
#make a column that just indicates whether the value is positive or negative (i.e. increase or decrease) so that different colors can be assigned
data_long$change <- data_long$value/abs(data_long$value)
data_long$change <- str_replace(data_long$change, "-1", "decrease")
data_long$change <- str_replace(data_long$change, "1", "increase")
data_long$change <- str_replace(data_long$change, "NaN", "no change")

#create a magnitude column to specify the size of each bubble
data_long$magnitude <- abs(data_long$value)

data_long$variable<-str_replace_all(data_long$variable,"mls","MLS")
data_long$variable <-factor(data_long$variable,levels=c("fosfomycin",      "triclosan",       "phenicol"  ,"quinolone",  "trimethoprim",    "rifamycin" ,      "sulfonamide"   , 
                                                        "peptide"  ,       "beta-lactam"  ,   "glycopeptide"   , "aminoglycoside" , "tetracycline"  ,  "MLS"   ,"nucleoside", 
                                                        "other"   ,"multidrug"  ))

VA_classChange<- subset(data_long, data_long$Plant=="Ozone BAC/GAC" & data_long$Relative == "Influent")

VA_classChange$`Treatment Stage`<-factor(VA_classChange$`Treatment Stage`,levels<-c("Primary","Act. Sludge","Secondary",
                                                                                    "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination",
                                                                                    "3 days", "6 months", "Background"))

VA_classChangeplot<-ggplot(VA_classChange, aes(x = `Treatment Stage`, y = variable, label = value)) +
  geom_point(aes(size = magnitude, colour = change)) + 
  scale_size(range = c(1,13)) +
  theme_bw() +
  labs(
    y = "Antibiotic Resistance Class", 
    x = "Ozone BAC/GAC",
    size = "magnitude", 
    colour = "change during treatment") +
  scale_colour_manual(values = c("green3", "red", "grey")) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 60, hjust = 1))

FL_classChange<- subset(data_long, data_long$Plant=="Denitrification/Filtration Chlorination"& data_long$Relative == "Influent")

FL_classChange$`Treatment Stage`<-factor(FL_classChange$`Treatment Stage`,levels=c("Primary","Act. Sludge","Secondary",
                                                                                   "Denitrification","Chlorination","Cl Storage","Short Dist.",
                                                                                   "Long Dist."))

FL_classChangeplot<-ggplot(FL_classChange, aes(x = `Treatment Stage`, y = variable, label = value)) +
  geom_point(aes(size = magnitude, colour = change)) + 
  scale_size(range = c(1,13)) +
  theme_bw() +
  
  labs(
    x = "Denitrification/Filtration \n Chlorination",
    size = "Log Difference \n from Influent", 
    colour = "change during treatment") +
  theme(axis.title.y = element_blank())+
  scale_colour_manual(values = c("green3", "red", "grey")) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 60, hjust = 1))

VA_classChange<- subset(data_long, data_long$Plant=="Ozone BAC/GAC" & data_long$Relative == "Secondary")

VA_classChange$`Treatment Stage`<-factor(VA_classChange$`Treatment Stage`,levels<-c(
                                                                                    "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination",
                                                                                    "3 days","6 months","Background"))

VA_classChangeplot_WR<-ggplot(VA_classChange, aes(x = `Treatment Stage`, y = variable, label = value)) +
  geom_point(aes(size = magnitude, colour = change)) + 
  scale_size(range = c(1,13)) +
  theme_bw() +
  labs(
    y = "Antibiotic Resistance Class", 
   
    size = "magnitude", 
    colour = "change during treatment") +
  scale_colour_manual(values = c("green3", "red", "grey")) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 60, hjust = 1))

FL_classChange<- subset(data_long, data_long$Plant=="Denitrification/Filtration Chlorination"& data_long$Relative == "Secondary")

FL_classChange$`Treatment Stage`<-factor(FL_classChange$`Treatment Stage`,levels=c(
                                                                                   "Denitrification","Chlorination","Cl Storage","Short Dist.",
                                                                                   "Long Dist."))

FL_classChangeplot_WR<-ggplot(FL_classChange, aes(x = `Treatment Stage`, y = variable, label = value)) +
  geom_point(aes(size = magnitude, colour = change)) + 
  scale_size(range = c(1,13)) +
  theme_bw() +
  
  labs(size = "Log Difference \n from Secondary",  
    colour = "change during treatment") +
  theme(axis.title.y = element_blank())+
  scale_colour_manual(values = c("green3", "red", "grey")) +
  theme(text = element_text(size=20), 
        axis.text.x = element_text(angle = 60, hjust = 1))

#ggsave("flclasschange.pdf", plot = FL_classChangeplot, scale = 1, width = 14, height = 8, units = "in", dpi = 600) #save pdf in working directory folder
#+theme(legend.position="none") 
pl<-list((VA_classChangeplot+theme(legend.position="none") ),(FL_classChangeplot),
         (VA_classChangeplot_WR+theme(legend.position="none") ),(FL_classChangeplot_WR))
#lay<-rbind(c(1,1,1,1,1,1,1,1,1,1,1),
#c(2,2,2,2,2,2,2,2,2,NA,NA))
pl2<-grid.arrange(grobs = pl,ncol = 2,label_size = 12, align = 'h',axis="l",widths = c(1,1.6)) #,layout_matrix = lay)
ggsave("ClassChanges_4panel_test2.pdf", plot = pl2, scale = 1, width = 14, height = 12, units = "in", dpi = 600) #save pdf in working directory folder

sul1<-subset(Samples.genetotals, gene == "SUL1")
sul1_av<-  ddply(sul1,.(Location,Plant,`Treatment Stage`, type.y), summarise,
                 GCtot_samp=sum(gc_samp),
                 freq_class_samp = mean(avClassfreqrel_samp),
                 sd_sul1= sd(avClassfreqrel_samp))
sul1_av$`Treatment Stage`<-factor(sul1_av$`Treatment Stage`,levels=c("Influent","Primary","Act. Sludge","Secondary",
                                                               "Denitrification", "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination",
                                                               "3 days","Background","Cl Storage","Short Dist.", "Long Dist."))

ggplot(data=sul1_av, aes(x = `Treatment Stage`))+
  geom_bar(aes(y =freq_class_samp,fill = `Treatment Stage`),stat = "identity",color="black",width=0.9)+ #factor(Resistance_Genes,levels = lvls)
  geom_errorbar(colour="black",aes(ymin=freq_class_samp,
                                   ymax = freq_class_samp+  sd_sul1)) +
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  #facet_grid(rows = vars(`Sampling Trip`)  ,scales = 'free_x',space = 'free')+
  theme(legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key.size = unit(0.5, 'lines'))+
  theme(axis.title.y=element_text(size=18))+
  theme(legend.title = NULL) +
  #theme(legend.position = "bottom")+
  theme(axis.text.x=element_text(angle = 60, hjust = 1, size=16,color="black"))+
  ylab("Count of Clinically Relevant Genes \n")+
  theme(axis.title.y=element_text(vjust=-1.5))+
  theme(axis.title.x=element_text(vjust=-1.5))+
  #  scale_y_continuous(breaks=seq(0,1.5,0.5),limits=c(0,1.5))+
  theme(legend.text = element_text(size=22,color="black")) +
  theme(axis.text.y=element_text( size=22,color="black")) + 
  facet_grid(rows = vars(Plant), scales = "free")


ColonyCounts<- read_excel("CDC Sampling Data_Final.xlsx",sheet = 1)

ColonyCounts$`CFU/100mL`<-as.numeric(ColonyCounts$`CFU/100mL`)
#ColonyCounts$`CFU/100mL point`<-as.numeric(ColonyCounts$`CFU/100mL point`)


sub<-subset(ColonyCounts, (Plant ==  "Denitrification-filtration/chlorination"
                           & Organism == "E.coli" & 
                             Resistance == "Resistant"))

a<-subset(sub,( Treatment == "Long Dist." )  )    
b<-subset(sub,(Treatment == "Chlorination" | Treatment == "Cl Storage" )   )
wilcox.test(a$`CFU/100mL`,b$`CFU/100mL`)


Summary<-ddply(ColonyCounts, .(Plant,Organism,Resistance,Treatment, Condition), summarize,
               Average =mean(log10(`CFU/100mL`),na.rm=TRUE), 
               Stdev = sd(log10(`CFU/100mL`),na.rm=TRUE)/sqrt(3),
               Difference = (Average-Stdev),
               DetectionLim= min(DetectionLimit,na.rm=TRUE))


Summary$Treatment<-factor(Summary$Treatment,levels=c("Influent","Primary","Act. Sludge","Secondary",
                                                     "Denitrification", "Floc Sed","Ozonation","BAC/GAC","UV","Chlorination","Cl Storage","Short Dist.","Long Dist.",
                                                     "3 days","6 months", "Background"))



Ecolitot<-Summary %>% subset((Organism == "E.coli"))
Ecoli<-Summary %>% subset((Organism == "E.coli"))
#Ecoli2<-Summary %>% subset((Organism == "E.coli"& `Sampling Trip`=="2"))

Ecolitot$Plant<- factor( Ecolitot$Plant, levels<- c("Ozone/BAC/GAC",
                                                    "Denitrification-filtration/chlorination"))
Ecoli$Plant<- factor( Ecoli$Plant, levels<- c("Ozone/BAC/GAC",
                                              "Denitrification-filtration/chlorination"))

Ecoliplot<-ggplot(Ecolitot, aes(x = Treatment,linetype = Resistance, group = 1))+
  geom_line(data = subset(Ecolitot,Resistance == "Sensitive"),aes(y = (DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Ecolitot,Resistance == "Resistant"),aes(y = (DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Ecoli,Resistance == "Sensitive"),aes(y = (Average)),linetype = "solid",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Ecoli,Resistance == "Resistant"),aes(y = (Average)),linetype = "solid",color = "black",na.rm=TRUE,size = 1)+
  #geom_line(data = Ecoli,aes(y = log10(DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_errorbar(data = Ecoli,linetype = "solid" ,width = 0.2, color ="black", mapping=aes( ymin=( Average-Stdev ),ymax=(Average +  Stdev)),
                size=1,
                position="identity")+
  geom_point(data = subset(Ecoli,Resistance == "Sensitive"),aes(y = (Average)),na.rm=TRUE,color = "#0b3d91",shape =21,fill ="#0b3d91",size = 5)+
  geom_point(data = subset(Ecoli,Resistance == "Resistant"),aes(y = (Average)),na.rm=TRUE,color = "#2bd6fb",shape =24,fill = "#2bd6fb",size = 5)+
  
  # geom_line(data = subset(Ecoli2,Resistance == "Sensitive"),aes(y = (Average)),na.rm=TRUE,color = "black",size = 1)+
  #  geom_line(data = subset(Ecoli2,Resistance == "Resistant"),aes(y = (Average)),na.rm=TRUE,color = "black",size = 1)+
  # geom_point(data = subset(Ecoli2,Resistance == "Sensitive"),aes(y = (Average)),na.rm=TRUE,color = "#0b3d91",shape = 21,fill =NA,size = 5)+
  #geom_point(data = subset(Ecoli2,Resistance == "Resistant"),aes(y = (Average)),na.rm=TRUE,color = "#2bd6fb",shape = 21,size = 5, fill = "#2bd6fb")+
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"))+
  theme(axis.title.x=element_blank(),
        strip.text = element_text(size=14))+
  labs(x=NULL,y="log(CFU/100mL)",size = 12)+
  scale_y_continuous(breaks=seq(-2,9,1),limits=c(-2,9))+
  facet_grid(cols=vars(Plant), scales= 'free')
#labs(x=NULL,y="CFU/100mL")
Ecoliplot
Efaectot<-Summary %>% subset((Organism == "E.faecium"))
Efaectot$Plant<- factor( Efaectot$Plant, levels<- c("Ozone/BAC/GAC",
                                                    "Denitrification-filtration/chlorination"))
Ecoli$Plant<- factor( Ecoli$Plant, levels<- c("Ozone/BAC/GAC",
                                              "Denitrification-filtration/chlorination"))

Efaecplot<-ggplot(Efaectot, aes(x = Treatment,linetype = Resistance, group = 1))+
  geom_line(data = subset(Efaectot,Resistance == "Sensitive"),aes(y = (DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Efaectot,Resistance == "Resistant"),aes(y = (DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Efaectot,Resistance == "Sensitive"),aes(y = (Average)),linetype = "solid",color = "black",na.rm=TRUE,size = 1)+
  geom_line(data = subset(Efaectot,Resistance == "Resistant"),aes(y = (Average)),linetype = "solid",color = "black",na.rm=TRUE,size = 1)+
  #geom_line(aes(y = log10(DetectionLim)),linetype = "dashed",color = "black",na.rm=TRUE,size = 1)+
  geom_errorbar(linetype = "solid" ,width = 0.2, color ="black",mapping=aes( ymin=(Average - Stdev),ymax=(Average +  Stdev)),
                size=1,
                position="identity")+
  geom_point(data = subset(Efaectot,Resistance == "Sensitive"),aes(y = (Average)),na.rm=TRUE,color = "#0b3d91",shape =21,fill ="#0b3d91",size = 5)+
  geom_point(data = subset(Efaectot,Resistance == "Resistant"),aes(y = (Average)),na.rm=TRUE,color = "#2bd6fb",shape =24,fill = "#2bd6fb",size = 5)+
  
  #  geom_line(data = subset(Efaec2,Resistance == "Sensitive"),aes(y = (Average)),color = "black",na.rm=TRUE,size = 1)+
  # geom_line(data = subset(Efaec2,Resistance == "Resistant"),aes(y = (Average)),color = "black",na.rm=TRUE,size = 1)+
  #geom_point(data = subset(Efaec2,Resistance == "Sensitive"),aes(y = (Average)),color = "#0b3d91",shape = 21,fill =NA,size = 5)+
  #geom_point(data = subset(Efaec2,Resistance == "Resistant"),aes(y = (Average)),color = "#2bd6fb",shape = 21,size = 5, fill = "#2bd6fb")+
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=12,color="black"),
        axis.text.y=element_text(size=12,color="black"))+
  theme(axis.title.x=element_blank(),
        strip.text = element_text(size=14))+
  labs(x=NULL,y="log(CFU/100mL)",size = 12)+
  scale_y_continuous(breaks=seq(-2,9,1),limits=c(-2,9))+
  facet_grid(cols=vars(Plant), scales= 'free')


#labs(x=NULL,y="CFU/100mL")
Efaecplot

ggsave("Fig5.pdf", Ecoliplot, Efaecplot)

RawData<- read_excel("Culture.xlsx")

#RawData<-subset(RawData, (Organism != "E. coli" & Organism != "Enterococcus spp." ))

long<-reshape2::melt(RawData, id = c("Location", "Sampling Event", "Organism"))
long<-subset(long, (Organism != "Enterococcus spp." & Organism != "E. coli")&
               #                        (`Sampling Event` =="S3" | `Sampling Event` =="S4" ) &
               (variable == "# Isolates Screened" | variable == "# Isolates Confirmed to Genus/Species" ))
long$value<-as.numeric(as.character(long$value))

long<-long %>% dplyr::group_by(Location,Organism, variable) %>% dplyr::summarise(sumTot =sum(value, na.rm=TRUE))

#long<-long %>% dplyr::group_by(Location,Organism, variable) %>% dplyr::summarise(sumTot =mean(value, na.rm=TRUE))



mycolors  <- c( "mediumblue",  
                "indianred2", "darkmagenta", "green", "pink", "darkorchid1", 
                "dark red", "yellow", "gray47", "dark green",  "darkorange1", "purple", "darkturquoise", 
                "gray67", "goldenrod2", "magenta2")
long$Location<- factor( long$Location, levels<- c("Ozone/BAC/GAC",
                                                  "Denitrification-filtration/chlorination"))
BarChart<-ggplot(long, aes(x = Organism, y =sumTot, fill = variable)) + 
  geom_bar(stat = "identity",color="black", position = "dodge",width=0.9)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  labs(y="Number of Isolates",size = 12)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=12,color="black"),
        #axis.title.y = element_blank(),
        legend.text = element_text(size = 12),
        legend.title=element_blank(),
        strip.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        axis.text.x=element_text(angle = 60, hjust = 1,size=12,color="black", face = "italic"))+
  theme(legend.position = "bottom")+ 
  
  #scale_y_continuous(breaks=seq(0,1200,200),limits=c(0,1200))+
  theme(strip.background = element_blank())+
  facet_grid(cols=vars( Location ),space = 'free') 
BarChart
#ggsave("FL_3_4_Bars.pdf", plot = BarChart, scale = 1, width =8 , height =11, units = "in", dpi = 600)

pl<-list((Efaecplot+theme(legend.position="none")),(Ecoliplot+theme(legend.position="none")) )
pl2<-grid.arrange(grobs = pl,ncol=1,nrow=3, label_size = 12, align = 'hv',axis="l", hjust=-2, vjust=1)
ggsave("Fig5.pdf", plot = pl2, scale = 1, width =8 , height =11, units = "in", dpi = 600) #save pdf in working directory folder

