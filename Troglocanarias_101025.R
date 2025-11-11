##############
# Functional diversity of the subterranean aquatic fauna of the Canary Islands 
# Script written by Alejandro Martinez
# First version: 09.06.2021
# Last version: 07.11.2025
###############

###### Index of R-Packages -----------------------------------------------------------

# Visualization of spatial objects
library(sf)
library(sp)
library(lme4)
library(raster) # Combine rasters and extract climatic variables from the raster

# Ecological analyses and regression
library(BAT) # Main package
library(performance)
library(modEvA)

# Elaboration of figures
library(ggplot2) 
library(dplyr) # Used to arrange the data
library(ggrepel)



### Working directory ---------------------------------------------------------------

setwd("/Users/amartinez/Dropbox/_Papers/_READY/submitted - Martinez - Canarias Troglobites/01 Manuscript/06 R1_Journal of Biogeography/00 analyses")

# Functions

source("functions/abspres.R")



####### PART 1: DATA PREPARATION ------------------------------------------------


# ####### Map and check the distribution of localities ----------------------------
# 
#     access <- read.csv2("data/localities_211004.csv", dec=".", sep=";") #
#     access <- access[access$included=="yes",]
#     head(access)
#     str(access)
# 
# # Georreference the dataset
#     
#     coordinates(access) <- c("Long.Y", "Lat.X")
#     crs.geo <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#     proj4string(access) <- crs.geo
# 
# # Plot interactive map
#     
#     mapview::mapView(access, legend=T,alpha=0.2,cex=4,zcol = 'type',
#                    popup = leafpop::popupTable(access),
#                   map.types=c("Esri.WorldImagery","OpenStreetMap.DE"))
#     rm(crs.geo)
#     
#       # Everything is correctly placed! Next!
# 
#     
####### Extract variables for each locality from the raster maps -----------------
    
#### Load the raster files
#
#  maps <- list.files(path = "maps/", pattern = "\\.asc$",
#                     all.files = TRUE, full.names = TRUE) 
#
#### Stack the raster into one map
#    
#     rasStack <- stack(maps)
#
#
#### Extract the information after the localities
#    
#    rasValue <- raster::extract(rasStack, access)
#    data.geo.localities <- cbind(access,rasValue); head(data.locality)
#
#    
#    write.table(data.geo.localities, "data/Dataset_localities.csv",sep=";",dec=".")
#
#    rm(maps,rasStack,rasValue, data.geo.localities)
#
#
# data.locality <- read.csv2("data/Dataset_localities.csv")
#    
#

####### Prepare the functional trait matrix ----------------------------------------

# Load raw trait dabases
    
    traits <- read.csv2("data/traits_211004.csv", dec=".", sep=";") 
    head(traits)

# Filter to keep the useful traits and set rownames
    
    traits02 <- traits[c("ID","L.max01","E.max01","eyes","Ant.max01","Leg.max01","pigmentation.or",
                         "trophic.niche", "ecology")] 
    
    psych::pairs.panels(traits02)
    
    row.names(traits02) <- traits02[,1]; traits02 <- traits02[,-1] # row.names
    head(traits02); rm(traits)
    

# Scale the traits
    
    par(mfrow=c(2,3))
    L.max01 <- scale(traits02$L.max01); hist(L.max01)
    E.max01 <- scale(traits02$E.max01); hist(E.max01)
    Leg.max01 <- scale(traits02$Leg.max01); hist(Leg.max01)
    Leg.max01 <- scale(traits02$Leg.max01); hist(Leg.max01)
    Ant.max01 <- scale(traits02$Ant.max01); hist(Ant.max01)
    dev.off()
    

# Build a new scaled trait matrix
    
    traits03 <- cbind(L.max01,E.max01,Ant.max01,Leg.max01,traits02$pigmentation.or,
                      traits02$trophic.niche, traits02$ecology)
    row.names(traits03) <- row.names(traits02)
    
    colnames(traits03) <- c("L.max01","E.max01","Ant.max01","Leg.max01",
                            "pigmentation.or","trophic.niche","ecology")

# Filter out species without occurrence data or incomplete traits
    
    traits04 <- as.data.frame(traits03[!row.names(traits03) %in% "Ar031", ]) #Ar031 has a lot of missing data
    #nooccurrence <- c("Co008","Co009","Co033")
    #traits04 <- as.data.frame(traits04[!row.names(traits04) %in% nooccurrence, ])
    
    #rm(nooccurrence)
    rm(L.max01,E.max01,Ant.max01,Leg.max01, traits02, traits03) 

    
    
# Calculation of gower distance matrix
    
    gower.mat <- BAT::gower(traits04)
    euc.pco <- labdsv::pco(gower.mat,k=4)
   
    SUM <- sum(euc.pco$eig[euc.pco$eig > 0])

    barplot(euc.pco$eig)
    (euc.pco$eig[1]+euc.pco$eig[2]+euc.pco$eig[3]+euc.pco$eig[4])/sum(euc.pco$eig[euc.pco$eig > 0]) 
    # The 0.89 of trait variance is explained by first 4 axes

# Produce the file for the analyses: traits05
    
    traits05 <- euc.pco$points # <<-- Trait matrix based on Gower distance
    colnames(traits05) <- c(paste(rep("axes_",4),1:4,sep=''))
  

    
####### Prepare a presence/absence community matrices -------------------------------

records <- read.csv2("data/records_simple_211004.csv") # Table with records
head(records)


### Summary statistics per cave

mean(table(records02$Loc.ID))
range(table(records02$Loc.ID))
sd(table(records02$Loc.ID))/sqrt(length(records02$Loc.ID))

### Summary statistics per island

mean(table(records02$Island))
range(table(records02$Island))
sd(table(records02$Island))/sqrt(length(records02$Island))

# Select the species with traits only
    
    names <- row.names(traits04) 
    records01 <- records[records$ID %in% names, ]
    records02 <- records01[records01$Included =="yes", ]
    rm(records,records01)

    
# Matrix 1: grouped by islands (e.g.: F,C,T)
    comm.isl <- abspres(data=records02, sites.col="Island", sp.col = "ID", keep.n=F)
    row.names(comm.isl) <- comm.isl[,1]
    comm.isl <- comm.isl[,-1]
    
# Ordering traits as comm before calculation
    traits05 <- traits05[match(colnames(comm.isl),rownames(traits05)),]
    


    
######## Hypothesis 1: Hypervolumes converge across islands  -------------------------


#### Calculating of functional hypervolumes and metrics -------------------------------        
      
  
      kernelFD.isl <- BAT::kernel.build(comm=comm.isl,trait=traits05,abund=FALSE,cores=6,
                                   method="gaussian")
  
        rich.isl <- kernel.alpha(kernelFD.isl)
        even.isl <- kernel.evenness(comm=kernelFD.isl)
        spec.isl <- as.vector(rowSums(comm.isl))
        names.isl <- c("C","F","G","H","P","T")
        
        Alpha.isl <- data.frame(island=names.isl,Species=spec.isl,
                                Richness=rich.isl,Even=even.isl)
        
          rm(rich.isl,even.isl,spec.isl,names.isl)
      
      
      # saveRDS(kernelFD.isl, file = "kernelFD.rds")

####### Plot hypervolumes for Figure 01 to show hv converge --------------------
    
library(RColorBrewer)
    
    names.isl <- c("C","F","G","H","P","T")

    # Remember the colour code of each island:
    ## F - 993300 
    ## T - 669900
    ## C - 66CC00
    ## G - 739E49
    ## P - 3399CC
    ## H - 66CCFF
    
islands.color.hv6 <- c("#66CC00","#993300","#739E49",
                           "#66CCFF","#3399CC","#669900")
   
library(alphahull);library(rgl)
    
    plot(kernelFD.isl,
         show.3d=FALSE,
         show.random=TRUE,
         show.density=TRUE,
         show.data=FALSE,
         show.centroid=TRUE,
         cex.centroid=1.5,
         centroid.alpha=0.5,
         point.alpha.min=0.2, 
         point.dark.factor=0.7,
         cex.random=0.3,
         cex.data=0.5,
         contour.lwd=0.5,
         num.points.max.random = 1000,
         col= islands.color.hv6,
         names=c(  paste("Axes 1\n[",round(((euc.pco$eig[1])/SUM),2),"%]",sep=''),
                   paste("Axes 2\n[",round(((euc.pco$eig[2])/SUM),2),"%]",sep=''),
                   paste("Axes 3\n[",round(((euc.pco$eig[3])/SUM),2),"%]",sep=''),
                   paste("Axes 4\n[",round(((euc.pco$eig[4])/SUM),2),"%]",sep='')  ))
    
    
####### Functional beta diversity: convergence relies on difference species -------------------------

    
# Beta functional diversity per island (F,C,T)  
    betaFD.isl <- kernel.beta(comm=kernelFD.isl)
    betaTD.isl <- BAT::beta(comm.isl, abund=F)
    

### Organize the data  into a matrix for the plot (this code needs cleaning)
   
# Functional metrics

    betaFDm <- data.frame(Total = as.vector(betaFD.isl$Btotal),
                          Repl = as.vector(betaFD.isl$Brepl),
                          Rich = as.vector(betaFD.isl$Brich),
                          cat = rep("functional",15))


# Taxonomic metrics
    
    betaTDm <-data.frame(Total = as.vector(betaTD.isl$Btotal),
                         Repl = as.vector(betaTD.isl$Brepl),
                         Rich = as.vector(betaTD.isl$Brich),
                         cat = rep("taxonomic",15))

    
#Append both metrics into a single data file
    
    beta <- rbind(betaFDm,betaTDm); rm(betaFDm,betaTDm)
    

# Beta diversity density plots
    
    pb1 <- (ggplot(data=beta, aes(x=Total, group=cat, fill=cat)) + 
              geom_density(alpha=.4) + 
              xlab("βtotal (βreplacement +  βrichness)") +
              ylab("Density")+
              scale_fill_manual(values = c("#4682B4", "#B47846"))+
              labs(title="a")+
              theme_bw()+
              theme(legend.position = c(0.3, 0.85), 
                    legend.box.background = element_rect(colour = "black"),
                    plot.title = element_text(face="bold")))
    
    pb2 <- (ggplot(data=beta, aes(x=Repl, group=cat, fill=cat)) + 
              geom_density(alpha=.4) + 
              xlab("βreplacement") +
              ylab("Density")+
              scale_fill_manual(values = c("#4682B4", "#B47846"))+
              labs(title="b")+
              theme_bw()+
              theme(legend.position = "none"))
    
    pb3 <- (ggplot(data=beta, aes(x=Rich, group=cat, fill=cat)) + 
              geom_density(alpha=.4) + 
              xlab("βrichness") +
              ylab("Density")+
              scale_fill_manual(values = c("#4682B4", "#B47846"))+
              labs(title="c")+
              theme_bw()+
              theme(legend.position = "none"))
    
    gridExtra::grid.arrange(pb1,pb2,pb3, ncol = 2, 
                            layout_matrix = cbind(c(1,1), c(2,3)))
    rm(pb1,pb2,pb3)
    
    rm(betaFD.isl,betaTD.isl)

   
            
###### Hypothesis 2: Analyses at species level contribution ------
 
  #contri.isl <- kernel.contribution(kernelFD.isl)
  #write.csv2(contri.isl,"data/contribution_island_jan_211004.csv")
      
# Load the species contribution matrix    
    
    contri.isl <- read.csv2("data/contribution_island_jan_211004.csv", dec=",", sep=";")
    row.names(contri.isl) <- contri.isl[,1];contri.isl <- contri.isl[,-1]
    as.data.frame(contri.isl)
    str(contri.isl)

# Distribution species contribution per island
    
    contri.isl.t <- as.data.frame(t(contri.isl))
    colnames(contri.isl.t)[1]
    contri.isl01 <- c(contri.isl.t)
      
    
    contr.vector <- data.frame()
          for (i in 1:ncol(contri.isl.t)){
            this.col <- as.numeric(contri.isl.t[,i])
            this.col.name <- rep((colnames(contri.isl.t)[i]),nrow(contri.isl.t))
            this.island <- data.frame(row.names(contri.isl.t),this.col.name,this.col)
            contr.vector <- rbind(contr.vector,this.island)}
    
    colnames(contr.vector) <- c("ID","Island","contribution")
    
    contr.vector$Island = factor(contr.vector$Island)
    contr.vector01 <- contr.vector[which (contr.vector$contribution > 0),]

    str(contr.vector); rm(contri.isl,contri.isl.t,contri.isl01)
    rm(this.col,this.col.name,names.isl,names)
    
## Table summary     
   
     std <- function(x) sd(x)/sqrt(length(x))

    summary <- contr.vector01 %>%
      group_by(Island) %>%
      summarise(
        Contribution_mean = mean(contribution)*1000,
        Contribution_se = std(contribution)*1000,
        Contribution_median = median(contribution)*1000,
        mad_contrib_rel    = mad(contribution)*1000,
        S                  = n()
      ) %>%
      ungroup()
    
    summary <- as.data.frame(summary)

    
    Alpha.isl <- merge(Alpha.isl, summary, by.x = "island", by.y = "Island")
 
    colnames(Alpha.isl)[1] <- "Island"
 
#h Colour palette
      islands.color.hv6 <- c("#66CC00","#993300","#739E49",
                               "#66CCFF","#3399CC","#669900")

      

# Plots                
    cb1 <- (ggplot(data=contr.vector[contr.vector$contribution > 0 & 
                                    contr.vector$Island == "P"|
                                    contr.vector$Island =="H",],
                      aes(x=contribution, group=Island, fill=Island)) + 
            geom_density(alpha=.2) + 
            xlab("Species contribution") +
            xlim(0,0.01) +
            ylab("Density")+
            ylim(0,1500)+
            scale_fill_manual(values = c("#66CCFF","#3399CC")) +
            labs(title="a")+
            theme_bw()+
            theme(legend.position = c(0.8, 0.7), 
                  legend.box.background = element_rect(colour = "black"),
                  plot.title = element_text(face="bold")))
    
    cb2 <- (ggplot(data=contr.vector[contr.vector$contribution > 0 & 
                                       contr.vector$Island == "T"|
                                       contr.vector$Island =="G" |
                                       contr.vector$Island == "C",],
                   aes(x=contribution, group=Island, fill=Island)) + 
              geom_density(alpha=.2) + 
              xlab("Species contribution") +
              xlim(0,0.01) +
              ylab("Density")+
              ylim(0,1500)+
              scale_fill_manual(values = c("#66CC00","#739E49", "#669900")) +
              labs(title="a")+
              theme_bw()+
              theme(legend.position = c(0.8, 0.7), 
                    legend.box.background = element_rect(colour = "black"),
                    plot.title = element_text(face="bold")))
    
    cb3 <- (ggplot(data=contr.vector[contr.vector$contribution > 0 & 
                                       contr.vector$Island == "F",],
                   aes(x=contribution, group=Island, fill=Island)) + 
              geom_density(alpha=.2) + 
              xlab("Species contribution") +
              xlim(0,0.01) +
              ylab("Density")+
              ylim(0,1500)+
              scale_fill_manual(values = c("#993300")) +
              labs(title="a")+
              theme_bw()+
              theme(legend.position = c(0.8, 0.7), 
                    legend.box.background = element_rect(colour = "black"),
                    plot.title = element_text(face="bold")))
      
    
gridExtra::grid.arrange(cb1, cb2, cb3, nrow = 3)
    
rm(cb1,cb2,cb3)  
rm(contr.vector,contr.vector01,summary, this.island) 


data.isl <- read.csv2("data/islands.csv", dec=".",sep=";")
data.isl1 <- merge(Alpha.isl, data.isl, by.x = "Island", by.y = "island",  all=F)  

data.isl1[is.na(data.isl1)] <- 0 ## <<-- Island richness dataset

m_size <- lm(Contribution_mean ~ Species, data = data.isl1)
summary(m_size)
data.isl1$resid_size <- resid(m_size)
data.isl1$AgeClass <- as.factor(c("Mature","Senescent","Mature","Young","Young","Mature"))
performance::check_model(m_size)


# Island colours (using your code)
cols_island <- c(
  "F" = "#993300",   # Fuerteventura
  "T" = "#669900",   # Tenerife
  "C" = "#66CC00",   # Gran Canaria
  "G" = "#739E49",   # La Gomera
  "P" = "#3399CC",   # La Palma
  "H" = "#66CCFF"    # El Hierro
)

ggplot(data.isl1, aes(x = AgeClass, y = resid_size, colour = Island)) +
  geom_point(size = 3,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_text_repel(aes(label = Island),
                  position = position_jitter(width = 0.05, height = 0),
                  show.legend = FALSE,
                  size = 3) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 3,
               fill = "white",
               colour = "black") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = cols_island, name = "Island") +
  theme_bw()+
  labs(x = "Island age group",
       y = "Mean contribution (residuals after richness)")
aov_out <- anova(lm(resid_size ~ AgeClass, data = data.isl1))
aov_out


set.seed(123)
perm_F <- replicate(10000, {
  perm_lab <- sample(data.isl1$AgeClass)
  anova(lm(data.isl1$resid_size ~ perm_lab))$`F value`[1]
})
obs_F <- aov_out$`F value`[1]
perm_p <- mean(perm_F >= obs_F)
perm_p


# Optional planned contrasts vs Mature (directional expectation)
data.isl1 %>%
  group_by(AgeClass) %>%
  summarise(mean_resid = mean(resid_size),
            sd_resid = sd(resid_size),
            n = dplyr::n())




###### Hypothesis 3a: Mantel test relationships --------------------------------------------

# Visualization of spatial objects
library(sf)
library(sp)
library(lme4)
library(raster) # Combine rasters and extract climatic variables from the raster

# Ecological analyses and regression
library(BAT) # Main package
library(performance)
library(modEvA)
library(relaimpo)

# Taxonomic tree manipulations
library(ape)

# Elaboration of figures
library(ggplot2) 
library(dplyr) # Used to arrange the data



##### Preparation of the initial matrices ---------------------------------------------- 

#records <- read.csv2("data/records_simple_211004.csv") # Table with records
#str(records)

# Matrix occurrences grouped by localities (eg. Cueva del Viento, C. el Llano)
    comm.loc <- abspres(data=records02, sites.col="Loc.ID", sp.col = "ID", keep.n=F)
    row.names(comm.loc) <- comm.loc[,1]
    comm.loc <- comm.loc[,-1]
    

# Add climatic variables (all correlated, we kept only altitude, code has been simplified)
    
data.raster <- read.csv2("data/Dataset_localities.csv", sep=";", dec=",")
          
          data.raster <- data.raster[data.raster$Loc.ID %in% row.names(comm.loc), ] 

         # data.raster.cont <- data.raster0[,c(14:21)]
         # psych::pairs.panels(data.raster.cont) # They are all correlated with altitude
          
          data.niche <- data.raster[,c("Loc.ID","altitude","Long.Y","Lat.X")]
          
          data.niche <- merge(records02,data.niche, by="Loc.ID")
          
          species_niche <- data.frame(species = as.factor(data.niche$ID),
                                 altitude = as.numeric(data.niche$altitude),
                                 longitude = as.numeric(data.niche$Long.Y),
                                 latitude = as.numeric(data.niche$Lat.X),
                                 Island = as.factor(data.niche$Island))
          
          species_island <- unique(data.frame(species = as.factor(data.niche$ID),
                                         island = data.niche$Island))
          
          str(species_niche)

rm(data.raster,data.niche)
      
      
##### Calculation of the centroid and altitudinal mean ----------------------------

centroid <- species_niche %>%
  group_by(species, Island) %>%
  dplyr::summarise(x = mean(longitude),
                   y = mean(latitude),
                   z = mean(altitude)) #limit area of a station is defined by its furthest record

centroid <- as.data.frame(centroid)
rownames(centroid) <- paste0(centroid$species,".",centroid$Island)
centroid$Island <- as.factor(centroid$Island)



#### Calculation of the geographical distance matrix  --------------------------------

d.geographic.per.isl <- list()

  for (i in 1:nlevels(centroid$Island)){
    
    centroid.i <- centroid[ which (centroid$Island == levels(centroid$Island)[i]),]
    centroid.i$Island <- droplevels(centroid.i$Island)
    
      d.geographic <- geodist::geodist(centroid.i, measure = "geodesic")
      
            colnames(d.geographic)<- rownames(centroid.i)
            rownames(d.geographic)<- rownames(centroid.i)
      
            d.geographic <- d.geographic[order(rownames(d.geographic)), ]
            d.geographic <- t(d.geographic)
            d.geographic <- d.geographic[order(rownames(d.geographic)), ]
            
               # d.geographic[1,] = NA
               # for(j in 2:ncol(d.geographic)){
               #   d.geographic[j,j:ncol(d.geographic)] = NA
               # }
            
            d.geographic.per.isl[[i]]  <- d.geographic 
            
           names(d.geographic.per.isl)[[i]] <- levels(centroid$Island)[i]
  }     

 
#rm(centroid,d.geographic, XX, YY) #cleaning
rm(centroid.i,d.geographic) 
  


##### Calculation of habitat distances (=altitude distance) ----------------------------------------------                  


  d.altitudinal.per.isl <- list()

    for (i in 1:nlevels(centroid$Island)){
      
      centroid.i <- centroid[ which (centroid$Island == levels(centroid$Island)[i]),]
      centroid.i$Island <- droplevels(centroid.i$Island)      

      d.altitudinal <- as.matrix(dist(centroid.i$z, method = "euclidean"))
      rownames(d.altitudinal) <- rownames(centroid.i)
      colnames(d.altitudinal) <- rownames(centroid.i)
      d.altitudinal <- d.altitudinal[order(rownames(t(d.altitudinal[order(rownames(d.altitudinal)), ]))), ]
            
            d.altitudinal <- d.altitudinal[order(rownames(d.altitudinal)), ]
            d.altitudinal <- t(d.altitudinal)
            d.altitudinal <- d.altitudinal[order(rownames(d.altitudinal)), ]
            
           # d.altitudinal[1,] = NA
           # for(j in 2:ncol(d.altitudinal)){
           #   d.altitudinal[j,j:ncol(d.altitudinal)] = NA
           # }  
            
            d.altitudinal.per.isl[[i]] <- d.altitudinal
            names(d.altitudinal.per.isl)[[i]] <- levels(centroid$Island)[i]
    }      
 
rm(d.altitudinal)

##### Calculation of functional distances ---------------------------------------------- 

# Arrange Gower matrix

d.functional <- as.matrix(gower.mat)
d.functional <- d.functional[order(rownames(d.functional)), ]
d.functional <- t(d.functional)
d.functional <- d.functional[order(rownames(d.functional)), ]

d.functional.per.isl <- list()

for (i in 1:nlevels(centroid$Island)){
  centroid.i <- centroid[ which (centroid$Island == levels(centroid$Island)[i]),]
  centroid.i$Island <- droplevels(centroid.i$Island)      
  
    species_island.T <- species_island[ which (species_island$island == levels(centroid.i$Island)),]
        this.dfun <- d.functional[ which (colnames(d.functional) %in% species_island.T$species),]
        this.dfun <- this.dfun[, which ( colnames(this.dfun) %in% species_island.T$species) ]

    d.functional.per.isl[[i]] <- this.dfun
        names(d.functional.per.isl)[[i]] <- levels(centroid$Island)[i]
}

rm(d.functional)

### Mantel tests

mantel <- data.frame()
islands <- c("C","F","G","H","P", "T")
mantel.perm <- 999

for (i in 1:length(islands)){
  print(paste0("Mantel test ",i," out of " ,length(islands)))
  mantel.alt <- vegan::mantel(d.functional.per.isl[[i]], d.altitudinal.per.isl[[i]], 
                method = "spearman", permutations = mantel.perm, na.rm = TRUE)
  
  mantel.geo <- vegan::mantel(d.functional.per.isl[[i]], d.geographic.per.isl[[i]], 
                method = "spearman", permutations = mantel.perm, na.rm = TRUE)  
  
  
  mantel.geoxalt <- vegan::mantel.partial(d.functional.per.isl[[i]], d.geographic.per.isl[[i]], d.altitudinal.per.isl[[i]],
                        method = "pearson", permutations = mantel.perm, 
                        strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores")) 
  
  mantel.altxgeo <- vegan::mantel.partial(d.functional.per.isl[[i]], d.altitudinal.per.isl[[i]], d.geographic.per.isl[[i]], 
                        method = "pearson", permutations = mantel.perm, 
                        strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores")) 
  
  mantel.this.island <- data.frame(island = rep(islands[i],4),
                            test = c("altitude","geography","geoXalt","altXgeo"),
                            MantelR = c(mantel.alt$statistic,
                                        mantel.geo$statistic, 
                                        mantel.geoxalt$statistic,
                                        mantel.altxgeo$statistic),
                            pvalue = c(mantel.alt$signif,
                                       mantel.geo$signif,
                                       mantel.geoxalt$signif,
                                       mantel.altxgeo$signif))
  
  mantel <- rbind(mantel, mantel.this.island)
  rm(mantel.this.island,mantel.alt,mantel.geo,mantel.altxgeo,mantel.geoxalt)
}


# mantel.alt <- mantel[which (mantel$test == "alr"),c(3,4)]
# colnames(mantel.alt) <- c("Mantel.alt.R","Mantel.alt.p")

Alpha.isl2 <- data.frame(Alpha.isl$Island,
                         Alpha.isl$Species,
                         Alpha.isl$Richness,
                         Alpha.isl$Even,
                         abs(mantel[which (mantel$test == "altitude"),c(3,4)])[1],
                         abs(mantel[which (mantel$test == "altitude"),c(3,4)])[2],
                         abs(mantel[which (mantel$test == "geography"),c(3,4)])[1],
                         abs(mantel[which (mantel$test == "geography"),c(3,4)])[2],
                         abs(mantel[which (mantel$test == "altXgeo"),c(3,4)])[1],
                         abs(mantel[which (mantel$test == "altXgeo"),c(3,4)])[2],
                         abs(mantel[which (mantel$test == "geoXalt"),c(3,4)])[1],
                         abs(mantel[which (mantel$test == "geoXalt"),c(3,4)])[2])
  

colnames(Alpha.isl2) <- c("island",
                          "Species",
                          "Richness",
                          "Even",
                          "MantelR.alt",
                          "Mantelp.alt",
                          "MantelR.geo",
                          "Mantelp.geo",
                          "MantelR.altXgeo",
                          "Mantelp.altXgeo",
                          "MantelR.geoXalt",
                          "Mantelp.geoXalt")



data.isl <- read.csv2("data/islands.csv", dec=".",sep=";")
data.isl <- merge(Alpha.isl2,data.isl, by="island",all=F)  

data.isl[is.na(data.isl)] <- 0 ## <<-- Island richness dataset

data.isl[,c(7:10)] <- abs(data.isl[,c(7:10)])

str(data.isl)



m1 <- plot(data.isl$MantelR.geoXalt ~ data.isl$Age)
m2 <- plot(data.isl$MantelR.geoXalt ~ data.isl$Distance)
m3 <- plot(data.isl$MantelR.altXgeo ~ data.isl$Age)
m4 <- plot(data.isl$MantelR.altXgeo ~ data.isl$Distance)


gridExtra::grid.arrange(m1, m2, m3, m4, nrow = 2)

mantel <- data.isl[c(1:15,21)]

colnames(mantel)[16] <- "Distance"


islands.color.hv6 <- c("#66CC00","#993300","#739E49",
                       "#66CCFF","#3399CC","#669900")
    



data.isl <- data.isl[ which (data.isl$island != "L"),]

# Plotting the R square
plot(data.isl$Distance, abs(data.isl$MantelR.altXgeo), pch=24)
points(data.isl$Distance, abs(data.isl$MantelR.geoXalt), pch=19)



# Generating the curves
mantel$MantelR.altXgeo2 <- (mantel$MantelR.altXgeo)^2
mantel$Distance2 <- (mantel$Distance)^2
mantel$Distance3 <- (mantel$Distance)^3

qR.eco <- lm(MantelR.altXgeo ~ Distance + Distance2, data = mantel)
qR.geo <- lm(MantelR.geoXalt ~ Distance + Distance2, data = mantel)

summary(qR.eco)
summary(qR.geo)


dist.values <- seq(0,600,1)
eco.predict <- predict(qR.eco, list(Distance=dist.values, Distance2=dist.values^2))
geo.predict <- predict(qR.geo, list(Distance=dist.values, Distance2=dist.values^2))
lines(dist.values, eco.predict, col="blue")
lines(dist.values, geo.predict, col="red")



###### Hypothesis 3b: Species co-occurrence --------------------------------------------


coocurrence.per.isl <- list()

for (i in 1:6) {
  matrix.i <- records02[ which (records02$Island == islands[i]),]
  co.rec.i <- crossprod(table(matrix.i[c("Loc.ID","ID")]))
  diag(co.rec.i) <- 0
  
  coocurrence.per.isl[[i]] <- co.rec.i
  names(coocurrence.per.isl)[[i]] <- islands[i]
}

# Calculate co-occurences per island


for (k in 1:6){
  matrix.k <- coocurrence.per.isl[[k]]
        for (i in 1:ncol(matrix.k)) {
              XX = c(); SPECIES = c()
              YY=cbind(rep(rownames(matrix.k)[i],ncol(matrix.k)), rownames(matrix.k))
              XX=rbind(XX,YY)
              SPECIES=c(SPECIES,rep(rownames(matrix.k)[i],ncol(matrix.k)))
              SPECIES=c(SPECIES,rep(islands[k],ncol(matrix.k)))
            }
  DB.corec <- data.frame(ID=paste0(XX[,1],XX[,2]),XX=XX[,1],YY=XX[,2],co.rec=as.vector(co.rec))
  DB.corec <- DB.corec[complete.cases(DB.corec[ ,4]),]
  
}



rm(XX,YY,SPECIES)

# Arrange Gower matrix

d.functional <- as.matrix(gower.mat)
d.functional <- d.functional[order(rownames(d.functional)), ]
d.functional <- t(d.functional)
d.functional <- d.functional[order(rownames(d.functional)), ]

# Calculate distances per pairs of species
XX = c(); SPECIES = c()
for (i in 1:ncol(d.functional)) {
  YY=cbind(rep(rownames(d.functional)[i],ncol(d.functional)), rownames(d.functional))
  XX=rbind(XX,YY)
  SPECIES=c(SPECIES,rep(rownames(d.functional)[i],ncol(d.functional)))
}

DB.fun <- data.frame(ID=paste0(XX[,1],XX[,2]),XX=XX[,1],YY=XX[,2],d.functional=as.vector(d.functional))
DB.fun <- DB.fun[complete.cases(DB.fun[ ,4]),]

rm(d.functional,XX,YY,SPECIES)

co.ocurrence <- merge(DB.corec,DB.fun,by="ID")


plot(co.ocurrence$co.rec,co.ocurrence$d.functional, data = co.ocurrence[which (co.ocurrence$co.rec > 0),])




###---------- Null modeling  ###################################################################


richness <- Alpha.isl[c("island","Species")] # richness per island
species <- row.names(traits05)

traits05 #functional trait matrix, entire species pool


## loop for 1000 replicates 

nreplicates = 100
for (k in 1:10){ 
  
  random.communities.all <- data.frame()   
  

  for (i in 1:nrow(richness)){
    
    richness.this.island <- richness[i,2]
    this.island <- richness[i,1]

    random.communities.island <-  data.frame()
    for (p in 1:nreplicates){
      random.community <- as.vector(sample(species, richness.this.island))
      random.community01 <- data.frame(species=random.community,
                                       replicate = rep(paste0(this.island,p),richness.this.island),
                                       island = rep(this.island,richness.this.island))
      random.communities.island <- rbind(random.community01,random.communities.island)
    }  
    random.communities.all <- rbind(random.communities.island,random.communities.all)
  }
  
  comm.random <- abspres(data=random.communities.all, sites.col="replicate", sp.col = "species", keep.n=F)
  
  row.names(comm.random) <- comm.random[,1]; comm.random <- comm.random[,-1]
  
  kernelFD.random <- kernel.build(comm=comm.random,trait=traits05,abund=FALSE,cores=3,method="gaussian")
  saveRDS(kernelFD.random, file = paste0("kernelFD.random.new",k,".rds"))
  
  
  
  # Extracting hypervolumes for each habitat
  HVr.C <- hypervolume::hypervolume_join(kernelFD.random@HVList[1:nreplicates])
  HVr.F <- hypervolume::hypervolume_join(kernelFD.random@HVList[(nreplicates+1):(2*nreplicates)])
  HVr.G <- hypervolume::hypervolume_join(kernelFD.random@HVList[(2*nreplicates+1):(3*nreplicates)])
  HVr.H <- hypervolume::hypervolume_join(kernelFD.random@HVList[(3*nreplicates+1):(4*nreplicates)])
  HVr.P <- hypervolume::hypervolume_join(kernelFD.random@HVList[(4*nreplicates+1):(5*nreplicates)])
  HVr.T <- hypervolume::hypervolume_join(kernelFD.random@HVList[(5*nreplicates+1):(6*nreplicates)])
  
  
  rich.random <- kernel.alpha(kernelFD.random)
  even.random <- kernel.evenness(comm=kernelFD.random)
  isla.random <- c(rep("C",nreplicates),rep("F",nreplicates),rep("G",nreplicates),
                   rep("H",nreplicates),rep("P",nreplicates),rep("T",nreplicates))
  
  Alpha.null <- data.frame(island=isla.random,Richness=rich.random,
                           Even=even.random)
  write.csv2(Alpha.null,paste0("Alpha.null.new",k,".csv"),sep=";", dec="." )
  print(paste(k, "of", 10))
}


  Alpha.null <- read.csv2("Alpha.null.new.combined.csv", dec = ".")

## Gran Canaria
Alpha.null.C <- Alpha.null[which(Alpha.null$island=="C"),]
ses(Alpha.isl[1,3], as.vector(Alpha.null.C$Richness))
ses(Alpha.isl[1,5], as.vector(Alpha.null.C$Even))


Null.C.R <- ggplot(Alpha.null.C, aes(x=Richness)) +
  geom_histogram(fill="#504C13") +
  geom_vline(aes(xintercept=Alpha.isl[1,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="C")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))

Null.C.E <- ggplot(Alpha.null.C, aes(x=Even)) +
  geom_histogram(fill="#504C13") +
  geom_vline(aes(xintercept=Alpha.isl[1,5]),color="black", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=1)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="C")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))

plot(Null.C.R)
plot(Null.C.E)

## Fuerteventura
Alpha.null.F <- Alpha.null[which(Alpha.null$island=="F"),]
ses(Alpha.isl[2,3], as.vector(Alpha.null.F$Richness))
ses(Alpha.isl[2,5], as.vector(Alpha.null.F$Even))


Null.F.R <- ggplot(Alpha.null.F, aes(x=Richness)) +
  geom_histogram(fill="#993300") +
  geom_vline(aes(xintercept=Alpha.isl[2,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="F")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))

Null.F.E <- ggplot(Alpha.null.F, aes(x=Even)) +
  geom_histogram(fill="#993300") +
  geom_vline(aes(xintercept=Alpha.isl[2,5]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="F")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))

plot(Null.F.R)

## La Gomera
Alpha.null.G <- Alpha.null[which(Alpha.null$island=="G"),]
ses(Alpha.isl[3,3], as.vector(Alpha.null.G$Richness))
ses(Alpha.isl[3,5], as.vector(Alpha.null.G$Even))


Null.G.R <- ggplot(Alpha.null.G, aes(x=Richness)) +
  geom_histogram(fill="#739E49") +
  geom_vline(aes(xintercept=Alpha.isl[3,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="G")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))


Null.G.E <- ggplot(Alpha.null.C, aes(x=Even)) +
  geom_histogram(fill="#739E49") +
  geom_vline(aes(xintercept=Alpha.isl[3,5]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="G")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))

plot(Null.G.E)

## El Hierro
Alpha.null.H <- Alpha.null[which(Alpha.null$island=="H"),]
ses(Alpha.isl[4,3], as.vector(Alpha.null.H$Richness))
ses(Alpha.isl[4,5], as.vector(Alpha.null.H$Even))

Null.H.R <- ggplot(Alpha.null.H, aes(x=Richness)) +
  geom_histogram(fill="#66CCFF") +
  geom_vline(aes(xintercept=Alpha.isl[4,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="H")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))

plot(Null.H.R)

Null.H.E <- ggplot(Alpha.null.H, aes(x=Even)) +
  geom_histogram(fill="#66CCFF") +
  geom_vline(aes(xintercept=Alpha.isl[4,5]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="H")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))
plot(Null.H.E)

## La Palma
Alpha.null.P <- Alpha.null[which(Alpha.null$island=="P"),]
ses(Alpha.isl[5,3], as.vector(Alpha.null.P$Richness))
ses(Alpha.isl[5,5], as.vector(Alpha.null.P$Even))


Null.P.R <- ggplot(Alpha.null.H, aes(x=Richness)) +
  geom_histogram(fill="#3399CC") +
  geom_vline(aes(xintercept=Alpha.isl[5,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="P")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))
plot(Null.P.R)


Null.P.E <- ggplot(Alpha.null.H, aes(x=Even)) +
  geom_histogram(fill="#3399CC") +
  geom_vline(aes(xintercept=Alpha.isl[5,5]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="P")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))
plot(Null.P.E)

## Tenerife
Alpha.null.T <- Alpha.null[which(Alpha.null$island=="T"),]
ses(Alpha.isl[6,3], as.vector(Alpha.null.T$Richness))
ses(Alpha.isl[6,5], as.vector(Alpha.null.T$Even))

Null.T.R <- ggplot(Alpha.null.T, aes(x=Richness)) +
  geom_histogram(fill="#669900") +
  geom_vline(aes(xintercept=Alpha.isl[6,3]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Richness)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional richness")+ylab("Density")+labs(title="T")+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(face="bold"))
plot(Null.T.R)


Null.T.E <- ggplot(Alpha.null.T, aes(x=Even)) +
  geom_histogram(fill="#669900") +
  geom_vline(aes(xintercept=Alpha.isl[6,5]),color="black", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(Even)), color="red", linetype="dashed", size=0.2)+
  xlab("Functional Evenness")+ylab("Density")+labs(title="T")+
  theme_bw()+ theme(legend.position = "none",plot.title = element_text(face="bold"))

plot(Null.T.E)

gridExtra::grid.arrange(Null.H.R,Null.P.R,Null.G.R,Null.T.R,Null.C.R,Null.F.R,
                        Null.H.E,Null.P.E,Null.G.E,Null.T.E,Null.C.E,Null.F.E,
                        ncol=6,nrow=2)


summary.null <- Alpha.null %>% 
   group_by(island) %>%
  summarise(Richness_mean = mean(Richness),
            Richness_se = std(Richness),
            Evenness_mean = mean(Even),
            Evenness_se = std(Even))



#### Supplementary figures: islands—diversity metrics relationships. Figure S01  --------------------


data.islands <- read.csv2("data/islands.csv", dec=".",sep=";")
data.isl <- merge(data.islands,Alpha.isl,by="island",all=T)  
str(data.isl); rm(data.islands)
data.isl[is.na(data.isl)] <- 0 ## <<-- Island richness dataset

#write.csv2(data.isl, "data/islands2.csv")
#data.isl <- read.csv2("data/islands2.csv", dec=".", sep=";")

# A bit of checking
par(mfrow=c(2,2))
hist(data.isl$Species)
hist(data.isl$Richness)
hist(data.isl$Dispersion)
hist(data.isl$Even)
dev.off()


#### Explore correlation between variables
library(dplyr) 

data.isl.cor <- data.isl[,c(2:10)]
psych::pairs.panels(data.isl.cor) # keep altitude, age, surface
data.isl01 <- data.isl[,-c(5:12)] ; rm(data.isl.cor, data.isl)
#data.isl02 <- data.isl01[-6,]


# Preparation of the colour palettes

isl.names <- c("H","F","C","G","P","L","T")



islands.color2 <- c("#01665E","#35978F","#80CDC1",
                    "#C7EAE5","#F6E8C3",
                    "#BF812D","#8C510A")


p1 <-  data.isl01 %>%
  arrange(desc(Surface)) %>%
  mutate(island = factor(island, island)) %>%
  ggplot(aes(x=Altitude, y=Species, size=Surface, color=Age)) +
  scale_colour_gradientn(colours = islands.color2) +
  geom_point(alpha=0.8) +
  ylab("Species richness") +
  scale_size(range = c(0.5, 5), name="Surface (m)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.title = element_text(face="bold"))

p2 <- data.isl01 %>%
  mutate(island = factor(island, island)) %>%
  ggplot(aes(x=Altitude, y=Richness, size=Surface, color=Age)) +
  scale_colour_gradientn(colours = islands.color2) +
  geom_point(alpha=0.8) +
  ylab("Functional richness") +
  scale_size(range = c(0.5,5), name="Surface (m)") +
  theme_bw() +
  theme(legend.position = "none",plot.title = element_text(face="bold"))

# p3 <- data.isl01 %>%
#       mutate(island = factor(island, island)) %>%
#       ggplot(aes(x=Altitude, y=Dispersion, size=Surface, color=Age)) +
#       scale_colour_gradientn(colours = islands.color2) +
#       geom_point(alpha=0.8) +
#       ylab("Functional dispersion") +
#       scale_size(range = c(0.5,5), name="Surface (m)") +
#       theme_bw() +
#       theme(legend.position = "none",plot.title = element_text(face="bold"))

p3 <- data.isl01 %>%
  mutate(island = factor(island, island)) %>%
  ggplot(aes(x=Altitude, y=Even, size=Surface, color=Age)) +
  scale_colour_gradientn(colours = islands.color2) +
  geom_point(alpha=0.8) +
  ylab("Functional evenness") +
  scale_size(range = c(0.5,5), name="Surface (m)") +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(face="bold"))

p4 <- data.isl01 %>%
  mutate(island = factor(island, island)) %>%
  ggplot(aes(x=Altitude, y=Contribution_mean, size=Surface, color=Age)) +
  scale_colour_gradientn(colours = islands.color2) +
  geom_point(alpha=0.8) +
  ylab("Species contribution") +
  xlab("Altitide") +
  scale_size(range = c(0.5,5), name="Functional Richness") +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(face="bold"))

# p6 <- data.isl01 %>%
#   mutate(island = factor(island, island)) %>%
#   ggplot(aes(x=Species, y=Dispersion, size=Surface, color=Age)) +
#   scale_colour_gradientn(colours = islands.color2) +
#   geom_point(alpha=0.8) +
#   ylab("Number of species") +
#   xlab("Functional Dispersion") +
#   scale_size(range = c(0.5,5), name="Functional Richness") +
#   theme_bw() + 
#   theme(legend.position = "none",
#         plot.title = element_text(face="bold"))

# p7 <- data.isl01 %>%
#   mutate(island = factor(island, island)) %>%
#   ggplot(aes(x=Species, y=Even, size=Surface, color=Age)) +
#   scale_colour_gradientn(colours = islands.color2) +
#   geom_point(alpha=0.8) +
#   ylab("Number of species") +
#   xlab("Functional Evennes") +
#   scale_size(range = c(0.5,5), name="Functional Eve") +
#   theme_bw() + 
#   theme(legend.position = "none",
#         plot.title = element_text(face="bold"))

lemon::grid_arrange_shared_legend(p1, p2, p4, p5, ncol = 2, nrow = 2, position='right')
rm(p1,p2,p3,p4,p5,p6,p7)

Alpha.isl$Contribution_mean1 <- Alpha.isl$Contribution_mean*10000
Alpha.isl$Contribution_se1 <- Alpha.isl$Contribution_se*10000




