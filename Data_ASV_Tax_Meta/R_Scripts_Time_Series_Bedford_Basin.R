library(tidyr )
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(ggplot2)
library(QsRutils)
library(ggplot2)
library(data.table)
library(ggpubr)
library(plyr)
library(dplyr)
library(phyloseq)
library(knitr)
library(tibble)
library(ggrepel)
library(BiocManager)
library(microbiome)
library(aweek)
library(knitr)
library(ggords)
library(tidyverse)
library(ggrepel)
library(pairwiseAdonis)
library(eulerr)


######## Import data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

######## Remove taxa
Bedford_no_mito = subset_taxa(Bedford, !Kingdom=="NA" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Kingdom=="Eukaryota" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Phylum=="NA" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Order=="Chloroplast" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria" )
#Bedford_no_mito = subset_samples(Bedford_no_mito, Depth < 50) # run this line if you want to remove or keep deep samples

############ have a look at the read distribution
#Make a data frame with a column for the read counts of each sample##### 
sample_sum_df <- data.frame(sum = sample_sums(Bedford_no_mito))
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth - Bedford Basin 16S rRNA ASV data") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

########### rarefaction curves
library(devtools)
devtools::install_github("gauravsk/ranacapa")

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}
p <- ggrare(Bedford_no_mito, step = 10, color = "Month", se = FALSE)

########## Set factor levels

sample_data(Bedford_no_mito)$Month = factor(sample_data(Bedford_no_mito)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))
sample_data(Bedford_no_mito)$Year. = factor(sample_data(Bedford_no_mito)$Year., levels=c("Year_2014", "Year_2015", "Year_2016", "Year_2017"))

######### Set colours
my.cols <- c("darkorchid","orange", "blue1",
             "green" ,"black", 
             "cyan", "red", "saddlebrown",
             "goldenrod", "violet",
             "gold1", "slateblue", "chartreuse")

#########################
#Ordination PCA CLR and Euclidian distance
#######################
set.seed(66)
Bedford_clr <- microbiome::transform(Bedford_no_mito, "clr")   
out.pcoa.logt <- ordinate(Bedford_clr, method = "RDA", distance = "euclidean")
evals <- out.pcoa.logt$CA$eig
p3<-plot_ordination(Bedford_clr, out.pcoa.logt, type = "Sample", 
                    color = "Month" ) 
#p3$layers <- p3$layers[-1]
p3+ggtitle("PCA Bedford 16S rRNA ASV")+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.7))+
  theme(legend.key=element_blank())+
  geom_point(size=3)+
  scale_colour_manual(values=my.cols)

######################################
#Alpha diversity and Correlations     #change sample.size = 5000, to your read depth of interest
########################################

Bedford_rar <-rarefy_even_depth(Bedford_no_mito, sample.size = 5000,
                               rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

results = estimate_richness(Bedford_rar, measures =c( 'Chao1', "Shannon"))
write.csv(results, "ASV_alpha.csv")
Richness<-read.csv("ASV_alpha.csv", header=TRUE)
names(Richness)[1] <- "X"
Richness$X <- gsub("[.]", "-", Richness$X)

d = sample_data(Bedford_rar)
d1<- data.frame(d)
write.csv(d1, "META_ASV_alpha.csv")
META<-read.csv("META_ASV_alpha.csv", header=TRUE)
#names(write.csv)[1] <- "X"

METADATA<-merge(Richness,META, by="X", all=TRUE)

#changing name of first column as dates
names(METADATA)[5] <- "Date_Merge"
#checking class of dates
class(METADATA$Date_Merge)
#converting as dataframe
METADATA$Date_Merge <- as.Date(METADATA$Date_Merge, format = "%d/%m/%Y") #%Y/%m/%d - %d/%m/%Y
#checking class again
class(METADATA$Date_Merge)

### stats for alpha
Chao1_Means<-compare_means(Chao1 ~ Month,  data = METADATA, p.adjust.method = "bonferroni")
Shannon_Means<-compare_means(Shannon ~ Month,  data = METADATA, p.adjust.method = "bonferroni")
ddply(METADATA, ~Year, plyr:::summarise, mean = mean(Chao1), sd = sd(Chao1),
      max = max(Chao1), min = min(Chao1))

### Plotting correlations

## Set factor levels
METADATA$Depth..m. = factor(METADATA$Depth..m., levels=c("1m", "5m", "10m"), 
                            labels=c("euphotic","euphotic", "euphotic"))
METADATA$Month = factor(METADATA$Month, levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))

## Set colours
my.cols1 <- c("darkorchid","orange", "blue1",
              "green" ,"black", 
              "cyan", "red", "saddlebrown",
              "goldenrod", "violet",
              "gold1", "slateblue", "chartreuse")

## ggplot
p<-ggplot(data=METADATA, aes(x=day_length , y=Chao1,  fill=Depth..m.)) 
# use #day_length #Chlorophyll.A #Temperature
p <- p + geom_point(aes(color=Depth..m.))+
  geom_point(aes(color=Depth..m.),shape = 21,size = 2,colour = "black")+
  scale_fill_manual(values=my.cols1) +
  scale_colour_manual(values=my.cols1)+
  theme(axis.title.x = element_text(size=16, vjust = 0),
        axis.title.y = element_text(size=16, vjust = 0),
        axis.text.y = element_text(size=22, vjust = 0.5),
        axis.text.x = element_text(size=22, vjust = 0, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Day length")+ylab("Diversity")+ggtitle("...")+
  scale_x_continuous(limits=c(8, 16), breaks=seq(8,16,2))

p
p+stat_smooth(method=lm)+
  stat_cor(method = "pearson", aes(color = Depth..m.), size=6, label.x = 12, p.accuracy = 0.001, r.accuracy = 0.01)+
  guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE)

###################################
# Alpha Diversity and time
################################

Bedford_rar <-rarefy_even_depth(carbom_no_mito, sample.size = 5000,
                               rngseed = 123, replace = TRUE, trimOTUs = TRUE, verbose=TRUE)

### for the rare microbiome, ASVs contributing less < 1% 
## if you rarefy to  5000 reads, then look at ASVs which contribute less than 50 reads. 50 is 1% of 5000 reads
## if you rarefy to 80000 reads, then change 50 to 800
#Bedford_rar= prune_taxa(taxa_sums(Bedford_rar) < 50, Bedford_rar)

sample_data(Bedford_rar)$Depth..m. = factor(sample_data(Bedford_rar)$Depth..m.,levels=c("1m", "5m", "10m"), 
                                           labels=c("euphotic","euphotic", "euphotic"))

p<-plot_richness(Bedford_rar,x="Week", color="Year.", measures=c("Chao1", "Shannon")) #Depth..m.
p + theme_bw()+geom_point(shape = 1,size = 2,colour = "black")+
  scale_color_manual(values=c("blue1", "orange","#06f5fb",
                              "black" ,"red", "#a7cc7b",
                              "#cba9e5", "violet", "#808a8a",
                              "orange", "#ff249d" ))+
  ggtitle("Bedford Basin") +
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=14, vjust = 0.3, angle = 0, hjust=1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"))+
  theme(panel.grid.major = element_blank())+
  stat_smooth(method = "loess",  aes(color = Depth..m.), formula = y ~ x, se = TRUE)

######################################
#Redundancy Analyses

####################################
set.seed(56)
#BB_RDA <- subset_samples(Bedford_no_mito, !Depth..m.=='60m') #subset depth; exclude 60m or leave it in if you want to see the difference between euphotic and aphotic zone
BB_RDA <- subset_samples(Bedford_no_mito, RDA=='Yes') # the RDA doesn't deal with NAs, an easy way around is to make a new column in your metadata file here called RDA and then just fill this column with "Yes" (you have all the data) or "No" (you have missing data) 

sample_data(BB_RDA)$Month = factor(sample_data(BB_RDA)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))

BB_RDA_CLR<-transform(BB_RDA, transform = "clr", target = "OTU")
# extract the OTU table
BB_RDA_CLR_abund<-data.frame(otu_table(BB_RDA_CLR))
BB_RDA_CLR_abund_t<-t(BB_RDA_CLR_abund)

# extract the sample data from the 'phyloseq' object
dat <- data.frame(sample_data(BB_RDA_CLR))
# select the continuous variables to use in the constraint (Y)
# then standardise
X<-dat[c(2,18:24)]
X2<-dat[,c(8,14)]
#X<-na.omit(X) 
X <- X %>% rename(Chla = Chlorophyll.A) 
X <- X %>% rename(NO3 = Nitrate) 
X <- X %>% rename(PO4 = Phosphate) 
X <- X %>% rename(Sal = Salinity) 
X <- X %>% rename(Si = Silicate) 
X <- X %>% rename(Temp = Temperature) 
X <- X %>% rename(DO = oxygen) 
X <- X %>% rename(DL = day_length) 
#etc...
#Now standarise/normailse your metadata
X1 <- decostand(X, method='standardize',na.rm =TRUE)
X3<-cbind(X1,X2)

ord <- rda(BB_RDA_CLR_abund_t)
#ord<-cca(BB_RDA_CLR_abund_t) #if you ant to a CCA
fit <- envfit(ord, X1, perm = 999, na.rm = TRUE)
fit
scores(fit, "vectors")## check the significance of your vectors; then you can only include the significant ones in your RDA

##abundance data prep (sorting according to metadata)
BB_RDA_CLR_t <- t(BB_RDA_CLR_abund)
BB_RDA_CLR.t.sort <- as.data.frame(BB_RDA_CLR_t)
BB_RDA_CLR.t.sort <- cbind(BB_RDA_CLR.t.sort, X3$Month)
BB_RDA_CLR.t.sort <- as.data.frame(BB_RDA_CLR.t.sort)
BB_RDA_CLR.t.sort<-na.omit(BB_RDA_CLR.t.sort) 
BB_RDA_CLR.t.sort <- with(BB_RDA_CLR.t.sort, BB_RDA_CLR.t.sort[order(X3$Month),])
BB_RDA_CLR.t.sort$`X3$Month` <-NULL
BB_RDA_CLR.t.sort <- as.matrix(data.matrix(BB_RDA_CLR.t.sort))
X3.wf <- with(X3, X3[order(Month),])
X3.wf$Month <- factor(X3.wf$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"))
#X3.wf$Depth..m. <- factor(X3.wf$Depth..m., levels=c("1m", "5m", "10m","60m"))
X3.wf.data <- X3.wf[,-which(names(X3.wf) %in% c("Month","Depth..m."))]
gr <- X3.wf$Depth..m.
grl <- factor(gr)

gr2 <- X3.wf$Month
curl.n <- factor(gr2)
#RDA
BB_RDA_CLR.t.sort.rda <- rda(
  BB_RDA_CLR.t.sort ~Chla+NO3+PO4+Sal+Si+Temp+DO+DL, #######here add your significant parameters
  # ,
  data = X3.wf.data) #+DL

print(BB_RDA_CLR.t.sort.rda)

invisible(hist(residuals(BB_RDA_CLR.t.sort.rda), main = ""))#statistics

BB_RDA_CLR.t.sort.rda.anova <- anova.cca(BB_RDA_CLR.t.sort.rda)
print(BB_RDA_CLR.t.sort.rda.anova)

inertia.rda.tot <- BB_RDA_CLR.t.sort.rda$tot.chi#statistics
inertia.rda.tot
inertia.rda.constrained <- BB_RDA_CLR.t.sort.rda$CCA$tot.chi
inertia.rda.constrained
inertia.rda.constrained.prop <- inertia.rda.constrained/inertia.rda.tot#statistics
print(inertia.rda.constrained.prop)

#par(xpd=TRUE)
my.cols <- c("darkorchid","orange", "blue1",
             "green" ,"black", 
             "cyan", "red", "saddlebrown",
             "violet", 
             "gold1", "slateblue", "chartreuse") #4th last ""goldenrod
#Sites
print(ggrda(BB_RDA_CLR.t.sort.rda,group = curl.n, spearrow = NULL, farrow = 0.1, fzoom = 3, 
            ellipse = F,  scaling = 2, spe = F,fcol="red", facol="red",obssize = 3,fsize = 5,sprotate=45)+
        scale_color_manual(name = "Groups",values = my.cols) +
        scale_shape_manual(name = "Groups",values = c(1,2,3,4,5,6,7,8,9,10,11,12)))+ #4th last "3" #for station names include: obslab = T, obssize = 2 ; c(1,7,4,5,6,15,15,11,25,24,3,5)
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.7))

my.cols <- c("darkorchid", "blue", "darkgreen", "orange")
#Estuaries
p<-ggrda(BB_RDA_CLR.t.sort.rda,group = grl, spearrow = NULL, farrow = 0.1, fzoom = 3.8, 
         ellipse = F , scaling = 2, spe = F,fcol="red", facol="red",obssize = 3,fsize = 4,sprotate=45)+
  scale_color_manual(name = "Groups",values = my.cols)+
  scale_shape_manual(name = "Groups",values = c(21,21,21,3,2,2,3,3,3,3,7,7))+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.7))
p 
p+geom_text_repel(mapping = aes(label = curl.n), size = 3, vjust = 1.5,colour = "black",max.overlaps = 10)
grl

##################################
#VENN diagram
############ 
library(eulerr)
library(venneuler)
library(microbiome)
library(microbiomeutilities)
library(MicrobiotaProcess)
library(VennDiagram)
library(ggvenn)

######## Import data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

######## Remove taxa
Bedford_no_mito = subset_taxa(Bedford, !Kingdom=="NA" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Kingdom=="Eukaryota" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Phylum=="NA" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Order=="Chloroplast" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria" )
#Bedford_no_mito = subset_samples(Bedford_no_mito, Depth < 50) # remove deep samples

###
sample_data(Bedford_no_mito)$Depth..m. = factor(sample_data(Bedford_no_mito)$Depth..m., levels=c("1m", "5m", "10m", "60m"), 
                                          labels=c("euphotic","euphotic", "euphotic", "aphotic"))

Bedford_no_mito = subset_samples(Bedford_no_mito, Year== "2017") 
sample_data(Bedford_no_mito)$Month <- factor(sample_data(Bedford_no_mito)$Month, levels=c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul", "Aug", "Sep","Oct", "Nov", "Dec"),
                                             labels=c("Winter","Spring","Spring","Spring","Summer","Summer","Summer","Autumn","Autumn","Autumn","Winter","Winter"))
Bedford_no_mito = subset_samples(Bedford_no_mito, Month== "Winter") 

vennlist <- get_vennlist(Bedford_no_mito,
                         factorNames="Depth..m.")
ggvenn(vennlist, c("euphotic", "aphotic"),   fill_color = c("blue", "orange") )# this will give you a plot with the unique and shared ASVs, which you can then use to compile seasonal trends and yearly trends

### ggvenn will make a figure with unique and shared ASVs for each season (e.g., winter comparison between euphotic and aphotic)
### if you do this comparsion for each year and for each season you end up with the below table.

#Year	  Seasons	Euphotic_only	  Shared	Aphotic_only	Unique_Euphotic_%	Shared_%	Unique_Aphotic_%
#Y_2014	Winter	961	            634	    295	          50.8	            33.5        15.6
#Y_2014	Spring	689	            747	    310	          39.5	            42.8	      17.8
#Y_2014	Summer	397	            426	    389	          32.8	            35.1      	32.1
#Y_2014	Autumn	635	            401	    326	          46.6	            29.4      	23.9
#Y_2015	Winter	1233          	1145	  321	          45.7            	42.4      	11.9
#Y_2015	Spring	764	            811	    288	          41	              43.5      	15.5
#Y_2015	Summer	724	            577	    397	          42.6            	34	        23.4

### The table above is then loaded - overlap_year_ASVs_Months.csv and used to make the violin plots shown in the supplementary material of the manuscript
Overlap <-read.csv('overlap_year_ASVs_Months.csv', row.names=1)

Overlap$Seasons <- factor(Overlap$Seasons, levels=c("Winter", "Spring","Summer", "Autumn"))

colnames(Overlap)
p <- ggplot(Overlap, aes(x=Seasons, y=Unique_Aphotic_.)) + geom_violin()
p +  geom_jitter(height = 0.2, width = 0.1,size=3, aes(colour = factor(Year))) + 
  scale_colour_manual(name="colour", values=c("red", "purple", "orange", "blue"))+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=16, vjust = 0.3, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=16),legend.key=element_blank())+
  ggtitle("Aphotic zone") +
  ylab("Unique ASVs (%)")

###### Another way to visualise the Venn diagram
my.cols <- c("blue","orange")
#library(VennDiagram) #https://github.com/yanlinlin82/ggvenn
venn.diagram(vennlist, main = "Autumn",height = 480 , width = 450 , resolution = 300,
             filename = "test_venn_Autumn_2017.png", alpha = 0.5, fontfamily = "serif", fontface = "bold",cex = 0.5,
             cat.cex = 0.5,  cat.default.pos = "text", cat.dist = c(0.12,0.1), margin = 0.1, lwd = 1,
             imagetype = "png",fill = my.cols)
#


############################################
#SUP05 - Ecotypes 1
##############################################

##Import data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

Bedford_no_mito = subset_taxa(Bedford, !Order=="Chloroplast" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria" )

#### 
BB_subset_CLR <- microbiome::transform(Bedford_no_mito, "clr") #run this one if you want show Ecotyptes
BB_subset_ra_Phylum = subset_taxa(BB_subset_CLR, Genus== "SUP05 cluster")
data_glom<- psmelt(BB_subset_ra_Phylum) # create dataframe from phyloseq object

data_glom$Depth <- factor(data_glom$Depth, 
                                       levels=c("1m", "5m","10m","60m"),
                                       labels=c("Euphotic","Euphotic","Euphotic", "Aphotic"))
names(data_glom)[1] <- "ASVs"

p<-ggplot(data=data_glom, aes(x=Week, y=Abundance, colour=Seq_id, shape=Depth..m.))+
  geom_jitter(size=2)+
  scale_colour_manual(values=c("blue", "magenta1", "gold2", "palegreen", "red"))+
  theme(legend.position="right",panel.grid.major = element_blank()) +
  theme_bw()+ 
  #scale_y_continuous(limits=c(0, 0.04))+
  ylab("CLR")+
  ggtitle("SUP05 Cluster in Bedford Basin")+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=14, vjust = 0.3, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+
p
############################
#SUP05 - Genus seasonal trend
############################

glom <- tax_glom(Bedford_no_mito, taxrank = 'Genus')
BB_subset_CLR <- microbiome::transform(glom, "clr") 
BB_subset_Genus = subset_taxa(BB_subset_CLR, Genus== "SUP05 cluster")

data_glom<- psmelt(BB_subset_Genus) # create dataframe from phyloseq object

data_glom$Depth..m. <- factor(data_glom$Depth..m., levels=c("1m", "5m", "10m", "60m"))
data_glom$Depth..m. <- factor(data_glom$Depth..m., 
                                       levels=c("1m", "5m","10m","60m"),
                                     labels=c("Euphotic","Euphotic","Euphotic", "Aphotic"))
names(data_glom)[1] <- "ASVs"

p<-ggplot(data=data_glom, aes(x=Week, y=Abundance, color=Depth..m., shape=Depth..m.))+
  geom_jitter(size=2)+
  scale_colour_manual(values=c("blue", "orange", "saddlebrown", "forestgreen", "red"))+
  theme(legend.position="none",panel.grid.major = element_blank()) +
  theme_bw()+ 
  ylab("CLR")+
  ggtitle("SUP05 Cluster in Bedford Basin")+
  theme(axis.title.x = element_text(size=16, vjust = 0.3),
        axis.title.y = element_text(size=16, vjust = 0.9),
        axis.text.y = element_text(size=16, vjust = 0.9),
        axis.text.x = element_text(size=14, vjust = 0.3, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))#+
p
p+stat_smooth(method="gam")
#


##############################################
#Indicator Species

#############################################

##Import data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

Bedford_no_mito = subset_taxa(Bedford, !Order=="Chloroplast" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria" )

#### 

BB_glom_Phylum <- tax_glom(Bedford_no_mito, taxrank = 'Phylum')
BB_glom_Phylum_CLR <- microbiome::transform(BB_glom_Phylum, "clr") 

glom_taxo<-tax_table(BB_glom_Phylum_CLR)
write.csv(glom_taxo, "TAX_Phylum_ASV_Bedford_Basin_Indispecies.csv")

Y <- otu_table(BB_glom_Phylum_CLR)
Y<-t(Y)
write.csv(Y, "Phylum_ASV_Bedford_Basin_Indispecies.csv")

dat <- data.frame(sample_data(BB_glom_Phylum_CLR))
write.csv(dat, "Metadata_Phylum_ASV_Bedford_Basin_Indispecies.csv")

#Open "Metadata_Phylum_ASV_Bedford_Basin_Indispecies.csv" file, copy the "Month" column (i.e, Jan, etc) and insert in this in the first column in the "Phylum_ASV_Bedford_Basin_Indispecies.csv" file. Rename this file to "Phylum_ASV_Bedford_Basin_Indispecies_Month.csv". This file is now used as the input for the Indicator analysis 
#install.packages("indicspecies")
library(indicspecies)
pc = read.csv("Phylum_ASV_Bedford_Basin_Indispecies_Month.csv", header= TRUE)

abund = pc[,3:ncol(pc)] # "3" refers to where the sequence data starts. The first column has the sample ids and the second the Months
Month = pc$Month

inv1_duleg = multipatt(abund, Month, func = "r.g", duleg = TRUE, control = how(nperm=9999))
#Indicator species analysis without site groups combinations duleg = TRUE,
summary(inv1_duleg) 
# The summary details the significnat indicator taxa shown here by sequence_ids. The "TAX_Phylum_ASV_Bedford_Basin_Indispecies.csv" file can then be used to add the full taxonomy.

###### Plotting Indicator taxa
##Import data
otu_mat <-"ASV_Table_Bedford_Basin.csv"
tax_mat <- "Taxonomy_species_Silva138.1_Bedford_Basin.csv"
samples_df <- "METADATA_Bedford_Basin.csv"
file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)

Bedford<-read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

Bedford_no_mito = subset_taxa(Bedford, !Order=="Chloroplast" )
Bedford_no_mito = subset_taxa(Bedford_no_mito, !Family=="Mitochondria" )
Bedford_no_mito = subset_samples(Bedford_no_mito, Depth..m.== "60m") # remove or keep aphotic samples

glom <- tax_glom(Bedford_no_mito, taxrank = 'Phylum')
BB_subset_CLR <- microbiome::transform(glom, "clr")
BB_subset_ra_Phylum = subset_taxa(BB_subset_CLR, Phylum== "Verrucomicrobiota")#Other Phylas e.g., Verrucomicrobiota  
data_glom<- psmelt(BB_subset_ra_Phylum) # create dataframe from phyloseq object

data_glom$Depth..m. <- factor(data_glom$Depth..m., levels=c("1m", "5m", "10m", "60m"))
data_glom$Depth..m. <- factor(data_glom$Depth..m., levels=c("1m", "5m","10m", "60m"),
                              labels=c("Euphotic","Euphotic","Euphotic", "Aphotic"))

#data_glom <- subset(data_glom, Abundance > 0)

p<-ggplot(data=data_glom, aes(x=Week, y=Abundance, fill=Depth..m.)) +#Week
  geom_jitter(aes(colour = Depth..m.),shape=21, size=2)+
  scale_fill_manual(values = c("orange", "darkorchid", "blue", "darkgreen","orange"))+ 
  scale_colour_manual(values=c("black", "black", "black", "black"))+
  theme(legend.position="right",panel.grid.major = element_blank()) +
  theme_bw()+ 
  ylab("CLR")+
  ggtitle("Nitrospinota in Bedford Basin")+
  theme(axis.title.x = element_text(size=22, vjust = 0.3),
        axis.title.y = element_text(size=22, vjust = 0.5),
        axis.text.y = element_text(size=22, vjust = 0.5),
        axis.text.x = element_text(size=22, vjust = 0.3, angle = 0, hjust=0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
p+stat_smooth(method="gam", colour="red")
############################


