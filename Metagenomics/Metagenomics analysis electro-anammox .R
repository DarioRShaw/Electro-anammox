##Bin extraction Electro-anammox

#Install packages
#Find more about mmgenome package in https://kasperskytte.github.io/mmgenome2/
install.packages("Rtsne") # Install Rtsne package from CRAN
if(!require(devtools)) install.packages("devtools") # If not already installed
devtools::install_github("jkrijthe/Rtsne") # more info in https://github.com/jkrijthe/Rtsne/tree/openmp

#check for remotes
if(!require(remotes))
  install.packages("remotes")

#install mmgenome2 using remotes
remotes::install_github("kasperskytte/mmgenome2")

#invoke libraries
library (ggplot2)
library(mmgenome2)
library(Rtsne)
library(readxl)
options(scipen = 8)


#Load data
mm <- mmload(
  assembly = "data/assembly.fasta",
  list(
    BRO_1_ANAMMOX = read.csv("data/MQ180323-7_.coverage.csv", header = T)[,c("Name","Average.coverage")],
    BRO_2_SP_ELECTRODE = read.csv("data/MQ180323-8_.coverage.csv", header = T)[,c("Name","Average.coverage")],
    GO_BRO = read.csv("data/MQ180323-9_.coverage.csv", header = T)[,c("Name","Average.coverage")],
    SCA_2_SP_ELECTRODE = read.csv("data/MQ180319-20_.coverage.csv", header = T)[,c("Name","Average.coverage")],
    GO_SCA = read.csv("data/MQ180319-21_.coverage.csv", header = T)[,c("Name","Average.coverage")],
    SCA_INOCULUM = read.csv("data/MQ180319-22_.coverage.csv", header = T)[,c("Name","Average.coverage")]
  ),
  essential_genes = "data/essential.csv", header = T,
  taxonomy = "data/tax.csv", header = T,
  verbose = TRUE,
  kmer_pca = FALSE,
  kmer_BH_tSNE = FALSE
)

#Load paired end connections
paired_ends <- read.csv("data/network.csv")

#general overview of the data frame
mm

#Basic statistics of the metagenomics data
#Basic metagenome stats
#Basic stats of the meta genome libraries and assembly. 
#**Number of scaffolds (#)** is the number of scaffolds in the combined meta genome assembly 
#**Mean GC content (%)** is the mean GC content of all scaffolds in the assembly weighted by 
#scaffold length. **N50 (bp)** is a median statistic that indicates that 50% of the entire assembly 
#is contained in scaffolds equal to or larger than this value (bp).**Length Total (Mbp)** is the total
#combined length of the meta genome assembly in Mbp. **Length Maximum (bp)** is the length of the largest
#scaffold in the assembly. **Length mean** is the mean length of scaffolds in the assembly. 

mmstats(mm)


## Extraction bin Candidatus BROCADIA Differential coverage plot MBR VS ELECTRODE.

#Differential coverage plot of the assembled meta genomic scaffolds (>5000 bp), 
#highlighting the initial extraction of scaffolds for **Bin 1**. 
#The size of the circles represent the length of the scaffolds. 
#Colours of the circles represent phylum level taxonomic classification, 
#based on the essential genes identified on the scaffolds, scaffolds with no colour 
#either contains no essential genes or could not be assigned a phylum level classification. 
#The x and y-axes show the sequencing coverage in the samples (log-scaled).

mmplot(mm, 
       x = "cov_BRO_1_ANAMMOX", 
       y = "cov_BRO_2_SP_ELECTRODE", 
       min_length = 1000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(10, 1000),
       y_limits =c(10, 1000),
       #locator = TRUE, #uncomment to use the locator and return a selection
       )+
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 19), nrow = 8)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")    
    
#Select the bin
selection <- data.frame(cov_BRO_1_ANAMMOX = c(95.1, 140, 189, 133, 86.6),
                        cov_BRO_2_SP_ELECTRODE = c(57.1, 59.8, 225, 361, 100))


#Plot with selection of the scaffold of the bin
mmplot(mm, 
       x = "cov_BRO_1_ANAMMOX", 
       y = "cov_BRO_2_SP_ELECTRODE", 
       min_length = 1000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(10, 300),
       y_limits =c(10, 400),
       #locator = TRUE, #uncomment to use the locator and return a selection
       selection= selection
)+
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 19), nrow = 8)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom") 

#Save the differential coverage plot
ggsave(filename="Brocadia MBR vs Electrode.pdf", width = 8, height = 8, dpi=300) 

#Initial scaffold extraction
mm_subset1  <- mmextract(mm,selection = selection)

#In mmgenome 1 use this code to the initial subset of the Brocadia bin:
#mm_subset1 <-mmextract(mm, selection,exclude = c(131356,104539,31946,76299,96588,82444,119489,60323,69577,102112,94224,102412,60333,72463,52709))

#Intial scaffold extraction stats
mmstats(mm_subset1)

#Plot the extracted scaffolds (Check if there are contaminats belonging to different phylum)
#In this case all the scaffolds belongs to Plantomycetes
mmplot(mm_subset1,
       x = "cov_BRO_1_ANAMMOX",
       y = "cov_BRO_2_SP_ELECTRODE",
       color_by = "phylum")

#Since there are no conatminant scaffolds we can proceed to extarct the bin directly:
Bin_Brocadia <- mmextract(mm_subset1,
                          selection = selection)

#See the statistics of the extarcted bin:
mmstats(Bin_Brocadia)

#Export .fa file of Bin Brocadia 
mmexport(Bin_Brocadia,
         assembly = assembly,
         file = "Bin_Brocadia.fa")

#If you desire further polishing of the Bin you can proceed with following steps
#See variables available for plotting
str(mm)

#Look for any eventual contaminants/See different profiles of coverage
#No contaminants Scaffolds belonging to different phylum in any of the samples with Brocadia:
mmplot_pairs(mm_subset1,
             variables = c("cov_BRO_1_ANAMMOX",
                           "cov_BRO_2_SP_ELECTRODE",
                           "cov_GO_BRO",
                           "gc"),
             color_by = "phylum",
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.4,
             size_scale = 0.7,
             textsize = 4) 


#Using paired-end network to extract repeats
#Until now we have just used coverage profiles to extract scaffolds related to our genome of interest.
#However, some scaffolds might be present in many copies (repeats) and hence have a much higher coverage
#than the rest of the genome. In addition, some scaffolds will by chance have a slightly different coverage
#profile than the rest of the genome and thereby also been missed. Here we generate a network plot of scaffolds
#connected by paired-end reads and use it to extract scaffolds connected to our genome bin of interest.

selection <- data.frame(x = c(5.651, 7.273, 9.761, 11.08, 12.349, 12.256, 9.798, 7.64, 5.967),
                        y = c(12.43, 14.255, 14.536, 12.246, 10.177, 9.289, 8.12, 7.869, 9.672))

p<- mmnetwork(mm_subset1, 
          network = paired_ends,
          color_by = "phylum",
          #locator = TRUE, #uncomment to mark a selection to highlight and extract
          min_connections = 1,
          selection = selection)
p

#In mmgenome 1 use this code for paired-end network extraction 
#mm_subset2<-mmextract_network(subset = mm_subset1,network = paired_ends,original = mm,nconnections =5,type = "direct")
#mmplot_network(data = mm_subset2,network = paired_ends,nconnections = 2,color = "phylum")

#Extarct slected subset from paired-end network connections 
mm_subset2 <- p$data_in_selection #there is no separate mmextract() function for network plots

mm_subset2

#Final trimming and evaluation
selection <- data.frame(cov_BRO_1_ANAMMOX = c(159.497, 204.892, 139.13, 82.411, 97.757),
           cov_BRO_2_SP_ELECTRODE = c(259.574, 142.073, 59.488, 66.512, 171.757))

mmplot(mm_subset2,
       x = "cov_BRO_1_ANAMMOX", 
       y = "cov_BRO_2_SP_ELECTRODE",
       x_scale = "log10",
       y_scale = "log10",
       x_limits = c(1, 250),
       y_limits = c(1, 300),
       #locator = TRUE,
       selection= selection,
       color_by = "phylum",
       fixed_size = 5)

#Save the inlet coverage plot of the exraction
ggsave(filename="Brocadia MBR vs Electrode inlet.pdf", width = 7.7, height = 7.5, dpi=300) 

#Final extarction of the further polished bin
Bin_Brocadia2 <- mmextract(mm_subset2,
                  selection = selection)

mmstats(Bin_Brocadia2)


#Plot hughlighting Bin Brocadia
mmplot(mm,
       x = "cov_BRO_1_ANAMMOX",
       y = "cov_BRO_2_SP_ELECTRODE",
       x_scale = "log10",
       y_scale = "log10",
       highlight_scaffolds = Bin_Brocadia)+
  scale_x_log10(limits = c(10, 1000), breaks = c(10,100,1000)) + 
  scale_y_log10(limits = c(10, 1000), breaks = c(10,100,1000)) +
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 7)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")


#Export .fa file of Bin Brocadia 2
mmexport(Bin_Brocadia2,
         assembly = assembly,
         file = "Bin_Brocadia2.fa")


#
#
#
#
#
#

## Extraction bin Candidatus SCALINDUA Differential coverage plot MBR VS ELECTRODE.

#Differential coverage plot of the assembled meta genomic scaffolds (>5000 bp), 
#highlighting the initial extraction of scaffolds for **Bin 2**. 
#The size of the circles represent the length of the scaffolds. 
#Colours of the circles represent phylum level taxonomic classification, 
#based on the essential genes identified on the scaffolds, scaffolds with no colour 
#either contains no essential genes or could not be assigned a phylum level classification. 
#The x and y-axes show the sequencing coverage in the samples (log-scaled).


#Plot with selection of the scaffold of the bin

selection <- data.frame(cov_SCA_2_SP_ELECTRODE = c(42.5, 174, 188, 41.8),
                        cov_GO_SCA = c(91.3, 97.3, 244, 229))


mmplot(mm, 
       x = "cov_SCA_2_SP_ELECTRODE", 
       y = "cov_GO_SCA", 
       min_length = 5000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(5, 300),
       y_limits =c(5, 300),
       #locator = TRUE, #uncomment to use the locator and return a selection
       selection= selection
)+
  scale_x_log10(limits = c(5, 300), breaks = c(1,10,100,1000)) + 
  scale_y_log10(limits = c(5, 300), breaks = c(1,10,100,1000)) +
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500)) +
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 7)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")

#Save differential coverage plot
ggsave(filename="Scalindua rGO vs Electrode.pdf", width = 8, height = 8, dpi=300) 


#Extarct selected subset Scalindua
mm_subset1  <- mmextract(mm,selection = selection)

#In mmgenome 1 use this code to the initial subset of the Brocadia bin:
#mm_subset1 <-mmextract(mm, selection,exclude = c(131356,104539,31946,76299,96588,82444,119489,60323,69577,102112,94224,102412,60333,72463,52709))

#Intial scaffold extraction stats
mmstats(mm_subset1)

#Plot the extracted scaffolds (Check if there are contaminats belonging to different phylum)
#In this case all the scaffolds belongs to Plantomycetes
mmplot(mm_subset1,
       x = "cov_SCA_2_SP_ELECTRODE",
       y = "cov_GO_SCA",
       color_by = "phylum")

#Using paired-end network to extract repeats

selection <- data.frame(x = c(2.045, 9.155, 7.556, 0.73, -4.598, -5.677, -3.594),
                        y = c(8.229, 2.525, -3.649, -4.846, -3.734, 4.265, 8.031))

p<- mmnetwork(mm_subset1, 
              network = paired_ends,
              color_by = "phylum",
              #locator = TRUE, #uncomment to mark a selection to highlight and extract
              min_connections = 1,
              selection= selection
              )
p

#Extract slected subset from paired-end network connections 
mm_subset2 <- p$data_in_selection #there is no separate mmextract() function for network plots

mm_subset2



#Final trimming and evaluation

selection <-data.frame(cov_SCA_2_SP_ELECTRODE = c(92.314, 207.175, 192.167, 84.033, 36.747, 50.584),
                       cov_GO_SCA = c(343.304, 284.49, 84.857, 66.309, 115.162, 262.037))

mmplot(mm_subset2,
       x = "cov_SCA_2_SP_ELECTRODE", 
       y = "cov_GO_SCA",
       x_scale = "log10",
       y_scale = "log10",
       x_limits = c(1, 300),
       y_limits = c(1, 400),
       #locator = TRUE,
       selection = selection,
       color_by = "phylum",
       fixed_size = 5)

#Final extarction of the further polished bin
Bin_Scalindua <- mmextract(mm_subset2,
                           selection = selection)

mmstats(Bin_Scalindua)


#Plot highlighting Bin Scalindua
mmplot(mm, 
       x = "cov_SCA_2_SP_ELECTRODE", 
       y = "cov_GO_SCA", 
       min_length = 5000,
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(5, 300),
       y_limits =c(5, 300),
       #locator = TRUE, #uncomment to use the locator and return a selection
       highlight_scaffolds = Bin_Scalindua
)+
  scale_x_log10(limits = c(5, 300), breaks = c(1,10,100,1000)) + 
  scale_y_log10(limits = c(5, 300), breaks = c(1,10,100,1000)) +
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500)) +
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19), nrow = 7)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom")

#Export .fa file of Bin Scalindua
mmexport(Bin_Scalindua,
         assembly = assembly,
         file = "Bin_Scalindua.fa")


#
#
#
#
#
#


## Differential coverage plot BROCADIA MBR VS rGO

selection <- data.frame(cov_BRO_1_ANAMMOX = c(128.533, 202.89, 155.548, 97.872, 74.525, 90.805),
                        cov_GO_BRO = c(251.058, 160.088, 113.434, 102.801, 169.35, 225.929))

mmplot(mm, 
       x = "cov_BRO_1_ANAMMOX", 
       y = "cov_GO_BRO", 
       min_length = 2000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(10, 300),
       y_limits =c(10, 300),
       selection = selection
       #locator = TRUE, #uncomment to use the locator and return a selection
)+
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 19), nrow = 8)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom") 

#Save differential coverage plot
ggsave(filename="Brocadia MBR vs rGO.pdf", width = 8, height = 8, dpi=300) 


#
#
#
#
#
#


#Differential coverage plot SCALINDUA MBR vs Electrode

selection <-data.frame(cov_SCA_INOCULUM = c(106.501, 138.165, 121.641, 92.218, 75.132, 81.189),
                       cov_SCA_2_SP_ELECTRODE = c(174.453, 137.343, 84.546, 65.658, 88.689, 140.188))

mmplot(mm, 
       x = "cov_SCA_INOCULUM", 
       y = "cov_SCA_2_SP_ELECTRODE", 
       min_length = 5000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(11, 300),
       y_limits =c(11, 300),
       #locator = TRUE, #uncomment to use the locator and return a selection,
       selection= selection
)+
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 19), nrow = 8)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom") 

#Save differential coverage plot
ggsave(filename="Scalindua MBR vs electrode.pdf", width = 8, height = 8, dpi=300) 


#Differential coverage plot SCALINDUA MBR vs rGO
selection <- data.frame(cov_SCA_INOCULUM = c(94.521, 129.891, 141.568, 100.326, 73.493, 71.099),
                        cov_GO_SCA = c(205.545, 169.749, 124.812, 98.936, 114.986, 170.913))

mmplot(mm, 
       x = "cov_SCA_INOCULUM", 
       y = "cov_GO_SCA", 
       min_length = 1000,
       color_by = "phylum",
       x_scale = "log10",
       y_scale = "log10",
       x_limits =c(11, 300),
       y_limits =c(11, 300),
       #locator = TRUE, #uncomment to use the locator and return a selection,
       selection= selection
)+
  scale_size_area(max_size = 20, breaks = c(10000, 100000, 500000), name = "Scaffold \nLength (kbp)", label =  c(10, 100, 500))+
  scale_color_discrete(name = "Taxonomic \nClassification") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 19), nrow = 8)) +
  guides(size = guide_legend(nrow = 3)) +
  theme(legend.position = "bottom") 

#Save differential coverage plot
ggsave(filename="Scalindua MBR vs rGO.pdf", width = 8, height = 8, dpi=300) 

