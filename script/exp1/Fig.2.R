# load the packages 
library(ggplot2)
library(ggmagnify)
library(dplyr)
library(ggpubr)


# load the output files of Qualimap
## mearged bam 
AG <- read.delim('../Data/AG.genome_results.txt', header = T)
head(AG)

PG <- read.delim('../Data/PG.genome_results.txt', header = T)
head(PG)

# make the data frame with presence/absence column
dat <- data.frame(AG$Contig, AG$Mean.coverage, PG$Mean.coverage)
colnames(dat) <- c('Contig', 'AG.Mean.coverage', 'PG.Mean.coverage')
head(dat)
str(dat)
dim(dat) # 38822     3


# plot the DE genes 2244
DEgene <- read.csv('../Data/DE_PGgene.csv')

head(DEgene)
dim(DEgene) # 2244    7
DEgene <- DEgene[,1]


# plot for PG specific, pl1
pl1 <- ggplot(dat, aes(x = AG.Mean.coverage, y = PG.Mean.coverage)) +
  geom_point(shape = 16) +
  labs(
    y = 'Mean coverage of PG reads on PG transcriptome',
    x = 'Mean coverage of AG reads on PG transcriptome') +
  geom_point(data = dat[dat$Contig %in% DEgene, ], aes(x = AG.Mean.coverage, y = PG.Mean.coverage), col = "blue", shape = 16) + 
  ggtitle('PG specific transcripts') +
  theme_test() +
  coord_cartesian(xlim = c(0, 2500000), ylim = c(0, 5e+5), clip = "off") +
  geom_magnify(from = c(0, 50000, 0, 100000), to = c(1000000, 2500000, 3e+5, 5e+5), axes = "xy",colour = 'gray')  +
  scale_x_continuous(labels = scales::scientific_format())
pl1



###############
# plot for known SFP genes, pl2
dat <- dat %>%
  mutate(Group = ifelse(Contig == '07_TRINITY_DN37245_c1_g1_i1.p1', 'LyAcp4.1', 
                        ifelse(Contig == '12_TRINITY_DN26589_c0_g1_i2.p1', 'LyAcp4.2',
                               ifelse(Contig == '03_TRINITY_DN8105_c1_g1_i1.p1', 'LyAcp5', 
                                      ifelse(Contig == '08_TRINITY_DN216_c0_g1_i7.p1', 'LyAcp7a.1', 
                                             ifelse(Contig == '11_TRINITY_DN343_c0_g1_i2.p2' ,'LyAcp7a.2',
                                                    ifelse(Contig == '02_TRINITY_DN508_c0_g4_i2.p1', 'LyAcp7b',
                                                           ifelse(Contig == '01_TRINITY_DN4850_c1_g1_i1.p1', 'LyAcp8a.1', 
                                                                  ifelse(Contig == '03_TRINITY_DN2221_c1_g3_i1.p2', 'LyAcp8a.2',
                                                                         ifelse(Contig == '01_TRINITY_DN271_c0_g1_i5.p1', 'LyAcp8b',
                                                                                ifelse(Contig == '02_TRINITY_DN114_c0_g1_i3.p1', 'LyAcp10.1',
                                                                                       ifelse(Contig == '05_TRINITY_DN102_c0_g1_i2.p1', 'LyAcp10.2',
                                                                                              NA))))))))))))

head(dat)
summary(as.factor(dat$Group))         

pl2 <- ggplot(dat, aes(x = AG.Mean.coverage, y = PG.Mean.coverage, color = Group)) +
  geom_point() +
  scale_color_manual(values = c('LyAcp4.1' = 'darkblue','LyAcp4.2' = 'lightgreen', 'LyAcp5' = 'red',
                                'LyAcp7a.1' = 'orange', 'LyAcp7a.2' = 'orange4',
                                'LyAcp7b' = 'blue', 'LyAcp8a.1' = 'black', 'LyAcp8a.2' = 'pink',
                                'LyAcp8b' = 'green','LyAcp10.1' = 'darkgreen', 'LyAcp10.2' = 'lightblue'),
                     na.value = NA,
                     breaks = c("LyAcp4.1", "LyAcp4.2","LyAcp5", "LyAcp7a.1", 
                                "LyAcp7a.2","LyAcp7b","LyAcp8a.1","LyAcp8a.2",
                                'LyAcp8b','LyAcp10.1', 'LyAcp10.2')) +
  ggtitle('Previosuly characterized SFP genes') +
  theme_test() +
  labs(x = 'Mean coverage of AG reads on PG transcriptome',
       y = 'Mean coverage of PG reads on PG transcriptome') +
  coord_cartesian(xlim = c(0, 2500000), ylim = c(0, 5e+5), clip = "off") +
  theme(legend.position = c(0.85, 0.8),legend.title = element_blank()) +
  scale_x_continuous(labels = scales::scientific_format())
pl2 


###########
#plot for DE genes 
# candidate SFP genes from DEexp and PGAG
DE.sfp <- read.csv('../Data/candidate.sfp.gene.name.csv')
head(DE.sfp)
DE.sfp <- DE.sfp$sfp.gene


pl3 <- ggplot(dat, aes(x = AG.Mean.coverage, y = PG.Mean.coverage)) +
  geom_point(shape = 16) +
  labs(
    y = 'Mean coverage of PG reads on PG transcriptome',
    x = 'Mean coverage of AG reads on PG transcriptome') +
  geom_point(data = dat[dat$Contig %in% DE.sfp, ], aes(x = AG.Mean.coverage, y = PG.Mean.coverage), col = "orange", shape = 16) + 
  ggtitle('PG specific-, DE transcripts') +
  theme_test() +
  coord_cartesian(xlim = c(0, 2500000), ylim = c(0, 5e+5), clip = "off") +
  geom_magnify(from = c(0, 50000, 0, 100000), to = c(1000000, 2500000, 3e+5, 5e+5), axes = "xy",colour = 'gray')  +
  scale_x_continuous(labels = scales::scientific_format())

pl3

# combine plots 
ggarrange(pl1, pl2, pl3, nrow = 1,labels = c('A','B','C'))

