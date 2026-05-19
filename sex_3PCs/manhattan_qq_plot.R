# QQ plot-----------------------------------------------------------------------

#install.packages("qqman")
library(qqman)

# Read in association results
gwas_results = read.table("./PAHmerge2_MPAPdiff_white_kiana_original_lambda_2024-03-05.assoc.linear", header=TRUE)

png('qq_3PCs.png', width=6, height=6, units ='in', res=300)
qq(gwas_results$P, main = "Q-Q plot: 3PCs")
dev.off()


##################################################################################################
#Figure 1 Code
##################################################################################################
## Manhattan plot with odd ratio's legend---------------------------------------
# load packages
library(readr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)
library(ggtext)

# Load data
dat <- read.table("./PAHmerge2_MPAPdiff_white_kiana_original_lambda_2024-03-05.assoc.linear", header=TRUE)

# Define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

# Define genes for annotating
topsnp = dat[which.min(dat$P),]

don <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  filter(CHR <= 23) %>%
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot,
          #	OR = exp(Beta)
  ) 

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )

label_snps <- c("kgp4351376", "kgp5910174", "kgp11844398")

don <- don %>%
  mutate(label = case_when(
    SNP == "kgp4351376" ~ "rs79362544 (<i>DAPL1</i>)<br>C allele, MAF=0.04<br>β=9.31 (SE 1.53), p=3.04×10<sup>-9</sup>",
    SNP == "kgp5910174"  ~ "rs8057488 (<i>SNX29</i>)<br>T allele, MAF=0.05<br>β=7.62 (SE 1.27), p=4.95×10<sup>-9</sup>",
    SNP == "kgp11844398" ~ "rs57871170 (<i>MAD1L1</i>)<br>A allele, MAF=0.04<br>β=7.45 (SE 1.30), p=2.00×10<sup>-8</sup>",
    TRUE ~ NA_character_
  ))

# Check if the SNPs exist in the data
cat("SNPs found in data:\n")
print(don %>% filter(SNP %in% label_snps) %>% select(SNP, CHR, BP, BPcum, P, label))

# Check how many labels are non-NA
cat("\nNon-NA labels:", sum(!is.na(don$label)), "\n")

g <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR %% 2)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c("navyblue", "red")) +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0, 12), breaks = seq(0,12,2)) +
  geom_hline(yintercept = -log10(sig), color = "red") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "red") +
  geom_richtext(
     data = subset(don, SNP %in% label_snps),
     aes(label = label),
     size = 3.5,
     nudge_y = 1.2,
     fill = "white",
     label.size = 0.2,
     label.padding = unit(0.3, "lines")
   ) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "gray20")
  ) +
  labs(x = "Chromosome", y = expression(paste("-log"[10]," (p-value)",sep="")))

g
tiff('manhattan_3PCs_figuref1_600_dpi_test.tiff', width=10, height=6, units='in', res=600, compression='lzw')

g
dev.off()

##################################################################################################
#Figure 7 Code
##################################################################################################
# Load data
dat <-read.table("./pahmerge_white_vasoresp_agesexpc123_coexcludenew0.3_05_202.assoc.logistic", header=TRUE)

# Define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

# Define genes for annotating
topsnp = dat[which.min(dat$P),]

don <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  filter(CHR <= 23) %>%
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot,
          #	OR = exp(Beta)
  ) 

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum, na.rm = T) + min(BPcum, na.rm = T) ) / 2 )


g <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Color by odd/even chromosome instead of BETA
  geom_point(aes(color=as.factor(CHR %% 2)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c("navyblue", "red")) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(0, 10),
                     breaks = seq(0,10,2)
  ) +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -log10(sig), color = "red") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "red") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "gray20")
  ) +
  labs(x = "Chromosome", y = expression(paste("-lo", g[10]," (p-value)",sep="")))

g
png('manhattan_3PCs_figuref7_600dpi_test.png', width=10, height=6, units='in', res=600)
g
dev.off()
