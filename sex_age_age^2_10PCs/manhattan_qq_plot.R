# QQ plot-----------------------------------------------------------------------

#install.packages("qqman")
library(qqman)

# Read in association results
gwas_results = read.table("./PAHmerge2_white_standard_outliersRemoved_10PCs.assoc.linear", header=TRUE)

png('qq_10PCs.png', width=6, height=6, units ='in', res=300)
qq(gwas_results$P, main = "Q-Q plot: 10PCs")
dev.off()


## Manhattan plot with odd ratio's legend---------------------------------------
# load packages
library(readr)
library(tidyr)
library(dplyr)
library(ggrepel)
library(data.table)

# Load data
dat <- read.table("./PAHmerge2_white_standard_outliersRemoved_10PCs.assoc.linear", header=TRUE) 


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

# Make the plot
g <- ggplot(don, aes(x=BPcum, y=-log10(P))) + # alpha=0.8
  
  geom_point( aes(color=BETA), size=1.3) + #alpha=0.8
  scale_color_gradientn(
    colours = c("navyblue","#1e90ff", "ghostwhite", "red", "darkred"),
  ) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0.01, 0),
                     limits = c(0, 10),
                     breaks = seq(0,10,2)
  ) +     # remove space between plot area and x axis
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -log10(sig), color = "red") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "red") +
  
  # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  # geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label= topgenes$Symbol), size=2,
  #                   segment.color = "transparent") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(color = "gray20")
  ) +
  labs(x = "Chromosome", y = expression(paste("-lo", g[10]," (p-value)",sep="")), color = "beta")
g

png('manhattan_10PCs.png', width=10, height=6, units ='in', res=300)
g
dev.off()




