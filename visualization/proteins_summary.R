library(ggplot2)
library(gridExtra)
library(scales)

proteins <- read.csv("../data/tables/proteins.txt",na.strings='')
colnames(proteins) <- c("disprotID","uniprotID","name","species","seqlength","numdisorder","discontent")

##############################################################################
# DLength distribution plot
p1 <- ggplot(proteins, aes(x=proteins$seqlength)) +                       
  geom_histogram(alpha = 0.9, fill="grey") +
  xlab("Sequence length") +
  ylab("Number of proteins in Disprot") +
  labs(title = "Sequence length distribution") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,550))

##############################################################################
# Disorder content plot
p2 <- ggplot(proteins, aes(x=proteins$discontent,y=..count../sum(..count..))) +
  geom_bar(alpha = 0.9, fill="grey",binwidth=5) +
  xlab("Disorder content(%)") +
  ylab("Number of proteins in Disprot") +
  labs(title = "Disorder content distribution") +
  theme_bw() +
  scale_x_continuous(limits=c(0, 105),expand = c(0,0)) +
  scale_y_continuous(labels = percent_format(),expand = c(0,0),limits=c(0,0.17))

##############################################################################
# Species distribution plot
topspecies = sort(table(proteins$species),decreasing = TRUE)[1:10]
speciesdf = data.frame(labels=labels(topspecies)[[1]],number = topspecies)

p3 <- ggplot(speciesdf, aes(x=reorder(labels, number),y=number)) +
  geom_bar(alpha = 0.9, fill="grey",stat="identity") +
  xlab("Species") +
  ylab("Number of proteins in Disprot") +
  labs(title = "Species distribution") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() 

##############################################################################
# Protein length vs number of disordered residues
p4 <- ggplot(proteins,aes(x=seqlength,y=numdisorder)) + geom_point(alpha = 1/4) + geom_smooth(method="lm") +
  theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(0,5000))  +
  labs(title = "Protein length vs number of disordered residues") +
  xlab("Protein sequence length") +
  ylab("Number of disordered residues")

##############################################################################
#Saving into a file
pdf("../results/disprot_summary.pdf",width=17,height=15)
grid.arrange(p1, p2, p3,p4, ncol=2,nrow=2,as.table =TRUE)
dev.off()