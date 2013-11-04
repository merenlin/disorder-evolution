library(ggplot2)
library(scales)

proteins <- read.csv("../data/tables/proteins.txt",na.strings='')
colnames(proteins) <- c("disprotID","uniprotID","name","species","seqlength","numdisorder","discontent")
##############################################################################
# Disorder content plot
ggplot(proteins, aes(x=proteins$discontent,y=..count../sum(..count..))) +
  geom_bar(fill=I("grey30")) +
  xlab("Disorder content") +
  ylab("Number of proteins in Disprot(%)") +
  labs(title = "Disorder content distribution") +
  scale_x_continuous(limits=c(0, 105),expand = c(0,0)) +
  scale_y_continuous(labels = percent_format(),expand = c(0,0),limits=c(0,0.17))

##############################################################################
# Species distribution plot
topspecies = sort(table(proteins$species),decreasing = TRUE)[1:10]
speciesdf = data.frame(labels=labels(topspecies)[[1]],number = topspecies)

ggplot(speciesdf, aes(x=reorder(labels, number),y=number)) +
  geom_bar(fill=I("grey30"),stat="identity") +
  xlab("Species") +
  ylab("Number of proteins in Disprot") +
  labs(title = "Species distribution") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() 

##############################################################################
# Protein length vs number of disordered residues
ggplot(proteins,aes(x=seqlength,y=numdisorder)) + geom_point(alpha = 1/4) + geom_smooth(method="lm") +
  theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(0,5000))  +
  xlab("Protein sequence length") +
  ylab("Number of disordered residues")
