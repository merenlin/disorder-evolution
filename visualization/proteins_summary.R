library(ggplot2)
library(gridExtra)
library(scales)

dataset <- 'Yeast'
proteins <- read.csv("../data/tables/mobidb/Saccharomyces+cerevisiae.csv", strip.white=TRUE, row.names=NULL)

##############################################################################
# DLength distribution plot
p1 <- ggplot(proteins, aes(x=proteins$seqlength)) +                       
  geom_histogram(alpha = 0.9, fill="grey") +
  xlab("Sequence length") +
  ylab("Number of proteins") +
  labs(title = paste("Sequence length distribution,",dataset)) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,550))

##############################################################################
# Disorder content plot
discontent <- (proteins$numdisorder / proteins$seqlength) * 100
proteins <- cbind(proteins,discontent)
p2 <- ggplot(proteins, aes(x=proteins$discontent)) +
  geom_bar(alpha = 0.9, fill="grey") +
  xlab("Disorder content(%)") +
  ylab("Number of proteins") +
  labs(title = paste("Disorder content distribution,",dataset)) +
  theme_bw() +
  scale_x_continuous(limits=c(0, 105),expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,350))

##############################################################################
# Species distribution plot
#topspecies = sort(table(proteins$species),decreasing = TRUE)[1:10]
#speciesdf = data.frame(labels=labels(topspecies)[[1]],number = topspecies)

#p3 <- ggplot(speciesdf, aes(x=reorder(labels, number),y=number)) +
#  geom_bar(alpha = 0.9, fill="grey",stat="identity") +
#  xlab("Species") +
#  ylab("Number of proteins in Disprot") +
#  labs(title = "Species distribution") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#  coord_flip() 

##############################################################################
# Protein length vs number of disordered residues
p4 <- ggplot(proteins,aes(x=seqlength,y=numdisorder)) + geom_point(alpha = 1/4) + geom_smooth(method="lm") +
  theme_bw() + scale_x_continuous(expand=c(0,0),limits=c(0,5000))  +
  scale_y_continuous(expand=c(0,0),limits=c(0,1500))  +
  labs(title = paste("Protein length vs number of disordered residues,",dataset)) +
  xlab("Protein sequence length") +
  ylab("Number of disordered residues")

##############################################################################
#Saving summary plots into a file
pdf(paste("../results/summary_",dataset,".pdf",sep=""),width=10,height=20)
grid.arrange(p1, p2, p4, ncol=1,nrow=3,as.table =TRUE)
dev.off()
##############################################################################
disfactor <- cut(proteins$discontent,breaks = c(-1,25,50,75,100), labels=c("0-25%","25-50%","50-75%","75-100%"))  

p1 <- ggplot(proteins, aes(x=proteins$numhomologs, fill = factor(disfactor))) + 
  geom_bar(binwidth=25) +
  xlab("Number of homologs (cuttoff 1e-4)") +
  ylab("Number of proteins in Disprot") +
  labs(title = "HHSearch results on Disprot proteins") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),limits=c(0,3300)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,35))


p2 <- ggplot(proteins, aes(x=(proteins$gaps)*100, fill = factor(disfactor))) +
  geom_bar(binwidth=2) +
  xlab("Average number of disordered gaps in the alignment") +
  ylab("Number of proteins in Disprot") +
  labs(title = "HHblits results on Disprot proteins") +
  theme_bw() +
  scale_x_continuous(expand = c(0,0),limits=c(0,100)) +
  scale_y_continuous(expand = c(0,0),limits=c(0,55))

pdf("../results/hhblits_summary.pdf",width=17,height=6)
grid.arrange(p1, p2, ncol=2,nrow=1,as.table =TRUE)
dev.off()

proteins_subset <- proteins[(proteins$numhomologs>1000)&(proteins$numdisorder>50)&(proteins$gaps<0.5),]
write.csv(proteins_subset, file = "../results/proteins_subset.txt")

