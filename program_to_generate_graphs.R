

library(ggplot2)
library(gridExtra)
all_ratio<-read.csv("overview_file1.csv")  
all_ratio<-all_ratio[complete.cases(all_ratio),] 
all_ratio<-all_ratio[all_ratio$chromosome == "chr2",] # which chromosome do you want to use?
all_ratio$rightward.fork<-rowMeans(cbind(all_ratio$delta.on.reverse, all_ratio$epsilon.on.forward)) #add rightward moving forks
region<-all_ratio[all_ratio$position > 4470300
                  & all_ratio$position < 4492700
                  ,] # input your region of interest in BP
region$position<-region$position/1000

al_ratio<-read.csv("overview_file2.csv")  
al_ratio<-al_ratio[complete.cases(al_ratio),] 
al_ratio<-al_ratio[al_ratio$chromosome == "chr2",] # which chromosome do you want to use?
al_ratio$rightward.fork<-rowMeans(cbind(al_ratio$delta.on.reverse, al_ratio$epsilon.on.forward)) #add rightward moving forks
region2<-al_ratio[al_ratio$position > 4470300
                  & al_ratio$position < 4492700
                  ,] # input your region of interest in BP
region2$position<-region2$position/1000

###################
topstrand<-ggplot(data = region) +
  xlab(NULL) + ylab(NULL) +
  ylim(0.05,1.05)+ theme_classic() +
  geom_hline(yintercept=0.5,aes(color="grey"),linetype="dashed") +
  expand_limits(y=c(0,1.0)) +
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks=NULL) +
  
  geom_line(aes(y=region[,3], x=region[,2]) ,colour = "blue") + 
  geom_line(aes(y=region[,5], x=region[,2]), colour = "red") +
  
  geom_line(data = region2, aes(y=region2[,3], x=region2[,2]) ,colour = "black") + 
  geom_line(data = region2, aes(y=region2[,5], x=region2[,2]), colour = "green4") +

  annotate("rect", xmin = 4485.664, ymin = 0.05, xmax = 4486.357, ymax = 1, fill = "orange", alpha=0.4) + #STE3
  annotate("rect", xmin = 4487.100, ymin = 0.05, xmax = 4489.804, ymax = 1, fill = "orange", alpha=0.4) + #STE2
  annotate("rect", xmin = 4491.853, ymin = 0.05, xmax = 4492.598, ymax = 1, fill = "orange", alpha=0.4) + #STE1
  annotate("rect", xmin = 4478.628, ymin = 0.05, xmax = 4478.830, ymax = 1, fill = "red", alpha=0.4) + #cenH
  annotate("rect", xmin = 4479.026, ymin = 0.05, xmax = 4479.324, ymax = 1, fill = "red", alpha=0.4) + #cenH
  annotate("rect", xmin = 4476.885, ymin = 0.05, xmax = 4482.899, ymax = 1, fill = "grey", alpha=0.4) #Tlh2

bottomstrand<-ggplot(data = region) + 
  xlab("distance (Kb)") + 
  ylab(NULL) +
  ylim(0.05,1.05)+ theme_classic() +
  geom_hline(yintercept=0.5,aes(color="grey"),linetype="dashed") +
  expand_limits(y=c(0,1.0)) +
  scale_y_continuous(breaks=NULL) +
  #scale_x_continuous(breaks=NULL) +
  geom_line(data = region, aes(y=region[,4], x=region[,2]) ,colour = "blue") + 
  geom_line(data = region, aes(y=region[,6], x=region[,2]), colour = "red") +
  
  geom_line(data = region2, aes(y=region2[,4], x=region2[,2]) ,colour = "black") + 
  geom_line(data = region2, aes(y=region2[,6], x=region2[,2]), colour = "green4") +

  annotate("rect", xmin = 4485.664, ymin = 0.05, xmax = 4486.357, ymax = 1, fill = "orange", alpha=0.4) + #STE3
  annotate("rect", xmin = 4487.100, ymin = 0.05, xmax = 4489.804, ymax = 1, fill = "orange", alpha=0.4) + #STE2
  annotate("rect", xmin = 4491.853, ymin = 0.05, xmax = 4492.598, ymax = 1, fill = "orange", alpha=0.4) + #STE1
  annotate("rect", xmin = 4478.628, ymin = 0.05, xmax = 4478.830, ymax = 1, fill = "red", alpha=0.4) + #cenH
  annotate("rect", xmin = 4479.026, ymin = 0.05, xmax = 4479.324, ymax = 1, fill = "red", alpha=0.4) + #cenH
  annotate("rect", xmin = 4476.885, ymin = 0.05, xmax = 4482.899, ymax = 1, fill = "grey", alpha=0.4) #Tlh2

grid.arrange(topstrand, bottomstrand, nrow=2)

