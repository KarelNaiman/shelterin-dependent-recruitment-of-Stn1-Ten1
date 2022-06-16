

puseq_origin_activity<-function(all_ratio_csv, distance, bs_bedgraph, genotype_name) {
  
  all_ratio<-read.csv(all_ratio_csv)
  all_ratio<-all_ratio[complete.cases(all_ratio),] 

  bs<-read.table(bs_bedgraph)
  bs<-bs[bs$V2 >= distance,]
  levels(bs$V1)<-gsub("Chromosome_", "chr", levels(bs$V1)) #replace Chromosome_1 with chr1
  
  all_ratio$avg_diff<-rowMeans(all_ratio[,c(12,13)])
  
  chro<-unique(all_ratio$chromosome)
  bins_needed<-ceiling(distance/300)
  activity<-data.frame(matrix(NA))
  
  for (chromo in chro[1:length(chro)]) {
    
    ori.chr<-all_ratio[all_ratio$chromosome == chromo,]
    
    sites.chr<-bs[bs$V1 == chromo,]
    
    activity.chr<-data.frame(matrix(NA, nrow = ((bins_needed*2)+1), ncol = length(sites.chr[,1])))
    
    if (length(sites.chr[,1]) == 0) {  activity.chr<-0 ; activity.chr<-as.data.frame(activity.chr)  } else { #if there are no sites on that chromosome
      
      for(i in 1:length(sites.chr[,1])) {
        
        mid_site<-((sites.chr[i,2] + sites.chr[i,3] ) /2 )
        diff<-as.vector(round(ori.chr[,2]-mid_site))
        a<-min(abs(diff))
        
        if (length(which(diff == a)) == 0 ) {pos_of_site<-(which(diff == (-a)) -1 )} else
          if (length(which(diff == a)) !=0 ) {pos_of_site<-(which(diff == (a)) -1 )}
        # when you take away the negatives, the distance fail safe, fails.
        if(pos_of_site-bins_needed < 0 ) {activity.chr[,i] <- NA } else if
        (pos_of_site < bins_needed || pos_of_site == bins_needed){activity.chr[,i] <- NA } else if
        (pos_of_site+bins_needed > length(ori.chr[,1]) ) {activity.chr[,i] <- NA } else
          
          
        { activity.chr[,i]<-ori.chr[c((pos_of_site-bins_needed):(pos_of_site+bins_needed)), 14 ] }  
        
      } #close for loop
    }#close else () statement
    
    activity<-cbind(activity,activity.chr)
    
    
  } # close chromo loop 
  activity$bins<-data.frame(bins=c(seq(from= -(bins_needed*300), to = 0, by = 300), seq(from =150, to=(bins_needed*300), by=300)))
  activity$means<-2#clumsyfix
  activity$genotype<-genotype_name
  activity
}#close function

program2<-function(all_ratio_csv1, all_ratio_csv2, bin_size, sites_bedgraph, genotype_name) {
  
  
  data_1<-puseq_origin_activity(all_ratio_csv1, bin_size, sites_bedgraph, genotype_name)
  data_2<-puseq_origin_activity(all_ratio_csv2, bin_size, sites_bedgraph, genotype_name)
  
  
  bs<-read.table(sites_bedgraph)
  bs<-bs[bs$V2 >= bin_size,]
  start<-2
  stop<-(length(bs[,1])+1)
  
  combined_data<-cbind( data_1[,start:stop], data_2[,start:stop] ) 
  
  #standard error of origin activity - standard deviation / square root of the number of samples
  
  combined_length<-(length(bs[,1]) * 2 )
  
  combined_data$standard_dev<-apply(combined_data,1,sd, na.rm=TRUE) #g[,1:884] #hetislands[,1:42]
  combined_data$mean<-apply(combined_data[1:combined_length],1,mean, na.rm=TRUE) #tazIndep[,1:28] #tazdep[,1:12]
  
  combined_data$SEM<-combined_data$standard_dev / (sqrt(combined_length)) 
  
  a<-cbind(data_1$bins, combined_data$mean, combined_data$SEM, data_1$genotype) 
  colnames(a)<-c("bins", "mean", "SEM" , "genotype")
  a }


##########################################Change this part######

#insert two different overview_file_csv files, for the same genotype
#insert bin size in bp, e.g. , for 50 kb insert 50000
#insert binding sites / heterochromatin islands file etc in bedgraph format
#insert genotype name 

wildtype<-program2("8-856-7-655-34WT_all_ratios.csv", "34wt_all_ratios.csv", 6000, "strong_Rif1_BSs.bedgraph", "wildtype")
rif1<-program2("9-YKP17-10-YKP19-34rif_all_ratios.csv", "9-YKP17-42-YKP19-34rif_all_ratios.csv", 6000, "strong_Rif1_BSs.bedgraph", "rif1")
stn1<-program2("3-1332-34-1333-34stn_all_ratios.csv", "34stn_all_ratios.csv", 6000, "strong_Rif1_BSs.bedgraph", "stn1")
stn1rif1<-program2("5-1338-6-1339-34stn-rif_all_ratios.csv", "34stn-rif_all_ratios.csv", 6000, "strong_Rif1_BSs.bedgraph", "stn1rif1")


all_data<-rbind(wildtype, stn1, rif1, stn1rif1) # combine the data in any combination


#use this to plot the final graph
#save graph as PDF and open in CorelDraw to modify

p<-ggplot(NULL,aes(x=bins, y=mean, colour=genotype)) + geom_line(data=all_data)+ 
  ylab("Average origin activity") + xlab("Distance away from rif1BS / bp") + theme_bw()+
  geom_errorbar(data=all_data, aes(ymin=mean-SEM, ymax=mean+SEM), position=position_dodge(width=150), alpha = 0.5) +
  geom_hline(yintercept = 0)+
  scale_colour_manual(values=c("blue3", "forestgreen", "purple", "orangered"))+
  ggtitle("Origin activity around strong Rif1 BSs") +
  ylim(-0.02, 0.025)


p











