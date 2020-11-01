data <-  read.csv("~/Affinity/Ovarian/ SIOL -Affinity.csv",header = TRUE)
data[lower.tri(data)] <- 0L

methylation.names <- names(data)
methylation.count <- length(methylation.names)

Finalmethynames <- NULL
Finalmethy1names <- NULL
Finalsimilarityvalues <- NULL
Finalgenenames <- NULL


for (i in 1:methylation.count) { 
 print(i)   
    for (j in 1:methylation.count)  {
      similarityValue <- abs(data[i,j])
        
          if((methylation.names[i] != methylation.names[j]) & (similarityValue >= 0.2))  {
          methylation.gene <- toString(methylation.names[i])
          methylation1.gene <- methylation.names[j]      
          genenames <- paste(methylation.gene, methylation1.gene, sep = ".")
          
          Finalmethynames <- rbind (Finalmethynames, methylation.gene)
          Finalmethy1names <- rbind (Finalmethy1names, methylation1.gene)     
          Finalgenenames <- rbind (Finalgenenames, genenames)
          
          Finalsimilarityvalues <- rbind (Finalsimilarityvalues, similarityValue)
          }
     }      
 }
 
Finaldata = data.frame(Finalmethynames=Finalmethynames, Finalmethy1names=Finalmethy1names, Finalgenenames =Finalgenenames, Finalsimilarityvalues = Finalsimilarityvalues)
write.table(Finaldata, "~/AllSimilarities/Ovarian/SIOL.csv", sep=",")



##################################### Overall Mean and MEdian calulation of selected affinities ##############################

data <-  read.csv("~/Affinity/AllSimilarities - Case I/Fused-AllSimilarities.csv",header = TRUE)

data <- data[rev(order(data$Finalsimilarityvalues)),]
head(data)
for( num in seq(10,80,by=10) )
{
	cutoff <- round(nrow(data)*num/100)
	top <- head(data, cutoff)
	median <- median(top[,4])
	mean<- mean(top[,4])
	result <- sprintf('%d\t%f\t%f', num, median, mean)
	print(result)
	write(result, file = "~/Affinity/Mean - Median/Fused-MeanMedian.txt", append = TRUE)
}
              

##################################### Top x % of Similarity values selection  #########################################


data <-  read.csv("~/Affinity/AllSimilarities - Case I/Fused-AllSimilarities.csv",header = TRUE)

data <- data[rev(order(data$Finalsimilarityvalues)),]

cutoff <- round(nrow(data)*60/100)
top <- head(data, cutoff)
head(top)
write.table(top, file = "~/Affinity/Top60%/Fused-top60%.csv", sep = ",")

              
