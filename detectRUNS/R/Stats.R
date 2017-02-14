#####################
## STATISTIC FOR RUNS
#####################

#' Function to reorder data frames by CHROMOSOME
#'
#' The data frame will be reordered according to chromosome:
#' from 1 to n, then X, Y, XY, MT
#' The data frame needs to have a column with name "CHROMOSOME"
#'
#' @param dfx data frame to be reordered (with column "CHROMOSOME")
#'
#' @details
#' Reorder results based on chromosome
#'
#' @return A reordered data frame by chromosome
#' @export
#'
#' @examples
#'

reorderDF <- function(dfx) {

  chr_order <- c((0:99),"X","Y","XY","MT")
  list_chr <- unique(dfx$CHROMOSOME)
  chr_order <- chr_order[chr_order %in% list_chr]
  # new_list_chr <- as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))
  ordered_dfx <- dfx[match(chr_order,dfx$CHROMOSOME),]

  return(ordered_dfx)
}


#' Function to found max position for each chromosome
#'
#'
#' @param mapfile_path Plink map file (for SNP position)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @return A data frame with the max position for chromosome
#' @export
#'
#' @examples
#'

find_Chromosome_length <- function(mapfile_path){

  #LANCIARE IL TEST PER VEDERE SE IL FILE ESISTE
  if(file.exists(mapfile_path)){
    # using data.table to read data
    mapFile <- data.table::fread(mapfile_path, header = F)
  } else {
    stop(paste("file", mapfile_path, "doesn't exists"))
  }

  mappa<-read.table(mapfile_path, header=F)
  name<-c("CHR","SNP_NAME","0","POSITION")
  colnames(mappa)<-name
  maps<-mappa[mappa$POSITION != 0, ] #delete chromosome 0

  #variable
  nchrom=unique(maps$CHR)

  #trova la lunghezza massima per ogni cromosoma
  line=1
  LengthGenome<-matrix(0,ncol=2,nrow=length(nchrom))

  for(i in nchrom){
    #print(paste("Read Chromosome:",i,sep=' '))
    LengthGenome[line,2]<-maps[which.max(maps$POSITION[maps$CHR==i]),4]-maps[which.min(maps$POSITION[maps$CHR==i]),4]
    LengthGenome[line,1]<-i
    line=line+1
  }

  LengthGenome=as.data.frame(LengthGenome, stringsAsFactors = FALSE )
  names<-c("CHROMOSOME","CHR_LENGTH")
  colnames(LengthGenome)<-names

  LengthGenome$CHR_LENGTH = as.numeric(as.vector(LengthGenome$CHR_LENGTH))


  #somma lunghezza genoma
  print(paste("Total genome length:",sum(LengthGenome$CHR_LENGTH),sep=' '))

  return(LengthGenome)
}



#' Plot N. of ROH by sum/mean
#'
#' Function to plot the number of times/percentage a SNP in in a run (population-specific signals)
#' Proportions on the y-axis, bps on the x-axis
#'
#' @param file_runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param path_map map file (.map) file location
#' @param method
#'
#' @return plot of n. of ROH by sum/mean
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @import utils
#'
#' @examples
#'
#'

plot_SumMean_ROH <- function(file_runs,path_map,method=c('sum','mean')){

  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  #Calcolo della lunghezza del genoma e max value nei cromosomi
  LengthGenome=find_Chromosome_length(path_map)

  names(file_runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  #start calculation by method
  if (method=="sum") {
    print("Faccio la somma")
    sum_ROH_genome <- ddply(file_runs,.(IND),summarize,sum=sum(LENGTH)/1000000)
    method="Sum"
  } else {
    print("Faccio la media")
    sum_ROH_genome <- ddply(file_runs,.(IND),summarize,sum=mean(LENGTH)/1000000)
    method="Mean"
  }

  #sum of ROH for Sample
  count_ROH_genome <- count(file_runs,"IND")
  sum_ROH_genome=merge(sum_ROH_genome,count_ROH_genome,by='IND')
  sum_ROH_genome=merge(sum_ROH_genome,unique(file_runs[,c("IND","GROUP")]),by='IND')
  head(sum_ROH_genome)

  #RESULTS!!!!!
  p <- ggplot(data=sum_ROH_genome, aes(x=sum, y=freq, colour=GROUP)) + geom_point()
  p <- p + xlab(paste(method," of ROH in Mbps" , sep='')) + ylab ("Number of ROH for Individual")
  p

}


#' Violin plot sum/mean ROH
#'
#' Function to plot violin plot
#'
#' @param file_runs a data.frame with runs per animal (breed, id, chrom, nSNP, start, end, length)
#' @param method
#'
#' @return Violin plot of n. of ROH by sum/mean
#' @export
#'
#' @importFrom grDevices dev.off pdf
#'
#' @import utils
#'
#' @examples
#'
#'

plot_Violin_ROH <- function(file_runs, method=c("sum","mean")){
  print("inizio a fare i conti")

  names(file_runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")


  #check method
  method <- match.arg(method)
  message(paste("You are using the method:", method))

  #start calculation by method
  if (method=="sum") {
    #use ddply
    mean_roh=ddply(file_runs,.(IND,GROUP),summarize,sum=sum(LENGTH/1000000))
    method="Sum"
  }else{
    mean_roh=ddply(file_runs,.(IND,GROUP),summarize,sum=mean(LENGTH/1000000))
    method="Mean"
  }

  #Violinplot
  p <- ggplot(data=mean_roh, aes(x=GROUP, y=sum, colour=GROUP))
  p <- p + geom_violin (aes(fill=GROUP)) + geom_boxplot(width=0.1)
  p <- p + ylab(paste(method," of ROH in Mbps" , sep=''))
  p

  return(p)

}


#' Function to calculated Froh genome-wide or chromosome-wide
#' Froh = (sum of all ROH for individual) / (Genome length covered by SNP)
#'
#'
#' @param mapFile Plink map file (for SNP position)
#' @param file_runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
#' @param genome_wide vector of TRUE/FALSE (Analisys genome-wide)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @import reshape2
#'
#' @return A data frame with the max position for chromosome
#' @export
#'
#' @examples
#' Froh_new <- Froh_inbreeding(file_runs = runs2, path_map = path_map2, genome_wide = TRUE)


Froh_inbreeding <- function(file_runs,path_map,genome_wide=TRUE){

  LengthGenome=find_Chromosome_length(mapfile_path = path_map)

  names(file_runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")

  #sum of ROH for Sample
  if (genome_wide) {
    print("calcolo Froh su tutto il genoma")

    #RESULTS!!!!!
    Froh <- ddply(file_runs,.(IND),summarize,sum=sum(LENGTH))
    Froh$Froh_genome =  Froh$sum/sum(LengthGenome$CHR_LENGTH)

  } else {
    print("calcolo Froh cromosoma by cromosoma")

    Froh_temp <- ddply(file_runs,.(IND,CHROMOSOME),summarize,sum=sum(LENGTH))
    Froh_temp=merge(Froh_temp,LengthGenome,by='CHROMOSOME')
    Froh_temp$Froh =  Froh_temp$sum/Froh_temp$CHR_LENGTH

    Froh=reshape2::dcast(Froh_temp,IND ~ CHROMOSOME ,value.var = "Froh")
    chr_order <-c((1:99),"X","Y","XY","MT")
    list_chr=unique(Froh_temp$CHROMOSOME)
    new_list_chr=as.vector(sort(factor(list_chr,levels=chr_order, ordered=TRUE)))
    new_list_chr1=paste("Chr_",new_list_chr,sep="")
    new_list_chr=c("IND",new_list_chr)

    #RESULTS!!!!!
    Froh <- Froh[new_list_chr]
    colnames(Froh) <- c('IND',new_list_chr1)
  }

  return(Froh)
}




#' Function to calculated Froh using a ROH-class
#'
#'
#' @param path_map Plink map file (for SNP position)
#' @param file_runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
#' @param Class group of length (in Mbps) by class (defaul: 0-2, 2-4, 4-8, 8-16, >16)
#'
#' @details
#' Create a data frame with the max position in map file (plink format)
#'
#' @return A data frame with the max position for chromosome
#' @export
#'
#' @examples
#'Froh_new_class <- Froh_inbreeding_Class(file_runs = runs2, path_map = path_map2, Class = 2)
#'


Froh_inbreeding_Class <- function(file_runs,path_map,Class=2){

  step_value=Class
  range_mb=c(0,0,0,0,0,99999)
  for (i in seq(from = 2 , to= length(range_mb)-1, by = 1) ){
    range_mb[i]=step_value
    step_value=step_value*2
  }

  #range_mb
  name_CLASS=c(paste(range_mb[1],"-",range_mb[2],sep=''),
               paste(range_mb[2],"-",range_mb[3],sep=''),
               paste(range_mb[3],"-",range_mb[4],sep=''),
               paste(range_mb[4],"-",range_mb[5],sep=''),
               paste(">",range_mb[5],sep=''),
               paste(">",range_mb[6],sep=''))


  print(paste("Class created:"  ,name_CLASS[0:5],sep=' '))

  names(file_runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")
  file_runs$MB <- file_runs$LENGTH/1000000
  file_runs$CLASS=cut(as.numeric(file_runs$MB),range_mb)
  levels(file_runs$CLASS) = name_CLASS
  file_runs$CLASS=factor(file_runs$CLASS)
  table(file_runs$CLASS)

  LengthGenome=find_Chromosome_length(mapfile_path = path_map)

  #sum of ROH for Sample
  print("calcolo Froh per Classe")

  Froh_Class=unique(file_runs[c('GROUP','IND')])
  for (i in range_mb[1:5]){
    print(paste("Class used: >",i,sep=''))

    head(file_runs)

    subset_roh <- file_runs[file_runs$MB >= i,]

    #if subset is empty (no runs for that class) skip/continue
    if(nrow(subset_roh)<1) next

    Froh_temp <- ddply(subset_roh,.(IND),summarize,sum=sum(LENGTH))
    Froh_temp[[paste("Froh_Class_",i,sep="")]] =  Froh_temp$sum/sum(LengthGenome$CHR_LENGTH)
    colnames(Froh_temp)[2]<- paste("Sum_Class_",i,sep="")

    Froh_Class=merge(Froh_Class,Froh_temp,by="IND",all=TRUE)
    #print(head(subset_roh))
  }


  #RESULTS!!!!!
  return(Froh_Class)

}



#' Summary for Runs file
#' Report a list data frame for all analisys
#'
#'
#' @param path_map Plink map file (for SNP position)
#' @param runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
#' @param Class group of length (in Mbps) by class (defaul: 0-2, 2-4, 4-8, 8-16, >16)
#'
#' @details
#' Report a list data frame for all analisys
#'
#' @return FILIPPO
#' @export
#'
#' @examples
#'
#'

runs_summary <- function(runs,mapFile,Class=2, snpInRuns=FALSE,genotype_path){
  print("inizio a controllare i file")
  print(paste("classe usata:",Class))

  n_class=Class

  result_Froh_genome_wide=Froh_inbreeding(file_runs = runs,
                                          path_map = mapFile,
                                          genome_wide = TRUE)
  result_Froh_chromosome_wide=Froh_inbreeding(file_runs = runs,
                                              path_map = mapFile,
                                              genome_wide = FALSE)

  result_Froh_chromosome_class=Froh_inbreeding_Class(file_runs = runs,
                                                     path_map = mapFile,
                                                     Class = n_class)

  names(runs) <- c("GROUP","IND","CHROMOSOME","COUNT","START","END","LENGTH")
  runs$MB=runs$LENGTH/1000000
  #step_value=2

  range_mb=c(0,0,0,0,0,99999)

  for (i in seq(from = 2 , to= length(range_mb)-1, by = 1) ){
    range_mb[i]=n_class
    n_class=n_class*2
  }

  #range_mb
  name_CLASS=c(paste(range_mb[1],"-",range_mb[2],sep=''),
               paste(range_mb[2],"-",range_mb[3],sep=''),
               paste(range_mb[3],"-",range_mb[4],sep=''),
               paste(range_mb[4],"-",range_mb[5],sep=''),
               paste(">",range_mb[5],sep=''),
               paste(">",range_mb[6],sep=''))

  print(paste("Class created:"  ,name_CLASS[0:5],sep=' '))
  runs$CLASS=cut(as.numeric(runs$MB),range_mb)
  levels(runs$CLASS) = name_CLASS
  runs$CLASS=factor(runs$CLASS)

  #RESULTS!!!!!
  summary_ROH_mean1 = ddply(runs,.(GROUP,CLASS),summarize,sum=mean(MB))
  summary_ROH_mean_class = dcast(summary_ROH_mean1,CLASS ~ GROUP ,value.var = "sum")
  levels(summary_ROH_mean_class$CLASS) = name_CLASS[0:5]
  summary_ROH_mean_class

  #RESULTS!!!!!
  summary_ROH_mean_chr1 = ddply(runs,.(GROUP,CHROMOSOME),summarize,sum=mean(MB))
  summary_ROH_mean_chr = reorderDF(dcast(summary_ROH_mean_chr1,CHROMOSOME ~ GROUP ,value.var = "sum"))
  # summary_ROH_mean_chr

  #RESULTS!!!!!
  summary_ROH_count = ddply(runs,.(CLASS,GROUP),nrow)
  names(summary_ROH_count)[3] <- "nRuns"
  summary_ROH_percentage <- summary_ROH_count[,c(1,2)]
  summary_ROH_percentage$pctRuns <- summary_ROH_count$nRuns/sum(summary_ROH_count$nRuns)

  #RESULTS!!!!!
  summary_ROH_count_chr =  ddply(runs,.(CHROMOSOME,GROUP),nrow)
  names(summary_ROH_count_chr)[3] <- "nRuns"
  summary_ROH_count_chr <- reorderDF(summary_ROH_count_chr)
  summary_ROH_percentage_chr <- summary_ROH_count_chr[,c(1,2)]
  summary_ROH_percentage_chr$pctRuns <- summary_ROH_count_chr$nRuns/sum(summary_ROH_count_chr$nRuns)

  result_summary <- list(summary_ROH_count_chr=summary_ROH_count_chr,
                          summary_ROH_percentage_chr=summary_ROH_percentage_chr,
                          summary_ROH_count=summary_ROH_count,
                          summary_ROH_percentage=summary_ROH_percentage,
                          summary_ROH_mean_chr=summary_ROH_mean_chr,
                          summary_ROH_mean_class=summary_ROH_mean_class,
                          result_Froh_genome_wide = result_Froh_genome_wide,
                          result_Froh_chromosome_wide = result_Froh_chromosome_wide,
                          result_Froh_chromosome_class= result_Froh_chromosome_class)

  if (snpInRuns){

    print("Calcolo degli SNP dentro le ROH")
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")
    #read map file
    if(file.exists(mapFile)){
      # using data.table to read data
      mappa <- data.table::fread(mapFile, header = F)
    } else {
      stop(paste("file", mapFile, "doesn't exists"))
    }
    names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
    mappa$x <- NULL

    #Start calculation % SNP in ROH
    print("Calculation % SNP in ROH") #FILIPPO
    all_SNPinROH <- data.frame("SNP_NAME"=character(),
                               "CHR"=integer(),
                               "POSITION"=numeric(),
                               "COUNT"=integer(),
                               "BREED"=factor(),
                               "PERCENTAGE"=numeric(),
                               stringsAsFactors=FALSE)

    # create progress bar
    total <- length(unique(runs$CHROMOSOME))
    print(paste('Chromosome founds: ',total)) #FILIPPO
    n=0
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom <- mappa[mappa$CHR==chrom,]
      snpInRuns <- snp_inside_ROH(runsChrom,mapKrom, genotype_path)
      all_SNPinROH <- rbind.data.frame(all_SNPinROH,snpInRuns)
      n=n+1
      setTxtProgressBar(pb, n)
    }
    close(pb)

    result_summary=append(result_summary,list(SNPinRun = all_SNPinROH))
    print("Calculation % SNP in ROH finish") #FILIPPO
  }



  return(result_summary)
}




#' Summary table  for Runs file
#' Report a list data frame for all analisys
#'
#' @param genotype_path Plink ped file (for SNP position)
#' @param mapfile_path Plink map file (for SNP position)
#' @param runs R object (dataframe) with results per chromosome: subsetted output from RUNS.run()
#' @param threshold value 0 to 1 (default 0.7)
#'
#' @details
#' Table
#'
#' @return FILIPPO
#' @export
#'
#' @examples
#'
#'

table_roh <- function(runs=NULL,SnpInRuns=NULL,genotype_path, mapfile_path, threshold = 0.5) {

  #set a threshold
  threshold_used=threshold*100
  print(paste('Threshold used:',threshold_used))

  #read map file
  if(file.exists(mapfile_path)){
    # using data.table to read data
    mappa <- data.table::fread(mapfile_path, header = F)
  } else {
    stop(paste("file", mapfile_path, "doesn't exists"))
  }
  names(mappa) <- c("CHR","SNP_NAME","x","POSITION")
  mappa$x <- NULL


  if(!is.null(runs) & is.null(SnpInRuns)){
    print('I found only Runs data frame. GOOD!')

    #change colnames in runs file
    names(runs) <- c("POPULATION","IND","CHROMOSOME","COUNT","START","END","LENGTH")

    #Start calculation % SNP in ROH
    print("Calculation % SNP in ROH") #FILIPPO
    all_SNPinROH <- data.frame("SNP_NAME"=character(),
                               "CHR"=integer(),
                               "POSITION"=numeric(),
                               "COUNT"=integer(),
                               "BREED"=factor(),
                               "PERCENTAGE"=numeric(),
                               stringsAsFactors=FALSE)

    # create progress bar
    total <- length(unique(runs$CHROMOSOME))
    print(paste('Chromosome founds: ',total)) #FILIPPO
    n=0
    pb <- txtProgressBar(min = 0, max = total, style = 3)

    #SNP in ROH
    for (chrom in sort(unique(runs$CHROMOSOME))) {
      runsChrom <- runs[runs$CHROMOSOME==chrom,]
      mapKrom <- mappa[mappa$CHR==chrom,]
      snpInRuns <- snp_inside_ROH(runsChrom,mapKrom, genotype_path)
      all_SNPinROH <- rbind.data.frame(all_SNPinROH,snpInRuns)
      n=n+1
      setTxtProgressBar(pb, n)
    }
    close(pb)
    print("Calculation % SNP in ROH finish") #FILIPPO

  }
  else if (is.null(runs) & !is.null(SnpInRuns)) {
    print('I found only SNPinRuns data frame. GOOD!')
    all_SNPinROH=SnpInRuns

  }
  else{
    stop('You gave me Runs and SNPinRuns! Please choose one!')
  }

  #consecutive number
  all_SNPinROH$Number <- seq(1,length(all_SNPinROH$PERCENTAGE))

  #final data frame
  final_table <- data.frame("GROUP"=character(0),"Start_SNP"=character(0),"End_SNP"=character(0),
                            "chrom"=character(0),"nSNP"=integer(0),"from"=integer(0),"to"=integer(0))


  #vector of breeds
  group_list=as.vector(unique(all_SNPinROH$BREED))

  for (grp in group_list){
    print(paste('controllo questa: ',grp))

    #create subset for group/thresold
    group_subset=as.data.frame(all_SNPinROH[all_SNPinROH$BREED %in% c(grp) & all_SNPinROH$PERCENTAGE > threshold_used,])

    #variable
    old_pos=group_subset[1,7]
    snp_pos1=group_subset[1,3]
    Start_SNP=group_subset[1,1]
    snp_count=0

    x=2
    while(x <= length(rownames(group_subset))) {

      snp_count = snp_count + 1
      new_pos=group_subset[x,7]
      old_pos=group_subset[x-1,7]
      chr_old=group_subset[x-1,2]
      chr_new =group_subset[x,2]

      diff=new_pos-old_pos

      if ((diff > 1) | (chr_new != chr_old) | x==length(rownames(group_subset))) {
        # print(paste("Group:",grp,'- Chr:',chr_old,'- n SNP in Runs:',snp_count)) #FILIPPO
        final_table <- rbind.data.frame(final_table,final_table=data.frame("Group"= group_subset[x-1,5],
                                                                           "Start_SNP"=Start_SNP,
                                                                           "End_SNP"=group_subset[x-1,1],
                                                                           "chrom"=group_subset[x-1,2],
                                                                           "nSNP"=snp_count,
                                                                           "from"=snp_pos1,
                                                                           "to"=group_subset[x-1,3]))

        #reset variable
        snp_count=0
        snp_pos1=group_subset[x,3]
        Start_SNP=group_subset[x,1]
      }

      #upgrade x value
      x <- x+1

    }
  }

  rownames(final_table) <- seq(1,length(row.names(final_table)))
  return(final_table)
}


