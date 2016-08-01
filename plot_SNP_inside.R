#script per fare i grafici degli SNP dentro le ROH
library('ggplot2')

all_breed=read.csv('prova',sep=';',header=TRUE)
head(all_breed)

pdf('name_file.pdf',height=12, width=20)  
head(all_breed)

for (a in sort(unique(all_breed$CHR))){ 
    cromo<-subset(all_breed,CHR==a)
    grafico=ggplot(data=cromo,aes(x=POSITION/1000000,y=PERCENTAGE,colour=BREED))  + geom_line() + 
    ggtitle(paste('chr',a,sep=' ')) + scale_y_continuous(limits = c(-0, 100)) + 
      scale_x_continuous(limits = c(-0, max(all_breed$POSITION/1000000)+1))
    print(grafico) 

} 


dev.off()



