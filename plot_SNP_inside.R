#script per fare i grafici degli SNP dentro le ROH
library('ggplot2')

all_breed=read.csv('/Users/Gabriele/Documents/ROHet/eteROH',sep=';')
head(all_breed)

#all_breed<-read.table('"+OUTroh+"',header=T,sep=';')
#pdf('name_file.pdf',height=12, width=20)  
head(all_breed)
seq(as.factor(table(all_breed$CHR)))
for (a in seq(as.factor(table(all_breed$CHR)))){ 
    cromo<-subset(all_breed,CHR==a)
    grafico=ggplot(data=cromo,aes(x=POSITION/1000000,y=PERCENTAGE,colour=BREED))  + geom_line() + 
    ggtitle(paste('chr',a,sep=' ')) + scale_y_continuous(limits = c(-0, 100))
    print(grafico) 
    } 

#dev.off()




