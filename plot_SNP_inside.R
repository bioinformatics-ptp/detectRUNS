#script per fare i grafici degli SNP dentro le ROH
library('ggplot2')

system2("python", 
        "SNP_inside.py --f detected.ROHom_x.csv --m test_ROHet.map --r testRuns.raw --o snpInRuns"
)


snpInRuns <- read.csv('snpInRuns',sep=';',header=TRUE)

pdf('name_file.pdf',height=12, width=20)  
head(all_breed)

for (chrom in sort(unique(snpInRuns$CHR))) { 
    
    krom <- subset(snpInRuns,CHR==chrom)
    
    p <- ggplot(data=krom, aes(x=POSITION/(10^6), y=PERCENTAGE, colour=BREED))
    p <- p + geom_line() +  ggtitle(paste('chr', chrom, sep=' '))
    p <- p + scale_y_continuous(limits = c(-0, 100))
    p <- p + scale_x_continuous(limits = c(-0, max(snpInRuns$POSITION/(10^6))+1))
    print(p) 

} 


dev.off()



