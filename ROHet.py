#Runs of Heterogosity
#per lanciare il programma (esempio):  
# python ROHet.py --file data/test_code12 --out TEST --missing 1 --homoz 1 --minSNP 10 --length 50  --distance 50


import subprocess as sub
import argparse
import sys 

################################################################################################################################
################################################################################################################################
################################################################################################################################

def try_file(files,roh):
    try:
        ff=open(files)
        if '.ped' in files:
            global recode
            first_line = ff.readline()
            geno=first_line.strip().split()[6:]
            genotype=[geno[x]+geno[x+1] for x in range(0,len(geno)-1,2)]
            geno= set(genotype)            

            recode={'00':'5'}
            for a in geno:
                if a=='00':continue
                if a=='NN':
                    recode[a]='5'
                    continue
                if a[0]==a[1]:
                    if roh:recode[a]='1'
                    else:recode[a]='2'
                elif a[0]!=a[1]:
                    if roh:recode[a]='2'
                    else:recode[a]='1'
                else:
                    print '*'*50,'\n *** ERRORE NEL FILE: '+files+' ***\n','*'*50
                    print '***        ERRORE SCONOSCIUTO NEI GENOTIPI!        ***'
                    sys.exit()

            return (recode)

    except:
        print '*'*50,'\n *** ERRORE NEL FILE: '+files+' ***\n','*'*50
        print '***        I FILE .ped e .map NON SONO STATI TROVATI!        ***'
        print '***    il File .ped e .map devono avere lo stesso nome       ***'
        print "*** Accetto solo 'Example/filename' o 'Example/filename.ped' ***"
        print '-'*50
        sys.exit()


def get_args():
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script for Runs of heterozygous')
    # Add arguments
    parser.add_argument('--files', type=str, help='File name (Example/filename or Example/filename.ped)', required=True)
    parser.add_argument('--roh', action = "store_true", required=False, help = "Runs of Homozygosity")
    parser.add_argument('--out', type=str, help='File OUTPUT name', required=False, default='example')
    parser.add_argument('--length',type=int, help='Minimun ROHet length in kb', required=False, default=1)
    parser.add_argument('--distance',type=int, help='Maximun distance between SNP in kb', required=False, default=100000)
    parser.add_argument('--missing', type=int, help='Maximun missing in ROHet', required=False, default=0)
    parser.add_argument('--homoz', type=int, help='Maximun homozygotes in ROHet', required=False, default=0)
    parser.add_argument('--minSNP', type=int, help='Minimun SNP in ROHet', required=False, default=10)
    
    # Array for all arguments passed to script
    args = parser.parse_args()

    # Assign args to variables
    if ('.ped' or '.map') in args.files:
        fileped=args.files[:-4]+'.ped'
        filemap=args.files[:-4]+'.map'
    else:
        fileped=args.files+'.ped'
        filemap=args.files+'.map'
    try_file(filemap,args.roh)
    try_file(fileped,args.roh)
    recode=try_file

    #Parameters
    length_min1= args.length
    distance= args.distance
    mini_SNP1=args.minSNP
    maxbuffer1=args.homoz
    maxmissing1=args.missing
    file_out=args.out
    roh=args.roh

    # Return all variable values
    return fileped,filemap,length_min1,maxmissing1,maxbuffer1,mini_SNP1,file_out,roh,distance

# Match return values from get_arguments()
# and assign to their respective variables
fileped,filemap, length_min1,maxmissing1,maxbuffer1,mini_SNP1,file_out,roh,distance= get_args()
n_roh=0

#Summary of parameters used
print '*'*50
print '*** Avvio del programma per calcolo delle Runs ***'
print '*'*50
print '\t\t*** PARAMETERS: *** '
if roh:print '\t+++ Runs of Homozygosity +++'
else: print '\t+++ Runs of Heterozygosity +++'
print '\t*** File ped: ',fileped
print '\t*** File map: ',filemap
print '\t*** Minimun length (kb): ',length_min1
print '\t*** Maximun distance between SNP (kb): ',distance
print '\t*** Maximun missing SNP in ROHet: ',maxmissing1
print '\t*** Maximun homozygous SNP in ROHet:',maxbuffer1
print '\t*** Minimun SNPs:',mini_SNP1
print '\t*** Output file name:',file_out
print '_'*50


#####################
### READ MAP FILE ###
#####################
def read_map(map_file,nchrom):
    chromosome=[]
    position=[]
    name=[]
    ##read map file 
    for mapp in open(map_file):
        chro,nam,xx,posiz=mapp.strip().split()
        #if not chro.isdigit() or int(chro)<1 or int(chro)>29: #cambiare i cromosomi!!!
        if not chro.isdigit() or not int(chro) in range(1,nchrom+1): #cambiare i cromosomi!!!
            chromosome.append('skip');position.append('skip')
            continue
        chromosome.append(chro)
        position.append(posiz)
        name.append(nam)
    return(chromosome,position,name)

############################
### RUNS OF HOMOZYGOSITY ###
############################
def ROHet(pedfile,list1,list2):
    
    ##creation of lists
    chromosome,position,name=list1
    mini_SNP1,maxbuffer1,maxmissing1,length_min1,file_out,distance=list2

    ##creation value
    breeds={}
    error=False
    mini_SNP=int(mini_SNP1)
    length_min=int(length_min1)
    maxbuffer=int(maxbuffer1)
    maxmissing=int(maxmissing1)
    #n_roh=0    

    #print chromosome,position,name
    #print mini_SNP,length_min,maxbuffer,maxmissing
    snp1=0

    ##output file
    output=open(file_out,'w')
    output.write('BREED;ANIMAL;CHROMOSOME;COUNT;START;END;LENGTH\n')

    ##def to write to output file
    def write_out(breed,animal,first,last,count,chrom):
        
        diff=int(last)-int(first)
        if count>mini_SNP and diff/1000. > length_min:
            output.write('%s;%s;%s;%s;%s;%s;%s\n'%(breed,animal,chrom,count,first,last,diff))
            global n_roh
            n_roh += 1

        return n_roh

    ##identification of ROH
    for en,line in enumerate(open(pedfile)):
        allpp=[]
        breed,ind,sire,dame,sex,phe,geno=line.strip().split(' ',6)  ## NOTE: This only works with SINGLE SPACE - to do (other delims)
        if not breeds.has_key(breed):breeds[breed]=0
        breeds[breed] += 1
        geno=geno.split()
        ##recode genotype 
        genotype=[recode.get(geno[x]+geno[x+1],'!') for x in range(0,len(geno)-1,2)]
        if genotype.count('!'): error=True 
        count1=0;last1=0;first1=0;lastcrom='1'
        for letter in range(len(genotype)):
            if chromosome[letter]=='skip':continue
            
            snp2=int(position[letter])/1000
            if (snp2-snp1) > distance: #Distance between SNPs
                write_out(breed,ind,first1,last1,count1,lastcrom)
                count1=0
            snp1=snp2

            if chromosome[letter]!=lastcrom:
                ##first write output
                write_out(breed,ind,first1,last1,count1,lastcrom)
                count1=0
                lastcrom=chromosome[letter]
            if count1==0:
                buff=0;first1=0; last1=0; missing=0;
                if genotype[letter]=='1':
                    count1+=1
                    first1=position[letter]
                    continue
                else:continue
            if genotype[letter]=='1':
                count1+=1
                last1=position[letter]
                continue
            elif genotype[letter]=='2': ##Control for heterozygous genotypes
                    buff+=1
                    if buff<=maxbuffer and missing<=maxmissing:
                        count1+=1;continue
                    elif buff > maxbuffer:
                        last1=position[letter-1]
                        ##second write output
                        write_out(breed,ind,first1,last1,count1,lastcrom)
                        count1=0
                        continue
            elif genotype[letter]=='5': ##Control for missings genotypes
                    missing+=1
                    if buff<=maxbuffer and missing<=maxmissing:
                        count1+=1;continue
                    elif missing > maxmissing:
                        last1=position[letter-1]
                        ##third write output
                        write_out(breed,ind,first1,last1,count1,lastcrom)
                        count1=0
    

    print '\t*** Ho trovato in tutto '+str(n_roh)+' ROH ***'
    print '\t*** Analisi Finita ***\n','*'*50
    return (breeds,error)


one=read_map(filemap,90)
list2=[mini_SNP1,maxbuffer1,maxmissing1,length_min1,file_out,distance]
list1=[one]

#avvio del programma
ROHet(fileped,one,list2)


