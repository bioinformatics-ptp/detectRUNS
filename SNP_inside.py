#Runs of Homozygosity

import subprocess as sub
import argparse
import sys


def try_file(files):
    try:
        if '.map' in files:
            ff=open(files)

        elif '.raw' in files:
            ff=open(files)
            global n_breed
            n_animal={}
            n_breed={}
            #chro=[]
            lista=[]
            for line in ff:
                breed,ids,sire,dam,sex,phe,geno=line.strip().split(' ',6) #AGGIUNGERE REPLACE PER I TAB
		if breed=='FID':
			continue
                if not n_animal.has_key(breed):n_animal[breed]=[]
                if not ids in lista:
                    n_animal[breed].append(ids)
                    lista.append(ids)
                else:continue
                    #print 'doppio'         
                
            for n in n_animal:
                n_breed[n]=len(n_animal[n])
            #print n_breed
            


        else:

            ff=open(files)
            print 'leggo il file:',files            
            global nchrom 
            #n_animal={}
            #n_breed={}
            chro=[]
            #lista=[]
            for line in ff:
                if 'CHROMOSOME' in line: continue
                #BREED,ANIMAL,CHROMOSOME,COUNT,START,END,LENGTH=line.strip().split(";")
                if 'chrom' in line:continue
                BREED,ANIMAL,CHROMOSOME,START,END,COUNT,LENGTH=line.strip().split(";")
                chro.append(int(CHROMOSOME))
                #if not n_animal.has_key(BREED):n_animal[BREED]=[]
                #if not ANIMAL in lista:
                #    n_animal[BREED].append(ANIMAL)
                #    lista.append(ANIMAL)
                #else:continue
                    #print 'doppio'

            #print n_animal
            #for n in n_animal:
                #print n_animal[n]
                #print len(n_animal[n]),n
                #n_breed[n]=len(n_animal[n])
            print n_breed
            nchrom=max(chro)


    except:
        print '*'*50,'\n *** ERRORE NEL FILE: '+files+' ***\n','*'*50
        print "***        I FILE delle reads NON E' STATO TROVATO!        ***"
        print "***    Controlla bene che il nome del file sia corretto!!  ***"
        print '-'*50
        sys.exit()


def get_args():
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script for calculated SNP inside ROH')
    # Add arguments
    parser.add_argument('--files', type=str, help='File name (.txt)', required=True)
    parser.add_argument('--map', type=str, help='Map file name (.txt)', required=True)
    parser.add_argument('--ped', type=str, help='Ped file name (.txt)', required=True)
    parser.add_argument('--out', type=str, help='File OUTPUT name', required=False, default='example')


    # Array for all arguments passed to script
    args = parser.parse_args()

    files=args.files
    map=args.map
    ped=args.ped

    try_file(args.ped)
    try_file(args.files)
    print 'letto file map'
    try_file(args.map)

    file_out=args.out

    # Return all variable values
    return files,map,file_out,ped


files,map,file_out,ped= get_args()


############################    
### COUNT SNP INSIDE ROH ###
############################
def snp_inside_ROH(dati_roh,map_file,save,n_breed,nchrom):

    save=open(save,'w')
    save.write('SNP_NAME;CHR;POSITION;COUNT;BREED;PERCENTAGE\n')

    animal={};marker={};final_SNP={};snp_name={}
    for val in open(map_file,'r'):
        chrom,name,x,pos=val.strip().split()
        if int(chrom) > nchrom or int(chrom) < 1:continue 
        if not marker.has_key(chrom):
            marker[chrom]=[]
            final_SNP[chrom]=[]
            snp_name[chrom]=[]
        marker[chrom].append(int(pos))
        final_SNP[chrom].append(0)
        snp_name[chrom].append(name)


        
    for bre in n_breed.keys():
        n_animal=int(n_breed[bre])
        for read in open(dati_roh,'r'):

            #if 'CHROMOSOME' in read:continue
            #breed,animal,chrom,conta,inizio,fine,differenza=read.strip().split(';')
            if 'chrom' in read:continue
            breed,animal,chrom,conta,inizio,fine,differenza=read.strip().split(';')
            #print inizio,fine
            if chrom=='30' or breed !=bre:continue ### THIS WORKS ONLY ON COW

            start=int(inizio);end=int(fine)
            for val,pos in enumerate(marker[chrom]):
                #print  val,pos,start,end
                if pos < start:continue
                if pos > end:break
                final_SNP[chrom][val] +=1    
                
        #print final_SNP

        for t in final_SNP:
            for n,y in enumerate(final_SNP[t]):
                #print snp_name[t][n],t,marker[t][n],final_SNP[t][n],bre,final_SNP[t][n],(n_animal),int(final_SNP[t][n]),int(n_animal)
                save.write('%s;%s;%s;%s;%s;%s\n' % (snp_name[t][n],t,marker[t][n],final_SNP[t][n],bre,int(final_SNP[t][n])*100/int(n_animal)))
                final_SNP[t][n]=0



snp_inside_ROH(files,map,file_out,n_breed,n_breed)
