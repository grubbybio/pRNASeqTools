#!/usr/local/bin/python3

## revFerno: Reverse maps miRferno results to genome co-ordinates
## Updated: version-1.0 12/30/2017
## Property of Meyers Lab at University of Delaware
## Author: kakrana@udel.edu
## Author: pupatel@udel.edu
##  ############

import sys,os,re,time,timeit,csv,glob,string,shutil,operator
import subprocess, multiprocessing
from multiprocessing import Process, Queue, Pool
from operator import itemgetter

#########################################################################
#### USER SETTINGS ######################################################

predTar         = "All.targs.parsed.csv"
Local           = "Y"                           ## 'N' means internal version and 'Y' means external version i.e. local files used instead of miRNA and fasta files from server
GenomeDB        = "xxxxxxxxxx"                  ## DB to get fastaOut file
# PAREdb = 'xxxxxxxx '                          ## Make sure that your Library is in the DB 
genomeFeature   = 1                             ## 0 for gene and 1 for inter; 2 for both
scoreThres      = 5                             ## Score cutoff to filter results especially for intergenic targets

#### STEPS ######################
getCoords       = 0                             ## Get genomic coordinates from genomeDB, if OFF then "gene_coords" file in current folder will be used
                                                ## If 0 then make sure the gene_coords file correspond to genomeFeature selected above or just turn this ON
revMap          = 1                             ## Reverse map miRferno coordinates

##### DEVELOPER OPTIONS ##########
nthread         = 6                             ## Need automatic calculation like nProc
nproc           = "Y"                           ## Used by parallel processing
generateFasta   = "N"                           ## Functionality not required in current version
dataserver      = "xxxxxxxxx" 

### FUNCTIONS ###################

def ConnectToDB(server, infile):

    ##infile values are '0' when you dont want to uplaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE bioinfo_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ('\nTrying to connect to mySQL server on %s' % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='xxxxxxxx', passwd='xxxxxxxx')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ('Connection Established\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

def extractCoords(con,db):
    cur= con.cursor()
    
    ###Filteration of RNA types not done at this stage becasue we need all the genes to get intergenic regions
    ###after filteration few genes will be missing and counted as intergenic region
    ###########################################################TEST LINE#################################################################
    #cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))###extra % escapes the % in query
    #####################################################################################################################################
    ##If we remove miRNA here than only gene entry will be removed but its sequence will be covered into intergenic
    if genomeFeature == 1: ##Intergenic
        print("Fetching intergenic coords - If script stucks here for more then a minute, then cancel - ctrl+c and rerun ")
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type NOT LIKE 'mirna' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))### Extra % escapes the % in query
    elif genomeFeature == 0: ##Protein Coding
        print("Fetching gene coords  - If script stucks here for more then a minute, then cancel - ctrl+c and rerun ")
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))### Extra % escapes the % in query

    genome_info = cur.fetchall() ### List with information of all the genes in genome

    ##Check if list is empty
    if not genome_info:
        print ('^^^Gene Coords query returned with an empty list..Exiting^^^')
        sys.exit()
    ####(1, 'c','AT1G01020', 5928, 8737, protein_coding)
    #print (genome_info)
    
    ##Find length of chromosomes to calculate intergenics
    cur.execute('SELECT chr_id, length FROM %s.chromosome_master' % (db))
    chromo_len = cur.fetchall()
    chromo_dict = dict(chromo_len)  ### Made a dict so that chromosome numer could be searched to get chromosome length
    # print ('These are the chromosomes: %s and their length' % (chromo_dict))
    
    genome_info_inter = genome_info ###This list will also hold intergenics
    
    ####GET INTERGENIC REGIONS AND APPEND TO THE LIST WITH GENE COORDS
    alist = []###list maintained to check if first gene on chromosome and strand shows up than intergenic is just the start of gene
    #for gene1, gene2 from izip(*[iter(genome_info)]*2):
    #for i,j in pairwise (genome_info):
    #for gene1, gene2 in it.izip(genome_info[1:], genome_info):
    #for i in range(0, int(len(genome_info))-1):###maybe for i in range(0, len(genome_info -1))
    

    for i in range(0, int(len(genome_info))+1): ## May 20-14 - Modifed from above after extensive trouble shooting - Now the last entry is read and both up and Down calculated
        #print (i)
        gene1 = (genome_info[i])
        gene2 = (genome_info[i+1])
        gene_type = 'inter' ###set to intergenic by default
        #print(gene1,gene2)
        
        ##Remove/skip redundant genes with same start and end......What about Overlapping genes????
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ##gene is same/overlapping consider next gene
            pass
        
        else:
            ##Calculate coordinates of intergenic regions
            ##Is this is first gene on chromosome and strand intergenic region is from position1 - Only chr_id and strand is checked. 
            if tuple(gene1[0:2]) not in alist:
                print ('\n------Caching gene coords for chromosome: %s and strand: %s------\n' % (gene1[0], gene1[1]))
                #print ('Gene1:%s\nGene2:%s' % (gene1,gene2))
                alist.append((gene1[0:2]))
                
                inter_start1 = 1
                inter_end1 = gene1[3]-1###1 nt before start of Gene1 a.k.a the first gene on chromosome in this case
                ##As two genes are read together, the upstream intergenic region gor gene2 must be calculated in same step
                
                ## If both the genes belong to same chr and strand i.e. chromosome has atleast two genes
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: 
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = gene2[3]-1###Till start of second gene
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene2[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])

                    #if inter_name1 == inter_name2:
                    #print ('\nLoop1 - First Gene and more than one gene on this chr and strand')
                    #print (gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type)
                    #print (gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type)
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                
                ## If the first two genes are not from same strand i.e. there is only one gene on this chromosme and strand- code added after Aspragus scaffolds had just one gene
                else: ## intergenic from end of chromosme/scaffold
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = chromo_dict[gene1[0]]###Till end of chromosome
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_down' % gene1[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])
                    #if inter_name1 == inter_name2:
                    
                    #print ('\nLoop2 - First gene on this chromosme and strand but also the only one')
                    #print (gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type)
                    #print (gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type)
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene typ
                
            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]:###If chr_id and strands are equal than find intergenic. These are gene on same chromosme and strand
                    inter_start = gene1[4]+1###End of Gene 1
                    inter_end = gene2[3]-1 ###1 nt before start of gene 2
                    if gene2[1] == 'w': ##Positive strand
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    #print ('\nLoop3 - Not the first gene on chr and strand')
                    #print (gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type)
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: ###That means gene1 is at end of one chromosome and gene 2 is begining of chromosome so we have to extract intergenic at end of one chromosome
                    inter_start = gene1[4]+1###End of gene1
                    inter_end = chromo_dict[gene1[0]]###End of chromosome searched using chromosome id of gene1 from chromosome dictionary
                    if gene1[1] == 'w':##Positive strand end of chromosme
                        inter_name = ('%s_down' % (gene1[2]))
                    else: ##Negative strand first intergenic of chromosme
                        inter_name = ('%s_up' % (gene1[2]))
                        
                    #print ('\nLoop4 - Not the first gene on chromosome and Strand AND the last gene on chromosome')
                    #print (gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)
                    #
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)) ## Chr_id, strand


    ## Additional check for scaffolded genomes, if there are no genes in a scffold it's whole seqeunce will be fetched as intergenic
    if genomeFeature == 1:
        for i in chromo_dict.keys():
            alen = chromo_dict[i]
            # print("Chr:%s | Length:%s" % (i,alen))
            if tuple((i,'c')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_c_all' % (i))
                genome_info_inter.append((i,'c',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

            if tuple((i,'w')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_w_all' % (i))
                genome_info_inter.append((i,'w',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

    
    
    ###Sort the list after adding intergenic regions on on basis of chr_id and strand that is essential while caching chromosme during slicing sequences
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
    #print(genome_info_inter_sort)
    ###Write all cooords for troubleshooting
    all_coords_out = open('all_coords', 'w')
    for i in genome_info_inter_sort:
        all_coords_out.write('%s,%s,%s,%s,%s,%s\n' % (i[0:]))
    all_coords_out.close()
            
    
    ###Filter list to remove unwanted types like miRNA,tRNA,rRNA,snoRNA,snRNA, short or no intergenic
    gene_coords_file = '../coords'####To check wheter coords are printed in chr_id and strand sorted or not
    coords_out = open(gene_coords_file, 'w')
    gene_coords = []## List that will hold genes to fetch, this removes unecessary RNAs and also fix miRNA double entry i.e as gene and miRNA
    for entry in genome_info_inter_sort:
        #print(entry[5])
        ###If RNA type is removed here than that region is not included in analysis but if RNA is removed in mySQL query than only gene is removed and region becomes intergenic
        #if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA'): ##Replace this IF with tuple check i.e miRNA in tuple
        if genomeFeature == 1: ## Inter
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'protein-coding' or entry [5] == 'protein_coding' or entry [5] == 'misc_RNA'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic regon b/w genes or too short than filter
                    #gene_coords.append(entry[:5])
                    #coords_out.write('%s,%s,%s,%s,%s\n' % (entry[0:5]))
                    gene_coords.append(entry[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
        else: ## Protein coding
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'inter'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic regon b/w genes or too short than filter
                    gene_coords.append(entry[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
            
            
    coords_out.close()
    #print(gene_coords)
    return gene_coords ### A list of selected gene coords


def RevMapCoord(ent): #### DOES THIS MODULE WORKS IF USER USES MULTIPLE LIBRARIES - Parallelize
    
    ## Create a dictionary from list of coords to be searched later
    ## Gene_coords structure: 1, 'c','AT1G01020', 5928, 8737, protein_coding   
    
    # print (ent)
    gene_name = ent[1]
    bind_site = ent[2].split('-')
    # cleave_site = int(ent[8])
    
    ## Reverse map co-ordinates ##########################################################
    if Local == 'N':
        print ('**Reverse mapping - Web analysis**')
        ###Check whether gene is from reverse or posiive strand by memebr ship test on dictionary
        if gene_name in coord_dict_wat:
            print ('Entry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
            geno_start = coord_dict_wat[gene_name][1]###Use for reverse mapping of postive genes
    
            #print('Looking for chr_id')
            chr_id = coord_dict_wat[gene_name][0]
            #print('chr_id found')
            strand = 'w' ## AVlilable in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_wat[gene_name][2] ## Gene type
            # new_cleave_site = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
            new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
            new_bind_site_end = (int(geno_start)-1)+int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
        else:
            print ('Entry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
            geno_end = coord_dict_crick[gene_name][2] ### Use for reverse mapping of negative genes
            #print('Looking for chr_id')
            chr_id = coord_dict_crick[gene_name][0]
            #print('chr_id found')
            strand = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_crick[gene_name][2] ##Gene type
            # new_cleave_site = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
            new_bind_site_end = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
            new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
   
    elif Local=='Y': ## 'Y' i.e Local analysis
        print ('**No Reverse mapping of Co-ordinates will be performed - Local analysis**')
        if gene_name in coord_dict_wat:
            print ('Entry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
            geno_start = coord_dict_wat[gene_name][1]###Use for reverse mapping of postive genes
    
            #print('Looking for chr_id')
            chr_id = coord_dict_wat[gene_name][0]
            #print('chr_id found')
            strand = 'w' ## AVlilable in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_wat[gene_name][2] ## Gene type
            # new_cleave_site = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
            new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
            new_bind_site_end = (int(geno_start)-1)+int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
        else:
            print ('Entry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
            geno_end = coord_dict_crick[gene_name][2] ### Use for reverse mapping of negative genes
            #print('Looking for chr_id')
            chr_id = coord_dict_crick[gene_name][0]
            #print('chr_id found')
            strand = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_crick[gene_name][2] ##Gene type
            # new_cleave_site = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
            new_bind_site_end = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
            new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)

    else:
        print("Do nothing")

    print("Rev Mapped: %s,%s,%s,%s" % (str(chr_id),strand,new_bind_site_start,new_bind_site_end))
    rev_mapped_entry = ("%s,%s,%s,%s,%s" % (','.join(ent),str(chr_id),strand,new_bind_site_start,new_bind_site_end))
    
    return rev_mapped_entry

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(lambda args: module(*args), alist2)

def PPResults(module,alist):
    npool = Pool(int(nproc))    
    res = npool.map_async(module, alist)
    if(res.get() == []):
        print("YEP!")
    results = (res.get())
    npool.close()           ### Added by Reza to close hanging 
    return results

### MAIN #######################
def main():

    ### 1. GET COORDS ####
    global con
    print ("\n\nStep 1/4: Getting coords from server or local file")
    if getCoords == 1:
        ## Get coords from GFF file or our genoem DB
        if Local == 'Y':
            if generateFasta == 'Y':
                coords = extractFeatures(genomeFile,gffFile) ## Extracts Coords from GFF3
                # fastaOut = getFASTALocal(genomeFile,coords) ##Creates FASTA file
                global tagLen ## Required later for tag2FASTA step as well
                # unambiguousBaseCounter(fastaOut,tagLen)
            else:
                print("\nThe input FASTA file is considered 'as is' for analysis\n")
                # fastaOut = genomeFile ### Make it better
        else:
            if generateFasta == 'Y':
                con = ConnectToDB(dataserver,0)
                coords = extractCoords(con,GenomeDB)###Get selected genes and all intergenic coords - all coords all also fetched too
                # fastaOut = getFASTAServer(con,GenomeDB,coords)###Extracts FASTA of supplied coords
                # unambiguousBaseCounter(fastaOut,tagLen)
            else:
                print("\nThe input FASTA file is considered 'as is' for analysis\n")
                con = ConnectToDB(dataserver,0)
                coords = extractCoords(con,GenomeDB) ##For reverse mapping
                # fastaOut = genomeFile ### Make it better
    else:
        ## Read coords from gene coords file
        fh_in = open("../coords","r")
        fileRead = fh_in.readlines()
        coords = [] ## List to store coords

        for i in fileRead:
            ent = i.strip("\n").split(",")
            # print("This is read entry",ent)
            coords.append((ent))
        print("Coords read from file:%s entries" % (len(coords)))
        print("Step 1/4: DONE!!\n\n")
        time.sleep(1)






    #### 2. PREPARE DICTIONARY OF COORDS ########################################
    print("Step 2/4: Preparing dictionary of coordinates")

    global coord_dict_wat, coord_dict_crick
    coord_dict_wat = {} ## Dictionary of genes at watson strand, +
    coord_dict_crick = {}###Dictionary of genes at crick strand, - strand
    shutil.rmtree('./revMapped', ignore_errors=True) ## AK Added 
    os.mkdir('./revMapped') ## AK added

    global nproc
    nproc ='1' ## Need better handling
    if Local == 'N':
        
        for i in coords:### gene_coords is a list in script, also written out as file of same name
            # print ("This is a coord:",i)
            #strand = i.split(',')[1] for file###TEST if reading from file
            strand = i[1]
            if strand == 'c':### if entry in reverse strand
                atuple = (i[0],i[4],i[5])
                coord_dict_crick[i[2]] = atuple ## Gene name as key and chr_id,strand, end and gene type as value
            else:
                atuple = (i[0],i[3],i[5])
                coord_dict_wat[i[2]] = atuple   ## Gene name as key and chr_id,strand, end and gene type as value
    
    elif Local == 'Y':

        for i in coords:### gene_coords is a list in script, also written out as file of same name
            # print ("This is a coord:",i)
            #strand = i.split(',')[1] for file###TEST if reading from file
            strand = i[1]
            if strand == 'c':### if entry in reverse strand
                atuple = (i[0],i[3],i[4],i[5])
                coord_dict_crick[i[2]] = atuple ## Gene name as key and chr_id,start, end and gene type as value
            else:
                atuple = (i[0],i[3],i[4],i[5])
                coord_dict_wat[i[2]] = atuple   ## Gene name as key and chr_id,start, end and gene type as value
                print (atuple)
    else:
        print("Do nothing")
        
    print("Strand dictionary made")
    print("Step 2/4: DONE!!\n\n")
    time.sleep(1)
    




    ##### 3. Read the scoring input extend file and change coordinates ##########
    print("Step 3/4: Reading predicted/validated targets in list")
    print("Reading '%s' file from folder" % predTar)
    
    ResFls = [] ## File list
    # ResFls = [file for file in os.listdir('./') if file.endswith ('targs.parsed.csv')]
    ResFls.append(predTar) ## To work on single file
    print(ResFls)
    acount = 0  ## All entries
    bcount = 0  ## That passed user provided score cutoff
    for afile in ResFls:
        fh_in = open('./%s' % afile, 'r')### PARE VALIDATED results
        header = fh_in.readline() ## Waste header
        ScoInpExt = [] ## list of final results or parlist
        for res in fh_in:
            acount+=1
            res_strp = res.strip('\n')
            ent =res_strp.split(',')
            if float(ent[5]) <= scoreThres:
                ScoInpExt.append(ent)
                bcount+= 1
            else:
                # print("Below user provides score threshold")
                pass
    
    print("Total entries in predicted file:%s | Passed user threshold:%s" % (acount,bcount))
    print("Step 3/4: DONE!!\n\n")
    time.sleep(2)

    ### 4.REVERSE MAP ##########
    
    ## Rev Map #######
    print("Step 4/4: Reverse mapping using coords dict and targets list")
    
    ## # TEST- SINGLE CORE - For troubleshooting ----##
    # print('\n***Entering RevMapCoord- Serial***\n')
    # ValidTarGeno = []
    # for i in ScoInpExt:
    #     print("\nEntry for reverse mapping:%s" % (i))
    #     z = RevMapCoord(i)
    #     ValidTarGeno.append(z)
    
    ##  NORMAL- PARALLEL MODE - Uncommnet test above for normal use  
    print('\n***Entering RevMapCoord- parallel***\n')
    print ('**Reverse mapping initiated**\n\n')
    ValidTarGeno = PPResults(RevMapCoord, ScoInpExt) ## Results are in form of list
    
    print ('Reverse mapping complete for:%s\n\n\n' % (afile))
    print("Step 4/4: Done!!\n\n")
    time.sleep(1)

    #### WRITE RESULTS ##############
    print ("Writing Results")
    
    revmapRes = './revMapped/%s_revmapped.csv' % (afile)
    fh_out = open(revmapRes, 'w')
    fh_out.write('%s,Chr,Strand,GenomicBindStart,GenomicBindEnd\n' % (header.strip('\n')))
    for i in ValidTarGeno: ## Write Results from list to file
        if i != 'E13-3-13':##Filter for error 13-13-13 that is small window abundance = 0 and ratio calculation error
            fh_out.write('%s\n' % (i))
        else:
            print (i)
    fh_in.close()
    #fh_in2.close()
    fh_out.close()

### RUN ########################
if __name__ == '__main__':

    if nproc == 'Y':
        nproc = int(multiprocessing.cpu_count()*0.80)
    else:
        nproc == int(nproc)
    start = time.time()
    main()
    end = time.time()
    print ('Complete run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit()


## v01
## Written to revmap internal mirFErno predicted targets

## v01 -> v02 [stable] [Nov-4-2015]
## Added compatibility to sPARTAv1.16 and above
## Modifed user input - now afile is provided for reverse mapping instead of lookup in folder
## Improved user experienece stuff, cleaned up the script

## v0.9 -> v1.0
## Fixed reverse mapping error (same as revFerno) - No valid triggers were being found at 'c' strand - In reverse mapping 
#### function for crick strand genomic start site was beng used instead of genomic end to compute 
#### new coordinates - Fixed by replacing geno_end    = coord_dict_crick[gene_name][1] with
#### geno_end    = coord_dict_crick[gene_name][2] -  see the difference in coordinates


