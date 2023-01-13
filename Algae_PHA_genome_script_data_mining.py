# -*- coding: utf-8 -*-
"""
Created on Fri May  7 12:34:42 2021

@author: Dan
"""

from Bio import SeqIO
from Bio import Entrez
from io import BytesIO
import os
Entrez.email = '***'
desired_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import seaborn as sns
import pandas as pd
import time
import json
import math
import re
import urllib
import numpy as np
#from ete3 import NCBITaxa, Tree, TreeStyle, NodeStyle, faces, AttrFace, ProfileFace, TextFace
#import ete3
ncbi = NCBITaxa()


mastercsvfile='AlgaePHA_CAZy_Families_67131_entries.csv'
mastercsvdf=pd.read_csv(mastercsvfile, header=0, index_col=0)

Entrez.email = '***'

####taxonomy_filters####
#target_taxa_families=['Streptomycetaceae','Nocardiaceae','Streptococcaceae','Bacillaceae','Paenibacillaceae', 'Alteromonadaceae','Halomonadaceae','Xanthomonadaceae','Vibrionaceae','Pseudomonadaceae','Bradyrhizobiaceae','Hyphomicrobiaceae','Rhizobiaceae','Rhodospirillaceae','Sphingomonadaceae','Caulobacteraceae','Rhodobacteraceae','Burkholderiaceae','Oxalobacteraceae','Comamonadaceae']
target_taxa_families_gt50=['Acetobacteraceae',
 'Aeromonadaceae',
 'Alteromonadaceae',
 'Bacillaceae',
 'Bradyrhizobiaceae',
 'Burkholderiaceae',
 'Caulobacteraceae',
 'Comamonadaceae',
 'Erythrobacteraceae',
 'Halomonadaceae',
 'Oxalobacteraceae',
 'Pseudomonadaceae',
 'Rhizobiaceae',
 'Rhodobacteraceae',
 'Rhodospirillaceae',
 'Sphingomonadaceae',
 'Vibrionaceae',
 'Xanthomonadaceae']

filtered_df=mastercsvdf[mastercsvdf['family'].isin(target_taxa_families_gt50)]


#handle=Entrez.esummary(db="genome", id='ALL71072.1', retmode="xml")
#record=Entrez.read(handle)
#print(record)
#
#
#handle=Entrez.esearch(db="genome", term='ALL71072.1', retmode="xml")
#record=Entrez.read(handle)
#print(record)

###test_accession below
#AAA24053.1


accesionlist_raw=list(filtered_df['Accession'])

def split_query_into_200_dict(inlist, splitby):
###this is nessessary as entrex will only accept 200 UIDs/accessions in one go\
### otherwise it times out randomly missing some queries    
    axn_split_dict={}
    index_tuples_check=[]
    maxval=math.ceil(len(inlist)/splitby)-1
    for index in list(range(0,maxval+1)):
        if not index==maxval:
            axn_split_dict[index]=inlist[index*splitby:(index+1)*splitby]
            index_tuples_check.append((index, index+1,index*splitby,(index+1)*splitby))
        elif index==maxval:
            axn_split_dict[index]=inlist[index*splitby:]
            index_tuples_check.append((index, index+1,index*splitby,'end'))
    return axn_split_dict

accession_split_dict=split_query_into_200_dict(accesionlist_raw, 200)

###get ids from accessions, this is required as elink etc return only id from now on\
###this makes it much easier to map data together

id_data_dict={}
errorids=[]
superfails=0
num_retries=4
for idx in accession_split_dict:
    esummaryquery=','.join(accession_split_dict[idx])
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for retrynum in list(range(0,num_retries)):
             try:
                 handlesum=Entrez.esummary(db="protein", id=esummaryquery, retmode="xml")
                 recordsum=Entrez.read(handlesum,validate = False)
                 for record in list(range(0, len(recordsum))):
                     id_data_dict[recordsum[record]['AccessionVersion']]=recordsum[record]['Id']
                 if not len(recordsum)==200:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordsum)) + '. SHOULD BE 200!')
                 time.sleep(2)  
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum))
                 time.sleep(2)
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions") 

###map the IDs to a new column in the master df
filtered_df['Gi']=filtered_df.Accession.map(id_data_dict)



####Link protein to nucleotide database with Retry data
genome_data_dict={}
dbs=['nuccore','genome','assembly']
globorecdict=dict.fromkeys(dbs,[])
globorecs=[[] for x in list(range(0, len(dbs)))]
erroraxns=[]
num_retries=4
superfails=0
for idx in accession_split_dict:
    queries=accession_split_dict[idx]
    ###the elink query must must be a list of strings or it returns info disordered!
    print("Starting to iterate over:  " +str(idx*200) +':' + str((idx+1)*200) + ' of around ' + str(len(accession_split_dict)*200))
    for dbidx in list(range(0, len(dbs))):
        for retrynum in list(range(0,num_retries)):
             try:
                 handle=Entrez.elink(dbfrom="protein", id=queries, db=dbs[dbidx], rettype='acc')
                 record=Entrez.read(handle)
                 globorecdict[dbs[dbidx]]=globorecdict[dbs[dbidx]]+record
                 if not len(record)==200:
                     print('Error with index: ' + str(idx) + ' databases: ' + str(dbs[dbidx]) + " RETURNED ONLY " + \
                           str(len(record)) + '. SHOULD BE 200!')
                 time.sleep(2)     
             except Exception:
                 print('Initial fail with index: ' + str(idx) + ' databases: ' + str(dbs[dbidx]) + " Attempt: " + str(retrynum))
                 time.sleep(2)
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx) + ' databases: ' + str(dbs[dbidx]))
                 erroraxns.append((dbs[dbidx], idx))
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
        else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx) + ' databases: ' + str(dbs[dbidx]))           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over erroraxns manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions") 

###Save dictionary as json
#with open('Genome, nuccore, assembly server responses_14112.json', 'w') as fp:
#    json.dump(globorecdict, fp)    
###Open dictionary as json
with open('Genome, nuccore, assembly server responses_14112.json', 'r') as fp:
    globorecdict2 = json.load(fp)
                              
####create a dictionary of elink data
elink_map_dict={}
errors=[]
for database in dbs:
    if database not in elink_map_dict:
         elink_map_dict[database]={}
    for querynum in list(range(0,len(globorecdict[database]))):
         try:
             queryID=globorecdict[database][querynum]['IdList'][0]
             if not globorecdict[database][querynum]['LinkSetDb'] == []:
                   returneddatabaseID=globorecdict[database][querynum]['LinkSetDb'][0]['Link'][0]['Id']
             else:
                   returneddatabaseID='NA'
         except:
             if globorecdict[database][querynum]['IdList'] == []:
                 errors.append((database, querynum))
                 print('Error: ' + str(database) +' ' + str(querynum))
                 ###retrieve ID from other databases
                 queryID=globorecdict[dbs[0]][querynum]['IdList'][0]
                 print('Retreived ID from nuccore')
         elink_map_dict[database][queryID]=returneddatabaseID        

for dbn in list(range(0,len(dbs))):
    filtered_df[dbs[dbn]]=filtered_df.Gi.map(elink_map_dict[dbs[dbn]])

#filtered_df.to_json('Genome, nuccore, assembly filtered dict_14112.json')
filtered_df=pd.read_json('Genome, nuccore, assembly filtered dict_14112.json')


unique_genomes=set(filtered_df[~(filtered_df['genome'].str.contains('NA'))]['genome'])

###Split the nuccore data into batches of 200, again
nuccore_split_dict=split_query_into_200_dict(list(filtered_df['nuccore']), 200)
##Identify whether the sequences in nuccore are genomic
#e.g. test genomic ID in nuccoire = 692433552
testquery='692433552'
handlesum=Entrez.esummary(db="nuccore", id=testquery, retmode="xml")
recordsum=Entrez.read(handlesum,validate = False)
if bool(re.search('genome',recordsum[0]['Title'].lower()))==True:
   print(recordsum)


genome_data_dict={}
genome_sum_rec_dict={}
errorids=[]
superfails=0
num_retries=4
recsumcum=0
for idx in nuccore_split_dict:
    #remove NAs
    esummaryquery=[x for x in nuccore_split_dict[idx] if not x=='NA']
    lenquery=len(esummaryquery)
    esummaryquery=','.join(esummaryquery)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for retrynum in list(range(0,num_retries)):
             try:
                 handlesum=Entrez.esummary(db="nuccore", id=esummaryquery, retmode="xml")
                 recordsum=Entrez.read(handlesum,validate = False)
                 for record in list(range(0, len(recordsum))):
                     genome_data_dict[recordsum[record]['Id']]=recordsum[record]['Title']
                     genome_sum_rec_dict[recordsum[record]['Id']]=recordsum[record]
                     recsumcum=recsumcum+1
                 if not len(recordsum)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordsum)) + '. SHOULD BE ' +str(lenquery) + '!')
                 time.sleep(2)  
                 print("Total records archived: " + str(recsumcum))
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum))
                 time.sleep(2)
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions") 

###map the IDs to a new column in the master df
filtered_df['Title nuccore']=filtered_df.nuccore.map(genome_data_dict)
#filtered_df.to_json('Genome, nuccore + TITLE, assembly filtered dict_14112.json')
filtered_df=pd.read_json('Genome, nuccore + TITLE, assembly filtered dict_14112.json')

genomedf=filtered_df[filtered_df['Title nuccore'].str.contains('genome',flags=re.IGNORECASE, na=False)]
nongenomedf=filtered_df[~(filtered_df['Title nuccore'].str.contains('genome',flags=re.IGNORECASE, na=False))]
#genomedf.to_csv("Nuccore only Genome dataframe.csv", sep=',')
#genomedf.to_csv("Nuccore only NON Genome dataframe.csv", sep=',')

genomedf=pd.read_csv("Nuccore only Genome dataframe.csv", sep=',')
nongenomedf=pd.read_csv("Nuccore only NON Genome dataframe.csv", sep=',')
#how many genomes?
print(len(set(genomedf['Title nuccore'])))
unique_genome_ids=list(set(genomedf['nuccore']))
###Test esearch for 
testID='12057211'
genometestid='100712'

#1049337857
#1707204610
#1540030467
#261373957
#1120763910]

###link genomes to all coding sequences within it nuccore_protein_cds
maxreturn=10000
handle=Entrez.elink(dbfrom="nuccore", id=testID, db='protein', rettype='ft')
record=Entrez.read(handle)

#turn the list of proteins into a single 
fetchquery=record[0]['LinkSetDb'][0]['Link']
fetchquerystr=[list(y.values())[0] for y in fetchquery]

###get the proteins title to search for keywords again
handle=Entrez.efetch(db="protein", id=fetchquerystr[0:190], rettype='fasta', retmode='xml')
record=Entrez.read(handle)
print(record)

####esearch for genomes within nuccore that contain gene of interest, remove shotgun
queryterm='PhaC[gene] AND Bacteria[Orgn] AND genome[Description]'# NOT shotgun[Description]'
handle=Entrez.esearch(db="nuccore", term=queryterm, retmode="xml", retmax=maxreturn)
record=Entrez.read(handle)
print(record['Count'])
for nuccoreid in record['IdList']:
    if nuccoreid in list(genomedf['nuccore']):






















        
#################################################################
####              START OF GENOME ONNLY ANALYSIS              ###    
#################################################################

##NB this part downloads genome nuccore IDS for all refseqs genomes containing PhaC...
# querystr='PhaC[Gene] AND Bacteria[Orgn] AND genome AND RefSeq'
# temphandle=Entrez.esearch(db='nuccore', term=querystr, retmode="xml", retmax=200000)
# temprecord=Entrez.read(temphandle)
# print(len(temprecord['IdList']))

# querystr='PhaC[Gene] AND Archaea[Orgn] AND genome AND RefSeq'
# temphandle=Entrez.esearch(db='nuccore', term=querystr, retmode="xml", retmax=200000)
# temprecord=Entrez.read(temphandle)
# print(len(temprecord['IdList']))

def esearch_database_for_querystring(genestr, search_term_string, database, maxret):
    querystr=str(genestr)+'[Gene] ' + search_term_string
    temphandle=Entrez.esearch(db=database, term=querystr, retmode="xml", retmax=maxret)
    temprecord=Entrez.read(temphandle)
    print(str(genestr) +': ' + str(temprecord['Count']) + ' returned from search out of retmax: ' + str(maxret))
    return temprecord

#### NON TAX FILTERED ###
ncbi_genomes_containing_PHA_genes_dict={}
taxacountdict={}
allids=[]
pha_gene_list=['PhaA', 'PhaC', 'PhaE', 'PhaG', 'PhaD', 'PhaF', 'FabG', 'PhaR', 'PhaZ', 'PhaC2', 'maoC', 'PhaB',\
     'PhaP', 'PhaJ', 'PhaQ', 'PhaI', 'Pct', 'Pcs', 'PrpE', 'PduP', 'Cat2', 'HadA']
loopcount=0
databasechoice='nuccore'
for gene in pha_gene_list:
        genecount=0
        addstr=' AND ' + 'Bacteria[Orgn] AND '+'genome AND RefSeq'
        temprec=esearch_database_for_querystring(gene, addstr, databasechoice, 1000000)
        allids.append(temprec['IdList'])
        ####Define the dictionary
        for UID in temprec['IdList']:
            if UID not in ncbi_genomes_containing_PHA_genes_dict:
                ncbi_genomes_containing_PHA_genes_dict[UID]={}
            #### give positive returns a value on 1, then the df can be filled with 0s as negative values at a later date
            if gene not in ncbi_genomes_containing_PHA_genes_dict[UID]:
                ncbi_genomes_containing_PHA_genes_dict[UID][str(gene)]=1
        loopcount+=1
        print("Found " + str(temprec['Count']) + ' ' + str(gene))
        time.sleep(2)


       
master_PHA_genome_df=pd.DataFrame.from_dict(ncbi_genomes_containing_PHA_genes_dict)
#master_PHA_genome_df.to_csv("22gene_Master_PHA_genome_df_from_nuccore_RefSeqs.csv", sep=',')
#master_PHA_genome_df=pd.read_csv("22gene_Master_PHA_genome_df_from_nuccore_RefSeqs.csv", sep=',', index_col=0, dtype='unicode')
#master_PHA_genome_df=pd.read_csv("Master_PHA_genome_df_from_genomes.csv", sep=',', index_col=0, dtype='unicode')


master_PHA_genome_df=master_PHA_genome_df.T
master_PHA_genome_df_C_C2only=master_PHA_genome_df
master_PHA_genome_df_C_C2only=master_PHA_genome_df_C_C2only.fillna(0)
master_PHA_genome_df_C_C2only['sum']=master_PHA_genome_df_C_C2only[pha_gene_list].astype(int).sum(axis=1)
master_PHA_genome_df_C_C2only[pha_gene_list]=master_PHA_genome_df_C_C2only[pha_gene_list].astype(float)
##filter the dataframe
##filter to contain only this with PhaC OR PhaC2
master_PHA_genome_df_C_C2only=master_PHA_genome_df_C_C2only.query('PhaC>0|PhaC2>0')
#master_PHA_genome_df_C_C2only.to_csv("22gene_Master_PHA_genome_df_from_nuccore_RefSeqs_C_C2.csv", sep=',')    
#master_PHA_genome_df_C_C2only=pd.read_csv("22gene_Master_PHA_genome_df_from_nuccore_RefSeqs_C_C2.csv", sep=',', index_col=0, dtype='unicode')

#################################################################
#### EXTRACT BIOSAMPLE AND BIOPROJECT FROM GENOME NUCCORE UID ###    
#################################################################
def split_query_into_200_dict(inlist, splitby):
###this is nessessary as entrex will only accept 200 UIDs/accessions in one go\
### otherwise it times out randomly missing some queries    
    axn_split_dict={}
    index_tuples_check=[]
    maxval=math.ceil(len(inlist)/splitby)-1
    for index in list(range(0,maxval+1)):
        if not index==maxval:
            axn_split_dict[index]=inlist[index*splitby:(index+1)*splitby]
            index_tuples_check.append((index, index+1,index*splitby,(index+1)*splitby))
        elif index==maxval:
            axn_split_dict[index]=inlist[index*splitby:]
            index_tuples_check.append((index, index+1,index*splitby,'end'))
    return axn_split_dict

PhaC_genomes_200_dict=split_query_into_200_dict(list(master_PHA_genome_df_C_C2only.index), 200)
###check where the genome is from, efetch -> biosample/bioproject, biosample->isolation source#
sourcedict={}
errorids=[]
superfails=0
num_retries=4
recsumcum=0
for idx in PhaC_genomes_200_dict:
    #remove NAs
    queryidlist=[x for x in PhaC_genomes_200_dict[idx] if not x=='NA']
    lenquery=len(queryidlist)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for UID in PhaC_genomes_200_dict[idx]:
         sourcedict[UID]={}
    for retrynum in list(range(0,num_retries)):
             try:
                 print('Fetching biosample...')
                 handleBS = Entrez.elink(dbfrom="nuccore", id=queryidlist, db='biosample')
                 time.sleep(1)
                 print('Fetching bioproject...')
                 handleBR = Entrez.elink(dbfrom="nuccore", id=queryidlist, db='bioproject')
                 print('Reading...')
                 recordsBS=Entrez.read(handleBS)
                 recordsBR=Entrez.read(handleBR)
                 print('Parsing...')
                 for recnum in list(range(0, len(recordsBS))):
                       if not recordsBS[recnum]['LinkSetDb'] == []:
                            sourcedict[float(recordsBS[recnum]['IdList'][0])]['biosampleUID']=recordsBS[recnum]['LinkSetDb'][0]['Link'][0]['Id']
                       if not recordsBR[recnum]['LinkSetDb'] == []:
                            sourcedict[float(recordsBR[recnum]['IdList'][0])]['bioprojectUID']=recordsBR[recnum]['LinkSetDb'][0]['Link'][0]['Id']
                       recsumcum=recsumcum+1
                 if not len(recordsBS)==lenquery or not len(recordsBR)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordsum)) + '. SHOULD BE ' +str(lenquery) + '!')
                 time.sleep(2)  
                 print("Total records archived: " + str(recsumcum))
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions") 
#####map the biosample UID and bioproject UIDs onto the main dataframe
bioSdict={}
bioRdict={}
for UID in sourcedict:
    if 'biosampleUID' in sourcedict[UID]:
         bioSdict[UID]=sourcedict[UID]['biosampleUID']
    if 'bioprojectUID' in sourcedict[UID]:
         bioRdict[UID]=sourcedict[UID]['bioprojectUID']

master_PHA_genome_df_C_C2only['Biosample_UID']=master_PHA_genome_df_C_C2only.index.map(bioSdict)
master_PHA_genome_df_C_C2only['Bioproject_UID']=master_PHA_genome_df_C_C2only.index.map(bioRdict)
#master_PHA_genome_df_C_C2only.to_csv("22gene_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_biosample.csv", sep=',')
#master_PHA_genome_df_C_C2only=pd.read_csv("22gene_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_biosample.csv", sep=',', index_col=0, dtype='unicode')

#################################################################
#### EXTRACT BIOSAMPLE AND BIOPROJECT FROM GENOME NUCCORE UID ###   
### PART 2 -> LINK KEYWORDS TO BIOSAMPLE ISOLATION SOURCE     ###
#################################################################
BioSampleDict=split_query_into_200_dict(list(master_PHA_genome_df_C_C2only['Biosample_UID']), 200)
keyword_df=pd.read_excel('Sun et al supplementary material/Table_1_Phylogenetic Distribution of Polysaccharide-Degrading Enzymes in Marine Bacteria.XLSX', sheet_name=0, header=0,engine="openpyxl") 
keyword_list=list(keyword_df['Word list'])

keyworddict={}
orgdict={}
errorids=[]
superfails=0
num_retries=4
recsumcum=0
for idx in BioSampleDict:
    #remove NAs
    queryidlist=[x for x in BioSampleDict[idx] if not x=='NA']
    lenquery=len(queryidlist)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for UID in BioSampleDict[idx]:
         keyworddict[UID]={}
         orgdict[UID]={}
    for retrynum in list(range(0,num_retries)):
             try:
                 handleBS=Entrez.efetch(db='biosample', id=queryidlist, rettype='docsum', retmode='xml')
                 recordBS=Entrez.read(handleBS)
                 for recnum in list(range(0, len(recordBS['DocumentSummarySet']['DocumentSummary']))):
                      kwlist=[]
                      if 'display_name="isolation source"' in str(recordBS['DocumentSummarySet']['DocumentSummary'][recnum]['SampleData']):
                             ret1=recordBS['DocumentSummarySet']['DocumentSummary'][recnum]['SampleData'].split('display_name="isolation source"')[1]
                             if not ret1.split(" ")[0]=='>missing</Attribute>':
                                isolation_source_string=ret1.split(" ")[0].replace("Attribute", '').strip("<>/")
                                for keyword in keyword_list:
                                      if keyword in isolation_source_string:
                                          kwlist.append(keyword)
                      uid=recordBS['DocumentSummarySet']['DocumentSummary'][recnum].attributes['uid']
                      org=recordBS['DocumentSummarySet']['DocumentSummary'][recnum]['Organism']
                      keyworddict[uid]=kwlist
                      orgdict[uid]=org
                      recsumcum=recsumcum+1
                 if not len(recordBS['DocumentSummarySet']['DocumentSummary'])==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordBS['DocumentSummarySet']['DocumentSummary'])) + '. SHOULD BE ' +str(lenquery) + '!')
                 time.sleep(2)  
                 print("Total records archived: " + str(recsumcum))
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions") 

master_PHA_genome_df_C_C2only['Biosample_keywords']=master_PHA_genome_df_C_C2only.Biosample_UID.map(keyworddict)
master_PHA_genome_df_C_C2only['Organism']=master_PHA_genome_df_C_C2only.Biosample_UID.map(orgdict)

###read and write the dataframe of filtering will not occur correctly
#master_PHA_genome_df_C_C2only.to_csv("22gene_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_biosample_KW.csv", sep=',')  
#master_PHA_genome_df_C_C2only=pd.read_csv("22gene_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_biosample_KW.csv", sep=',', index_col=0, dtype='unicode')

#arc_PHA_master_C_C2_only.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_version_biosample_KW.csv", sep=',')
#arc_PHA_master_C_C2_only=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_version_biosample_KW.csv", sep=',', index_col=0, dtype='unicode')

##remove non-matching keywords
master_PHA_genome_df_C_C2only_kwmatches=master_PHA_genome_df_C_C2only[master_PHA_genome_df_C_C2only['Biosample_keywords'].str.len() > 2]

#master_PHA_genome_df_C_C2only_kwmatches.to_csv("22gene_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered.csv", sep=',')  
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22gene_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered.csv", sep=',', index_col=0, dtype='unicode')

#arc_PHA_master_C_C2_only_kwmatches.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_version_biosample_KW_filtered.csv", sep=',')
#arc_PHA_master_C_C2_only_kwmatches=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_version_biosample_KW_filtered.csv", sep=',', index_col=0, dtype='unicode')

#c1=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered.csv", sep=',', index_col=0, dtype='unicode')
#c2=pd.read_csv("Master_PHA_genome_df_from_nuccore_new_C_C2_RefSeqs_version_biosample_KW_filtered.csv", sep=',', index_col=0, dtype='unicode')




t=sns.heatmap(master_PHA_genome_df_C_C2only_kwmatches[pha_gene_list].astype(float))
            
#################################################################
####      Count number of proteins in genome                  ###   
###                                                           ###
################################################################# 
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv('Master_PHA_genome_df_from_nuccore_new_C_C2_Keyword_matches.csv', sep=',', index_col=0,dtype='unicode')
genomeiddict=split_query_into_200_dict(list(master_PHA_genome_df_C_C2only_kwmatches.index), 100)  
   
protcountdict={}
errorids=[]
superfails=0
num_retries=4
recsumcum=0
genomefetchcount=0
for idx in genomeiddict:
    #remove NAs
    queryidlist=[x for x in genomeiddict[idx] if not x=='NA']
    lenquery=len(queryidlist)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for UID in genomeiddict[idx]:
         protcountdict[UID]=0
    for retrynum in list(range(0,num_retries)):
             try:
                 handle=Entrez.elink(dbfrom="nuccore", id=queryidlist, db='protein')
                 record=Entrez.read(handle)
                 time.sleep(0.5)
                 genomefetchcount+=lenquery
                 print('Returned proteins for ' +str(genomefetchcount) + ' genomes: index ' + str(idx)+ ' of ' + str(len(genomeiddict))+ '.')
                 if not len(record)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordsum)) + '. SHOULD BE ' +str(lenquery) + '!')
                 for recordnum in list(range(0, len(record))):
                     if not record[recordnum]['LinkSetDb']==[]:
                          UID=int(record[recordnum]['IdList'][0])
                          if UID not in queryidlist:
                              print('unknown uid')
                          proteinlist=[list(x.values())[0] for x in record[recordnum]['LinkSetDb'][0]['Link']]
                          protcount=len(proteinlist)
                          protcountdict[UID]=protcount
                     else:
                         protcountdict[UID]=0
                 if not len(record)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(record)) + '. SHOULD BE ' +str(lenquery) + '!')
                 time.sleep(2)  
                 print("Total records archived: " + str(recsumcum))
             except Exception:
                 print('Initial fail with index: ' + str(idx) + ': UID - ' + str (UID) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
             genomefetchcount+=lenquery
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')           
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200)) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions")     
    
master_PHA_genome_df_C_C2only_kwmatches['Number CDS']=master_PHA_genome_df_C_C2only_kwmatches.index.map(protcountdict)    
#master_PHA_genome_df_C_C2only_kwmatches.to_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count.csv", sep=',')  
#################################################################
####      Retrieve assembly UID for raw genome reads          ###   
###                                                           ###
#################################################################    
    
asemblydict={}
errorids=[]
superfails=0
num_retries=4
recsumcum=0
genomefetchcount=0
for idx in genomeiddict:
    #remove NAs
    queryidlist=[x for x in genomeiddict[idx] if not x=='NA']
    lenquery=len(queryidlist)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered!
    for UID in genomeiddict[idx]:
         asemblydict[UID]=0
    for retrynum in list(range(0,num_retries)):
             try:
                 handle=Entrez.elink(dbfrom="nuccore", id=queryidlist, db='assembly')
                 record=Entrez.read(handle)
                 time.sleep(0.5)
                 genomefetchcount+=lenquery
                 print('Returned proteins for ' +str(genomefetchcount) + ' genomes: index ' + str(idx)+ ' of ' + str(len(genomeiddict))+ '.')
                 if not len(record)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(recordsum)) + '. SHOULD BE ' +str(lenquery) + '!')
                 for recordnum in list(range(0, len(record))):
                     if not record[recordnum]['LinkSetDb']==[]:
                          UID=int(record[recordnum]['IdList'][0])
                          if UID not in queryidlist:
                              print('unknown uid')
                          asemblydict[UID]=record[recordnum]['LinkSetDb'][0]['Link'][0]['Id']
                     else:
                         asemblydict[UID]=0
                 if not len(record)==lenquery:
                     print('Error with index: ' + str(idx) + " RETURNED ONLY " + \
                           str(len(record)) + '. SHOULD BE ' +str(lenquery) + '!')
                 time.sleep(2)  
                 print("Total records archived: " + str(recsumcum))
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else: 
                 if retrynum>0:
                     print('Success!')
                 break    
             genomefetchcount+=lenquery
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')           
    print("Iterated over: "  +str(idx*len(genomeiddict[idx])) +':' + str((idx+1)*len(genomeiddict[idx]))) 
    time.sleep(1)
print("FINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions")     
    
master_PHA_genome_df_C_C2only_kwmatches['Assembly UID']=master_PHA_genome_df_C_C2only_kwmatches.index.map(asemblydict)   
 
#master_PHA_genome_df_C_C2only_kwmatches.to_csv("22gene_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly.csv", sep=',')  
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22gene_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly.csv", sep=',', index_col=0, dtype='unicode')
 

#################################################################
####             Retreive Assembly data for QC                ###   
###                                                           ###
################################################################# 
##to get assemblies, we must find the ftp link using esummary
###modified after https://dmnfarrell.github.io/bioinformatics/assemblies-genbank-python
import gzip, shutil
import pandas as pd
from Bio import Entrez
import urllib
import os
import json
import time
Entrez.email='daniel.leadbeater@york.ac.uk'
ftpdict={}
N50dict={}
assemblykeydict={}
num_retries=4
completed_ids=[]
bad_ftpdict={}
bad_ftps=[]
genomecount=0
###change to true to download
download=False
final_folder_name='Downloaded_assembly_fastas'
master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly.csv", sep=',', index_col=0, dtype='unicode')
idx=0
idxf=0
idxn=0
for UID in list(master_PHA_genome_df_C_C2only_kwmatches['Assembly UID']):
    ###make folder to dump fastas in
    cwd=os.getcwd()
    #check folders exist to in case script needs re-running
    if not os.path.exists(os.path.join(cwd,final_folder_name)):
        os.makedirs(final_folder_name)           
    if UID not in completed_ids:
        for retrynum in list(range(0,num_retries)):
            try:
                   #get ftp using full report from esummary
                   esummary_handle = Entrez.esummary(db="assembly", id=UID, report="full")
                   esummary_record = Entrez.read(esummary_handle)
                   ftpurl = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
                   N50=esummary_record['DocumentSummarySet']['DocumentSummary'][0]['ContigN50']
                   print('Doing N50s : ' + str(idxn))
                   idxn+=1
                   if ftpurl == '':
                       idxf+=1
                       if UID not in bad_ftpdict:
                           bad_ftpdict[UID]='Fail'
                           bad_ftps.append(UID)
                           print('\nBad URL - UID::' + str(UID))
                           idxf+=1
                       continue  
                   label = os.path.basename(ftpurl)
                   ftpdict[UID]=ftpurl
                   assemblykeydict[UID]=label
                   N50dict[UID]=N50
                   link = os.path.join(ftpurl+'/'+label+'_genomic.fna.gz')
                   ##attempt to download the .gz assembly fasta
                   if download==False:
                       break
                   if download==True:
                       urllib.request.urlretrieve(link, f'{label}.fna.gz')
                   print("Downloaded: " + str(label) + ' unzipping...')
                   ##unzip the fasta and store in a single directory
                   with gzip.open(os.path.join(cwd, label +'.fna.gz'), 'r') as f_in, open(os.path.join(cwd,final_folder_name, label).strip('.gz')+'.fna', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                   print('Unzipped ' + str(label) + '.fna')
                   completed_ids.append(UID)
                   print ('\nCompleted ' + str(idx) + ' genome (index ' +str(idx) +') downloads of ' + str(len(list(master_PHA_genome_df_C_C2only_kwmatches['Assembly UID']))))
                   print('Failed downloadeds so far = ' + str(len(bad_ftps)) + '. Empty urls: '+str(idxf))
                   time.sleep(0.5)
            except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
            else: 
                 if retrynum>0:
                     print('Success!')
                 break
    idx+=1
print('Bad URLS dict:\n\n')
print(bad_ftpdict)

print('\n\nBad FTP list:\n\n')
print(bad_ftps)

with open('badURLsdict.json', 'w') as f:
    json.dump(bad_ftpdict, f)
    
with open('ftplinkdict.json', 'w') as f:
    json.dump(ftpdict, f)    

with open('assemblykeydict.json', 'w') as f:
    json.dump(assemblykeydict, f)              

with open('N50dict.json', 'w') as f:
    json.dump(N50dict, f)  

# with open('22g_N50dict.json', 'r') as fp:
#     N50dict = json.load(fp)
    
# with open('22g_ftplinkdict.json', 'r') as fp:
#     ftpdict = json.load(fp)
    
# with open('22g_assemblykeydict.json', 'r') as fp:
#     assemblykeydict = json.load(fp)    
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2.csv", sep=',', index_col=0, dtype='unicode')

master_PHA_genome_df_C_C2only_kwmatches['FTPlink']=master_PHA_genome_df_C_C2only_kwmatches['Assembly UID'].map(ftpdict)
master_PHA_genome_df_C_C2only_kwmatches['Assembly_accession']=master_PHA_genome_df_C_C2only_kwmatches['Assembly UID'].map(assemblykeydict)
master_PHA_genome_df_C_C2only_kwmatches['N50']=master_PHA_genome_df_C_C2only_kwmatches['Assembly UID'].map(N50dict)
#master_PHA_genome_df_C_C2only_kwmatches.to_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F1.csv", sep=',')    
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F1.csv", sep=',', index_col=0, dtype='unicode')

#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50.csv", sep=',', index_col=0, dtype='unicode')
#master_PHA_genome_df_C_C2only_kwmatches.to_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50.csv", sep=',')

filtered_df.Accession.map(id_data_dict)
#####Import and read data output from CheckM QA###
qa_output='qa_output.tsv'
qa_out_df=pd.read_csv(qa_output, sep='\t', index_col=0, dtype='unicode')

master_PHA_genome_df_C_C2only_kwmatches['CheckM_completeness']=master_PHA_genome_df_C_C2only_kwmatches['Assembly_accession'].map(qa_out_df['Completeness'])
master_PHA_genome_df_C_C2only_kwmatches['CheckM_contamination']=master_PHA_genome_df_C_C2only_kwmatches['Assembly_accession'].map(qa_out_df['Contamination'])

arc_qa_output='Analysis_for_paper/arc_qa_output.tsv'
arc_qa_out_df=pd.read_csv(arc_qa_output, sep='\t', index_col=0, dtype='unicode')

arc_PHA_master_C_C2_only_kwmatches['CheckM_completeness']=arc_PHA_master_C_C2_only_kwmatches['Assembly_accession'].map(arc_qa_out_df['Completeness'])
arc_PHA_master_C_C2_only_kwmatches['CheckM_contamination']=arc_PHA_master_C_C2_only_kwmatches['Assembly_accession'].map(arc_qa_out_df['Contamination'])

#master_PHA_genome_df_C_C2only_kwmatches.to_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M.csv", sep=',') 
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M.csv", sep=',', index_col=0, dtype='unicode')
 
#arc_PHA_master_C_C2_only_kwmatches.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F2.csv", sep=',') 
#arc_PHA_master_C_C2_only_kwmatches=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F2.csv", sep=',', index_col=0, dtype='unicode')

#################################################################
####   Retreive genome fastas for Amphora and annotation      ###   
###                                                           ###
################################################################# 
import gzip, shutil
import pandas as pd
from Bio import Entrez, SeqIO
import urllib
import os
import json
import time
Entrez.email='daniel.leadbeater@york.ac.uk'
UID_gb_axn_dict={}
num_retries=4
completed_ids=[]
genomecount=0
final_folder_name='arc_Downloaded_genome_fastas'
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly.csv", sep=',', index_col=0, dtype='unicode')
idx=0
errorids=[]
#check folders exist to in case script needs re-running  
if not os.path.exists(os.path.join(cwd,final_folder_name)):
    os.makedirs(final_folder_name)  
for UID in list(arc_PHA_master_C_C2_only_kwmatches.index):
    ###make folder to dump fastas in
    cwd=os.getcwd()  
    for retrynum in list(range(0,num_retries)):
            try:
                   #get ftp using full report from esummary
                   print('Fetching...'+str(UID))
                   handlefetch=Entrez.efetch(db="nuccore", id=UID, rettype='fasta', retmode='text')
                   print('Fetched...')
                   genome=SeqIO.read(handlefetch, 'fasta')
                   print('Genome length = ' + str(len(genome.seq)))
                   f=open(os.path.join(final_folder_name, str(genome.id)+'.fasta'),'w')
                   print('Writing to... ' + str(os.path.join(final_folder_name, str(genome.id)+'.fasta')))
                   SeqIO.write(genome, f, 'fasta')
                   print('Written... appending...')
                   f.close()
                   UID_gb_axn_dict[UID]=str(genome.id)
                   idx+=1
                   print("Downloaded " + str(genome.id) + " index: " + str(idx) + ' of ' +str(len(list(arc_PHA_master_C_C2_only_kwmatches.index))) + '\n')
                   time.sleep(1)
            except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum) + ' Trying again in ' + str(2+(2*retrynum)) +'s')
                 time.sleep(2+(2*retrynum))
                 if retrynum==num_retries:
                     print('FAILURE AFTdos2unix DoER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                     errorids.append(idx)
                 time.sleep(2)
            else: 
                 if retrynum>0:
                     print('Success!')
                 break
    idx+=1
print('Full errors: ' + str(errorids))

# with open('arc_UID_gb_to_axn.json', 'w') as f:
#     json.dump(UID_gb_axn_dict, f)

with open('bac_UID_gb_to_axn.json', 'r') as fp:
    bac_UID_gb_to_axn = json.load(fp)
arc_PHA_master_C_C2_only_kwmatches['Genome_accesion']=arc_PHA_master_C_C2_only_kwmatches.index.map(UID_gb_axn_dict)

master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion']=master_PHA_genome_df_C_C2only_kwmatches.index.map({int(k): v for (k, v) in bac_UID_gb_to_axn.items()})

#master_PHA_genome_df_C_C2only_kwmatches.to_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F3.csv", sep=',') 
#arc_PHA_master_C_C2_only_kwmatches.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F3.csv", sep=',') 


#################################################################
####                     Check AMPHORA OUTPUT                 ###   
###                                                           ###
#################################################################    
###create dict
marker_gene_count_dict={}
for key in list(arc_PHA_master_C_C2_only_kwmatches['Genome_accesion']):
    marker_gene_count_dict[key]={}
    marker_gene_count_dict[key]['count']=0
    marker_gene_count_dict[key]['markers']=[]

amphora_output_directory='../amphora_out_arc'
for pepfile in os.listdir(amphora_output_directory):
        marker=str(pepfile)[:-4]
        for seq in SeqIO.parse(os.path.join(amphora_output_directory,pepfile), 'fasta'):
            ##isolate the genome accession
            if str(seq.description.split(' ')[0].rsplit('_',1)[0]) in marker_gene_count_dict:
                print('1')
                marker_gene_count_dict[str(seq.description.split(' ')[0].rsplit('_',1)[0])]['count']=marker_gene_count_dict[str(seq.description.split(' ')[0].rsplit('_',1)[0])]['count']+1
                marker_gene_count_dict[str(seq.description.split(' ')[0].rsplit('_',1)[0])]['markers'].append(marker)

count_dict_map={}
markers_dict_map={}
unique_marker_count_map={}
for axn in marker_gene_count_dict:
    count_dict_map[axn]=marker_gene_count_dict[axn]['count']
    markers_dict_map[axn]=marker_gene_count_dict[axn]['markers']
    unique_marker_count_map[axn]=len(set(marker_gene_count_dict[axn]['markers']))

master_PHA_genome_df_C_C2only_kwmatches['Amphora_total_markers']=master_PHA_genome_df_C_C2only_kwmatches.Genome_accesion.map(count_dict_map)
master_PHA_genome_df_C_C2only_kwmatches['Amphora_markers']=master_PHA_genome_df_C_C2only_kwmatches.Genome_accesion.map(markers_dict_map)
master_PHA_genome_df_C_C2only_kwmatches['Amphora_unique_markers']=master_PHA_genome_df_C_C2only_kwmatches.Genome_accesion.map(unique_marker_count_map)

#master_PHA_genome_df_C_C2only_kwmatches.to_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F4.csv", sep=',') 

arc_PHA_master_C_C2_only_kwmatches['Amphora_total_markers']=arc_PHA_master_C_C2_only_kwmatches.Genome_accesion.map(count_dict_map)
arc_PHA_master_C_C2_only_kwmatches['Amphora_markers']=arc_PHA_master_C_C2_only_kwmatches.Genome_accesion.map(markers_dict_map)
arc_PHA_master_C_C2_only_kwmatches['Amphora_unique_markers']=arc_PHA_master_C_C2_only_kwmatches.Genome_accesion.map(unique_marker_count_map)

#arc_PHA_master_C_C2_only_kwmatches.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F4.csv", sep=',') 

#################################################################
####                     GET TAXONOMIC INFO                   ###   
###                                                           ###
################################################################# 
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F4.csv", sep=',', index_col=0, dtype='unicode')
#arc_PHA_master_C_C2_only_kwmatches=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F4.csv", sep=',', index_col=0, dtype='unicode')

from ete3 import NCBITaxa
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()


desired_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
def get_desired_ranks_output_dict_from_existing_dict(taxdict, desired_ranks):
    no_ids=[]
    for UID in taxdict:
        taxid=taxdict[UID]['taxid']
        try:
            lineage = ncbi.get_lineage(taxid)
            lineage2ranks = ncbi.get_rank(lineage)
            lineagetranslated=ncbi.get_taxid_translator(lineage)
            ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
            for x in desired_ranks:
                 if x in list(ranks2lineage.keys()):
                       taxdict[UID][x]={}
                 else:
                       taxdict[UID][x]='NA' 
            for taxrank, taxid in ranks2lineage.items():
                  if taxrank in desired_ranks:
                       taxdict[UID][taxrank]=lineagetranslated[taxid]
                  else:
                      taxdict[UID][taxrank]='NA'
            taxdict[UID]['lineage']=lineagetranslated
        except:
            #print('Failed to gather taxid for id: ' + str(UID))
            for x in desired_ranks:
                taxdict[UID][x]='NA'
            taxdict[UID]['lineage']='NA'
            no_ids.append(UID)
    print("Number of taxids not retrievable: " + str(len(no_ids)))
    return taxdict

##check over missed taxids
#t=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches.Genome_accesion.isin(validated_pha_count_df_filt.index)]
taxid_dict={}
for idx in list(range(0, int(len(master_PHA_genome_df_C_C2only_kwmatches.index)/1000)+1)):
        if idx==0:
             query_list=list(master_PHA_genome_df_C_C2only_kwmatches.index)[idx:(idx+1)*1000]
             print("Fetching handle")
             handle=Entrez.elink(dbfrom="nuccore", id=query_list, db='taxonomy')
             print("Returned handle")
             record=Entrez.read(handle)
             print("Read handle")
             for Z in list(range(0, len(record))):
                 if int(Z)%200==0:
                      print('Iterated over ' + str(Z) + ' indexes to get taxids')
                 tid=int(record[Z]['IdList'][0])
                 taxid_dict[tid]={}
                 if not record[Z]['LinkSetDb'] ==[]:
                     taxid_dict[tid]['taxid']=record[Z]['LinkSetDb'][0]['Link'][0]['Id']
                 else:
                     taxid_dict[tid]['taxid']='NA'
             print("Iterated over " + str((idx+1)*1000) + " taxids")
        if idx==list(range(0, int(len(master_PHA_genome_df_C_C2only_kwmatches.index)/1000)+1))[-1]:
             query_list=list(master_PHA_genome_df_C_C2only_kwmatches.index)[idx*1000:]
             print("Fetching handle")
             handle=Entrez.elink(dbfrom="nuccore", id=query_list, db='taxonomy')
             print("Returned handle")
             record=Entrez.read(handle)
             print("Read handle")
             for Z in list(range(0, len(record))):
                 if int(Z)%200==0:
                      print('Iterated over ' + str(Z) + ' indexes to get taxids')
                 tid=int(record[Z]['IdList'][0])
                 taxid_dict[tid]={}
                 if not record[Z]['LinkSetDb'] ==[]:
                     taxid_dict[tid]['taxid']=record[Z]['LinkSetDb'][0]['Link'][0]['Id']
                 else:
                     taxid_dict[tid]['taxid']='NA'     
             print("Iterated over " + str(len(master_PHA_genome_df_C_C2only_kwmatches.index)-(idx*1000)) + " taxids")                     
        else:
             query_list=list(master_PHA_genome_df_C_C2only_kwmatches.index)[idx*1000:(idx+1)*1000]
             print("Fetching handle")
             handle=Entrez.elink(dbfrom="nuccore", id=query_list, db='taxonomy')
             print("Returned handle")
             record=Entrez.read(handle)
             print("Read handle")
             for Z in list(range(0, len(record))):
                 if int(Z)%200==0:
                      print('Iterated over ' + str(Z) + ' indexes to get taxids')
                 tid=int(record[Z]['IdList'][0])
                 taxid_dict[tid]={}
                 if not record[Z]['LinkSetDb'] ==[]:
                     taxid_dict[tid]['taxid']=record[Z]['LinkSetDb'][0]['Link'][0]['Id']
                 else:
                     taxid_dict[tid]['taxid']='NA'             
             print("Iterated over " + str((idx+1)*1000) + " taxids")
            
taxid_dict=get_desired_ranks_output_dict_from_existing_dict(taxid_dict, desired_ranks)

taxdf=pd.DataFrame.from_dict(taxid_dict, orient='index')
master_PHA_genome_df_C_C2only_kwmatches[['phylum','class','order','family','genus','species','lineage']]=taxdf[['phylum','class','order','family','genus','species','lineage']]


#master_PHA_genome_df_C_C2only_kwmatches.to_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax.csv", sep=',') 
#arc_PHA_master_C_C2_only_kwmatches.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F5.csv", sep=',') 
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F5.csv", sep=',', index_col=0, dtype='unicode')
#arc_PHA_master_C_C2_only_kwmatches=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F5.csv", sep=',', index_col=0, dtype='unicode')


#################################################################
####                     Quality_filtering                    ###   
###                                                           ###
################################################################# 
master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F6.csv", sep=',', index_col=0, dtype='unicode')
arc_PHA_master_C_C2_only_kwmatches=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F5.csv", sep=',', index_col=0, dtype='unicode')

master_PHA_genome_df_C_C2only_kwmatches['N50']=master_PHA_genome_df_C_C2only_kwmatches['N50'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['CheckM_completeness']=master_PHA_genome_df_C_C2only_kwmatches['CheckM_completeness'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['CheckM_contamination']=master_PHA_genome_df_C_C2only_kwmatches['CheckM_contamination'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['Amphora_unique_markers']=master_PHA_genome_df_C_C2only_kwmatches['Amphora_unique_markers'].astype(float,errors='ignore')


master_PHA_genome_df_F5_F=master_PHA_genome_df_C_C2only_kwmatches.query('N50>50000&CheckM_completeness>=95&CheckM_contamination<=5&Amphora_unique_markers==31')

arc_PHA_master_C_C2_only_kwmatches['N50']=arc_PHA_master_C_C2_only_kwmatches['N50'].astype(float,errors='ignore')
arc_PHA_master_C_C2_only_kwmatches['CheckM_completeness']=arc_PHA_master_C_C2_only_kwmatches['CheckM_completeness'].astype(float,errors='ignore')
arc_PHA_master_C_C2_only_kwmatches['CheckM_contamination']=arc_PHA_master_C_C2_only_kwmatches['CheckM_contamination'].astype(float,errors='ignore')
arc_PHA_master_C_C2_only_kwmatches['Amphora_unique_markers']=arc_PHA_master_C_C2_only_kwmatches['Amphora_unique_markers'].astype(float,errors='ignore')


arc_master_PHA_genome_df_F5_F=arc_PHA_master_C_C2_only_kwmatches.query('N50>50000&CheckM_completeness>=95&CheckM_contamination<=5&Amphora_unique_markers>=90')


#################################################################
####                     VALIDATE PHAC GENES                  ###   
###                                                           ###
################################################################# 
#master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv('Master_PHA_genome_df_from_nuccore_new_C_C2_Keyword_matches.csv', sep=',', index_col=0,dtype='unicode')
idqueseries=list(master_PHA_genome_df_F5_F.index)
PHA_check_dict={}
PHA_info_dict={}
PhaC_errors={}
errorids=[]
superfails=0
num_retries=4
gencount=0
genomeparsecount=0
GOI=['phaa', 'phac', 'phae', 'phag', 'phad', 'phaf', 'fabg', 'phar', 'phaz', 'phac2', 'maoc', 'phab',\
     'phap', 'phaj', 'phaq', 'phai', 'pct', 'pcs', 'prpe', 'pdup', 'cat2', 'hada']
genome_fasta_folder='High_quality_genome_fastas'
genomes_only=True
if genome_fasta_folder not in os.listdir():
    os.mkdir(genome_fasta_folder)
for idx in idqueseries:
    #remove NAs
    #esummaryquery=[x for x in genomeiddict[idx] if not x=='NA']#
    esummaryquery=idqueseries
    lenquery=len(esummaryquery)
    #esummaryquery=','.join(esummaryquery)
    ###NB esummary must have the query string formatted like this
    ### however the elink query must must be a list of strings or it returns info disordered
    PHA_check_dict[idx]={}
    PHA_info_dict[idx]={}
    PhaC_errors[idx]={}
    for retrynum in list(range(0,num_retries)):
             if genomes_only==True:
                 rettypestr='fasta'
                 retmodestr='text'
                 filetype='.fasta'
                 readas='fasta'
             if not genomes_only==True:
                rettypestr='gbwithparts'
                retmodestr='text'
                filetype='.genabnk'
                reasas='gb'
             try:
                 fn=os.path.join(os.getcwd(),genome_fasta_folder,str(idx)+filetype)
                 if str(idx)+'.genbank' not in os.listdir(os.path.join(os.getcwd(),genome_fasta_folder)):
                      print('Fetching genome: ' + str(gencount) + str('. UID: ' + str(idx)))
                      handlefetch=Entrez.efetch(db="nuccore", id=idx, rettype=rettypestr, retmode=retmodestr)   
                      print("Fetched... reading...")
                      t=SeqIO.read(handlefetch, readas)
                      print("Read... writing...")
                      SeqIO.write(t,fn, readas)
                      print('Written... parsing data...')
                      ###quality check the data downloaded and check it's not behind a data wall e.g.
                      ##
                 if genomes_only==True:
                     break
                 for seq in SeqIO.parse(fn, 'gb'):
                     for f in seq.features:
                        #loop over features in feature table for only genes
                        if f.type == "CDS":
                           #if a "gene", check for gene name throught the gene key in the subsequent dict
                           if "gene" in f.qualifiers:
                              if f.qualifiers['gene'][0].lower() in GOI:
                                  #print(f.qualifiers)
                                  gene=f.qualifiers['gene'][0]
                                  locus=f.qualifiers['locus_tag'][0]
                                  description=f.qualifiers['product'][0]
                                  startstop=f.location
                                  #print(gene, description)
                                  if gene+'_count' not in PHA_check_dict[idx]:
                                      PHA_info_dict[idx][gene+'_count']=0
                                      PHA_check_dict[idx][gene]=1
                                      #give the gene a name like PhaC_0, the next will be PhaC_1 etc
                                  gene_num=gene+'_'+str(PHA_info_dict[idx][gene+'_count'])
                                  PHA_info_dict[idx][gene_num]={}
                                  PHA_info_dict[idx][gene_num]['locus']=locus
                                  PHA_info_dict[idx][gene_num]['definition']=description
                                  PHA_info_dict[idx][gene_num]['position']=startstop
                                  PHA_info_dict[idx][gene+'_count']+=1
                                  if gene.lower() == 'phac' and  'multicomponent' in description:
                                      PhaC_errors[idx]='Error with PhaC'                    
                 time.sleep(0.5)
                 print("Total genomes archived: " + str(gencount+1) + ' of: ' + str(len(idqueseries)))
                 gencount+=1
             except Exception:
                 print('Initial fail with index: ' + str(idx) + " Attempt: " + str(retrynum))
                 time.sleep(2)
                 if retrynum==num_retries:
                     print('FAILURE AFTER ' + str(retrynum) + "ATTEMPTS: " +  'index: ' + str(idx))
                 errorids.append(idx)
                 time.sleep(2)
             else:
                 if retrynum>0:
                     print('Success!')
                 break    
    else:
            superfails+=1
            print('ESCAPING FAILURE AFTER ' + str(retrynum) + " ATTEMPTS: " +  'index: ' + str(idx)+'\n\n')          
    print("Iterated over: "  +str(idx*200) +':' + str((idx+1)*200))
    time.sleep(1)
print("\nFINSIHED\nTHERE WERE " + str(superfails) + " FAILURES\n")
if not superfails==0:
    print(" - you MUST go over errorids manually and append the errors - these"\
          "have failed because of servers issues, NOT because of faulty acessions")


PHA_count_only_dict={}
for uid in PHA_info_dict:
    PHA_count_only_dict[uid]={}
    for data in PHA_info_dict[uid].keys():
        if data.endswith('_count'):
            PHA_count_only_dict[uid][data]=PHA_info_dict[uid][data]
            if PHA_info_dict[uid][data]>1:
                print(uid, data)


validated_pha_count_df=pd.DataFrame.from_dict(PHA_count_only_dict).fillna(0).T 
master_PHA_genome_df_F5_F=pd.concat([master_PHA_genome_df_F5_F,validated_pha_count_df], axis=1)



#master_PHA_genome_df_F5_F.to_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F7.csv", sep=',') 
#master_PHA_genome_df_F5_F=pd.read_csv("ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_F7.csv", sep=',', index_col=0, dtype='unicode')

#arc_master_PHA_genome_df_F5_F.to_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F7.csv", sep=',') 
#arc_master_PHA_genome_df_F5_F=pd.read_csv("Archaea_Master_PHA_genome_df_from_nuccore_RefSeqs_F7.csv", sep=',', index_col=0, dtype='unicode')

validated_pha_count_df=master_PHA_genome_df_F5_F[[x for x in list(master_PHA_genome_df_F5_F.columns) if x.endswith('_count')]]
validated_pha_count_df.columns=[x.strip('_count') for x in list(validated_pha_count_df.columns) if x.endswith('_count')]       


def make_clustermap_color_map(taxlevel, dataframe, name):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    num_unique_tax=len(set(list(dataframe[taxlevel])))
    #sns.palplot(sns.color_palette("husl", 12))
    colrs=list(sns.color_palette("Paired",14))+list(sns.color_palette("hls", 12)[::-1])
    c=0
    legend_handles_group=[]
    for x in set(list(dataframe[taxlevel])):
        rowcoldict[x]=colrs[c]
        legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=rowcoldict[x], label=x, lw=1.5))
        c+=1
    ###use this color dict to create a list of colours from the dataframe
    col_list=[]
    for taxa in list(dataframe[taxlevel]):
        col_list.append(rowcoldict[taxa])
    df_colors=pd.Series(col_list, dataframe.index, name=str(name))
    return df_colors, legend_handles_group    

import scipy.spatial as sp, scipy.cluster.hierarchy as hc

cm_row_cols, leg_handles=make_clustermap_color_map('order', master_PHA_genome_df_F5_F, '')

#row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average') for x in (validated_pha_count_df.astype(float).values, validated_pha_count_df.astype(float).values.T))
#k=sns.clustermap(validated_pha_count_df.astype(float), row_linkage=row_linkage, col_linkage=col_linkage, metric="euclidean", row_colors=cm_row_cols, edgecolor='white', linewidth=1,dendrogram_ratio=(0.1,0.1))    
k=sns.clustermap(validated_pha_count_df.astype(float), metric="Euclidean",\
                 row_colors=cm_row_cols, edgecolor='white', linewidth=1,\
                     dendrogram_ratio=(0.08,0.08),yticklabels=False)
plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.setp(k.ax_heatmap.xaxis.get_majorticklabels(), fontsize=16)
k.ax_col_dendrogram
k.cax.set_position([1, .2, .03, .45]) #([x, y, width, height])    
legend2=k.ax_heatmap.legend(loc=(-0.075,-0.15),handles=leg_handles, ncol=4, fontsize=12, facecolor='gray', framealpha=0.1, frameon=True, fancybox=True, edgecolor='black')#title='Taxonomy'
k.ax_heatmap.add_artist(legend2)      
for a in k.ax_row_dendrogram.collections:
    a.set_linewidth(2)
for a in k.ax_col_dendrogram.collections:
    a.set_linewidth(2)      
    

####correlation
corr=validated_pha_count_df.astype(float).corr()
# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(corr, dtype=bool))
f, ax = plt.subplots(figsize=(11, 9))
# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})


#################################################################
####          Parse dbCAN output into dict to json            ###   
###                                                           ###
#################################################################    
    

# import os
# import pandas as pd
# from Bio import SeqIO    
# import re
# import json
   
# #construct a dictionary like this:
# #{genomeuid:{CDS:annotation1,annotation2, start, stop}, CDS2...}
# CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'
# CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)
# CAZy_keyword_list=list(CAZy_familiy_df['CAZy family'])
# #regex_cazy_keys=['^'+x+'$' for x in CAZy_keyword_list]  
# #regex_cazy_str='|'.join(regex_cazy_keys)  
# genome_gene_annot_dict={}
# #construct a dictionary like this:
# #{genomeuid:{CDS:annotation1,annotation2, start, stop}, CDS2...}
# genome_annot_count_dict={}
# ##construct a dictionary like this: for heatmap plots etc
# #genomuid:{cazy:count}
# #path=os.path.join(os.getcwd(), 'dbCAN_output_bacteria') 
# path=os.path.join(os.getcwd(), 'dbCAN_output_archaea')    
# for genomedir in os.listdir(path):
#      genelist=[]
#      genomeUID=str(genomedir).strip('output_')  
#      genome_gene_annot_dict[genomeUID]={}
#      genome_annot_count_dict[genomeUID]={}
#      for uniquecazy in CAZy_keyword_list:
#                   genome_annot_count_dict[genomeUID][uniquecazy]=0
#      pathtooverview=os.path.join(path, genomedir)
#      ##identify gene of interest from the hmmer and hmm outputs
#      ##use the CDS to keep track of where the gene is
#      ###check hmme-csvr output first
#      ovrdf=pd.read_csv(os.path.join(pathtooverview, 'overview.txt'), delimiter='\t')
#      ###determine the total cazy count per genome for heatmaps etc
#      ovrdf['Merged_annots']=ovrdf[['HMMER', 'Hotpep', 'DIAMOND']].agg(','.join, axis=1).str.split('\(|,|\+')
#      for idx in list(range(0, len(ovrdf))):
#          repeats=[]  
#          for y in ovrdf['Merged_annots'].iloc[idx]:
#              if y in CAZy_keyword_list and y not in repeats:
#                  #print(y,ovrdf.iloc[idx]['Gene ID'])
#                  repeats.append(y)
#                  genome_annot_count_dict[genomeUID][y]=genome_annot_count_dict[genomeUID][y]+1
#                  #continue so if hmm, diam and hotpep = the same, count doesn't return 3
#                  continue        
#      hmmer_mask=ovrdf['HMMER'].str.replace('+','(').str.split('(').apply(lambda x: any([k in x for k in CAZy_keyword_list]))
#      diam_mask=ovrdf['DIAMOND'].str.replace('+','(').str.split('(').apply(lambda x: any([k in x for k in CAZy_keyword_list]))
#      hotpep_mask=ovrdf['Hotpep'].str.replace('+','(').str.split('(').apply(lambda x: any([k in x for k in CAZy_keyword_list]))
#      hmmerdf=ovrdf[hmmer_mask]
#      diamdf=ovrdf[diam_mask]
#      hotpepdf=ovrdf[hotpep_mask]
#      for idx in list(range(0, len(hmmerdf))):
#           if not hmmerdf.iloc[idx]['Gene ID'] in genome_gene_annot_dict[genomeUID]:
#               genome_gene_annot_dict[genomeUID][hmmerdf.iloc[idx]['Gene ID']]={}
#               genome_gene_annot_dict[genomeUID][hmmerdf.iloc[idx]['Gene ID']]['CAZy']=[]
#               genelist.append(hmmerdf.iloc[idx]['Gene ID'])
#           for obj in re.split('\(|,|\+', hmmerdf.iloc[idx]['HMMER']):
#               if obj in CAZy_keyword_list and obj not in genome_gene_annot_dict[genomeUID][hmmerdf.iloc[idx]['Gene ID']]['CAZy']:
#                    genome_gene_annot_dict[genomeUID][hmmerdf.iloc[idx]['Gene ID']]['CAZy'].append(obj)
#      for idx in list(range(0, len(diamdf))):                   
#           if not diamdf.iloc[idx]['Gene ID'] in genome_gene_annot_dict[genomeUID]:
#               genome_gene_annot_dict[genomeUID][diamdf.iloc[idx]['Gene ID']]={}
#               genome_gene_annot_dict[genomeUID][diamdf.iloc[idx]['Gene ID']]['CAZy']=[]
#               genelist.append(diamdf.iloc[idx]['Gene ID'])
#           for obj in re.split('\(|,|\+', diamdf.iloc[idx]['DIAMOND']):
#               if obj in CAZy_keyword_list and obj not in genome_gene_annot_dict[genomeUID][diamdf.iloc[idx]['Gene ID']]['CAZy']:
#                    genome_gene_annot_dict[genomeUID][diamdf.iloc[idx]['Gene ID']]['CAZy'].append(obj)
#      for idx in list(range(0, len(hotpepdf))):                      
#           if not hotpepdf.iloc[idx]['Gene ID'] in genome_gene_annot_dict[genomeUID]:
#               genome_gene_annot_dict[genomeUID][hotpepdf.iloc[idx]['Gene ID']]={}
#               genome_gene_annot_dict[genomeUID][hotpepdf.iloc[idx]['Gene ID']]['CAZy']=[]
#               genelist.append(hotpepdf.iloc[idx]['Gene ID'])
#           for obj in re.split('\(|,|\+', hotpepdf.iloc[idx]['Hotpep']):
#               if obj in CAZy_keyword_list and obj not in genome_gene_annot_dict[genomeUID][hotpepdf.iloc[idx]['Gene ID']]['CAZy']:
#                    genome_gene_annot_dict[genomeUID][hotpepdf.iloc[idx]['Gene ID']]['CAZy'].append(obj)                 
#     ##get start and stop locations from the CDS file
#      for seq in SeqIO.parse(os.path.join(pathtooverview,'uniInput'), 'fasta'):
#         if seq.id in genelist:
#             #start=seq.description.split('#')[1]
#             #stop=seq.description.split('#')[2]
#             genome_gene_annot_dict[genomeUID][seq.id]['start']=seq.description.split('#')[1]
#             genome_gene_annot_dict[genomeUID][seq.id]['stop']=seq.description.split('#')[2]    


# with open('CAZy_data_dict_arc.json', 'w') as fp:
#     json.dump(genome_gene_annot_dict, fp)   
    
# with open('Genome_CAZy_count_dict_arc.json', 'w') as fp:
#     json.dump(genome_annot_count_dict, fp)  

#################################################################
####       Parse dbCAN output into dict to json               ###   
###                           END                             ###
#################################################################    
CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)

def make_algae_cazy_clustermap_color_map(input_df):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':'#669900'}
    legend_handles=[]
    legrepeats=[]
    for idx in list(range(0, len(input_df))):
        rowcoldict[input_df.iloc[idx]['CAZy family']]=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()]
        if not CAZy_familiy_df.iloc[idx]['Algae'] in legrepeats:
            legend_handles.append(patches.Patch(color=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()], label=str(CAZy_familiy_df.iloc[idx]['Algae'])))
            legrepeats.append(str(CAZy_familiy_df.iloc[idx]['Algae']))
    return pd.DataFrame.from_dict([rowcoldict]), legend_handles
    ###use this color dict to create a list of colours from the dataframe
    #col_list=[]
   # for xxx in list(dataframe[taxlevel]):
    #     col_list.append(rowcoldict[taxa])
    # df_colors=pd.Series(col_list, dataframe.index, name=str(taxlevel))
    # return df_colors, legend_handles_group    


cazy_cols, cazy_leg_handles=make_algae_cazy_clustermap_color_map(CAZy_familiy_df)


cazy_genome_countdf=pd.read_json('Genome_CAZy_count_dict_bac.json')  
cazy_genome_countdf=cazy_genome_countdf.T  
###get rid of columns that have no data
cazy_genome_countdf_filt=cazy_genome_countdf[cazy_genome_countdf.columns[cazy_genome_countdf.sum(axis=0)>0]]
#t=sns.clustermap(cazy_genome_countdf[cazy_genome_countdf.sum(axis=0)>0], col_colors=cazy_cols.T)   
t=sns.clustermap(cazy_genome_countdf_filt, col_colors=cazy_cols.T, col_cluster=False, cmap='YlGnBu')   
   


#################################################################
####                DATA ANALYSIS STARTS HERE                 ###
####       Merge PHA and dbCAN data into single dataframe     ###
###                           FINAL ANALYSIS                  ###
#################################################################  
master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M.csv", sep=',', index_col=0, dtype='unicode')

master_PHA_genome_df_C_C2only_kwmatches['N50']=master_PHA_genome_df_C_C2only_kwmatches['N50'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['CheckM_completeness']=master_PHA_genome_df_C_C2only_kwmatches['CheckM_completeness'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['CheckM_contamination']=master_PHA_genome_df_C_C2only_kwmatches['CheckM_contamination'].astype(float,errors='ignore')
master_PHA_genome_df_C_C2only_kwmatches['Amphora_unique_markers']=master_PHA_genome_df_C_C2only_kwmatches['Amphora_unique_markers'].astype(float,errors='ignore')


master_PHA_genome_df_F5_F=master_PHA_genome_df_C_C2only_kwmatches.query('N50>50000&CheckM_completeness>=95&CheckM_contamination<=5&Amphora_unique_markers==31')

#master_PHA_genome_df_F5_F.to_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M_filtered.csv", sep=',') 

CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)

cazy_genome_countdf=pd.read_json('Genome_CAZy_count_dict_small.json')  
cazy_genome_countdf=cazy_genome_countdf.T  
###get rid of columns that have no data
cazy_genome_countdf_filt=cazy_genome_countdf[cazy_genome_countdf.columns[cazy_genome_countdf.sum(axis=0)>0]]
#validated_pha_count_df=master_PHA_genome_df_F5_F[[x for x in list(master_PHA_genome_df_F5_F.columns) if x.endswith('_count')]]
validated_pha_count_df=master_PHA_genome_df_F5_F[pha_gene_list]
#validated_pha_count_df.columns=[x.strip('_count') for x in list(validated_pha_count_df.columns) if x.endswith('_count')]
##map genome accession onto pha count dataframe  
validated_pha_count_df['Genome accession']=validated_pha_count_df.index.map(master_PHA_genome_df_F5_F['Genome_accesion'])
###merge the dataframes using the genomaccession from pha df to the index of the cazy df  
###merge is pd.merge(left df, right df, how left, how right)
merged = pd.merge(validated_pha_count_df, cazy_genome_countdf, left_on='Genome accession', right_index=True)
#        .reindex(columns=['id', 'store', 'address', 'warehouse']))  
merged=merged.drop(['Genome accession'], axis=1) 
#filter genomes to contain activity greater than gtac
gtac=1
gtmask=merged[cazy_genome_countdf_filt.columns].sum(axis=1)>gtac
merged=merged[gtmask]
validated_pha_count_df=validated_pha_count_df[gtmask]
   
validated_pha_count_df=validated_pha_count_df.set_index('Genome accession')
###filter out columns with no data
validated_pha_count_df=validated_pha_count_df[validated_pha_count_df.columns[validated_pha_count_df.astype(float).sum(axis=0)>0]]
#k=sns.clustermap(merged.astype(float), col_cluster=False, cmap='YlGnBu')#, col_colors=cazy_cols.T)      
####merged dataframe looks messy, plot the plots separately, but use the dendrogram from plot 1 to make sure the axes
####align between the plots and re order the index




#########
import matplotlib.gridspec  
from matplotlib.transforms import Bbox  
import matplotlib.colors as colors
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import colors

def truncate_colormap(cmap, minval=0.1, maxval=0.9, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def make_algae_cazy_clustermap_color_map(input_df):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':(0.4, 0.6, 0)}
    legend_handles=[]
    legrepeats=[]
    legend_labels=[]
    for idx in list(range(0, len(input_df))):
        rowcoldict[input_df.iloc[idx]['CAZy family']]=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()]
        if not CAZy_familiy_df.iloc[idx]['Algae'] in legrepeats:
            legend_handles.append(patches.Patch(facecolor=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()],\
                                                label=str(CAZy_familiy_df.iloc[idx]['Algae']), edgecolor='black', lw=1))
            legrepeats.append(str(CAZy_familiy_df.iloc[idx]['Algae']))
            legend_labels.append(str(CAZy_familiy_df.iloc[idx]['Algae']))
    return pd.DataFrame.from_dict([rowcoldict]), legend_handles, legend_labels

def make_clustermap_color_map(taxlevel, dataframe, name):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    num_unique_tax=len(set(list(dataframe[taxlevel])))
    #sns.palplot(sns.color_palette("husl", 12))
    colrs=list(sns.color_palette("Paired",14))+list(sns.color_palette("hls", 12)[::-1])
    c=0
    legend_handles_group=[]
    for x in set(list(dataframe[taxlevel])):
        rowcoldict[x]=colrs[c]
        legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=rowcoldict[x], label=x, lw=1))
        c+=1
    ###use this color dict to create a list of colours from the dataframe
    col_list=[]
    col_list_white=[]
    for taxa in list(dataframe[taxlevel]):
        col_list.append(rowcoldict[taxa])
        col_list_white.append((1,1,1,1))
    df_colors=pd.Series(col_list, dataframe.Genome_accesion, name=str(name))
    df_colors_white=pd.Series(col_list_white, dataframe.Genome_accesion, name=str(name))
    return df_colors_white,df_colors, legend_handles_group    

f1cm_row_cols_axn_white,f1cm_row_cols_axn, f1leg_handles=make_clustermap_color_map('order', master_PHA_genome_df_F5_F, '')

f2cazy_colsx, f2cazy_leg_handlesx, f2cazy_leg_labels=make_algae_cazy_clustermap_color_map(CAZy_familiy_df)

#make pha dict grouping genes by function
pha_function_dict={'PhaA':'Synthesis intermediate', 'PhaC':'Synthesis terminal', 'PhaE':'Synthesis terminal',\
                   'PhaG':'Synthesis intermediate', 'PhaD':'Phasin regulator', 'PhaF':'Phasin', 'PhaP':'Phasin',\
                   'PhaZ':'Depolymerase', 'PhaB':'Synthesis intermediate', 'PhaR':'Synthesis regulator', 'maoC':'Synthesis intermediate',\
                   'PhaC2':'Synthesis terminal','PhaI':'Phasin', 'FabG':'Synthesis intermediate',\
                   'PhaJ':'Synthesis intermediate', 'PhaQ':'Phasin regulator', 'PhaP':'Phasin',\
                   'Pcs':'Synthesis intermediate', 'PrpE': 'Synthesis intermediate'}
    
pha_function_dict={'PhaA':'Synthesis intermediate', 'PhaC':'Synthesis terminal', 'PhaE':'Synthesis terminal',\
                   'PhaG':'Synthesis intermediate', 'PhaD':'Synthesis regulator', 'PhaF':'Phasin', 'PhaP':'Phasin',\
                   'PhaZ':'Depolymerase', 'PhaB':'Synthesis intermediate', 'PhaR':'Synthesis regulator', 'maoC':'Synthesis intermediate',\
                   'PhaC2':'Synthesis terminal','PhaI':'Phasin', 'FabG':'Synthesis intermediate',\
                   'PhaJ':'Synthesis intermediate', 'PhaQ':'Synthesis regulator', 'PhaP':'Phasin',\
                   'Pcs':'Synthesis intermediate', 'PrpE': 'Synthesis intermediate'}    
    
    
pha_col_colors={}
phaset=[]
for name in validated_pha_count_df.columns:
    if name in pha_function_dict and not name in phaset:
        pha_col_colors[name]=''
        phaset.append(name)
functiondict={}
c=sns.xkcd_palette(['coral', 'blue', 'lavender', 'pumpkin', 'light grey',  'dandelion'])  
count=0
for p in phaset:
    if not pha_function_dict[p] in functiondict:
        functiondict[pha_function_dict[p]]=c[count]
        count+=1
for n in pha_col_colors:
    pha_col_colors[n]=functiondict[pha_function_dict[n]]  
pha_col_colors_ser=pd.Series(pha_col_colors)
pha_col_colors_lower={x.lower():pha_col_colors[x] for x in pha_col_colors}
col_legends=[]
for x in functiondict:
    col_legends.append(patches.Patch(edgecolor='black', facecolor=functiondict[x], label=x, lw=1))






sns.set(font_scale = 1.5)
sns.set(style="white")
hmcmap='GnBu'
norm = colors.TwoSlopeNorm(vmin=0, vcenter=2, vmax=9)
kwargs={'norm':norm}
clustcmap=truncate_colormap(plt.get_cmap('bone_r'), 0.05, 1)
##intensify the lower values in the heatmap to make them all visible
edgelw=0.5
xticklabelsize=10
##bbox [[xmin, ymin], [xmax, ymax]]. -keep all but xmin
cbar_bbox=Bbox([[0.5, 0.5], [0.7, 0.6]])
k=sns.clustermap(validated_pha_count_df.astype(float), metric="Euclidean",\
                 row_colors=f1cm_row_cols_axn_white, edgecolor='black', linewidth=edgelw,\
                     dendrogram_ratio=(0.16,0.06),yticklabels=True, cmap=clustcmap,\
                     cbar_kws = {'use_gridspec': False, 'orientation': 'horizontal',"ticks":[0, 1]}) #col_colors=pha_col_colors_ser
plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=8)
plt.setp(k.ax_heatmap.xaxis.get_majorticklabels(), fontsize=xticklabelsize, rotation=90)
legend2=k.ax_heatmap.legend(loc=(1.525,-0.175),handles=f1leg_handles, ncol=3, fontsize=10, facecolor='gray',columnspacing=1, labelspacing=0.25,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1, title='Taxonomy')#title='Taxonomy'
legend2col=k.ax_heatmap.legend(loc=(-0.25,-0.19),handles=col_legends, ncol=1, fontsize=10, facecolor='gray',columnspacing=1, labelspacing=0.25,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1.5, title='PHA function')
k.ax_heatmap.add_artist(legend2) 
k.ax_heatmap.add_artist(legend2col)     
for a in k.ax_row_dendrogram.collections:
    a.set_linewidth(2)
for a in k.ax_col_dendrogram.collections:
    a.set_linewidth(2)  
k.ax_heatmap.set_xticklabels(k.ax_heatmap.get_xmajorticklabels())
k.ax_heatmap.tick_params(right=False, bottom=False, pad=-5) 
k.ax_heatmap.set_ylabel('')
#########k.ax_heatmap.set_ylabel(labels=labelpad=0)
##extract new order of ylabels to transfer over to the heatmap
dgramyticks = [lbl.get_text() for lbl in k.ax_heatmap.get_yticklabels()]
##extract the clustermap heatmap position, to ensure fig 2 is the same size
##get position returns [[xmin, ymin], [xmax, ymax]]. -keep all but xmin
f1position=k.ax_heatmap.get_position().get_points()
# set the gridspec to only cover half of the figure
k.gs.update(left=0.05, right=0.45)  
gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.42)
# create axes within this new gridspec
ax2 = k.fig.add_subplot(gs2[0])  
#make an axis for heatmap cbar
cax = inset_axes(ax2,
                 width="20%",  # width: 40% of parent_bbox width
                 height="3%",  # height: 10% of parent_bbox height
                 loc='lower left',
                 bbox_to_anchor=(0.41, 1.030, 1.5, .82),#([x, y, width, height])
                 bbox_transform=ax2.transAxes,
                 borderpad=0)
##define the order of the heatmap by category of algae
ordercathm=['Brown', 'Brown/Green', 'Green', 'Red']
ordered_heatmap_columns=[]
for category in ordercathm:
    for n in list(cazy_genome_countdf_filt.columns):
        if not n in ordered_heatmap_columns:
            if n in list(CAZy_familiy_df.loc[CAZy_familiy_df['Algae'] == category]['CAZy family']):
                ordered_heatmap_columns.append(n)
k2=sns.heatmap(cazy_genome_countdf_filt.reindex(index=dgramyticks, columns=ordered_heatmap_columns), ax=ax2, cbar_ax=cax,\
               cbar_kws = {'use_gridspec': False, 'orientation': 'horizontal', "ticks":[0, 3, 6, 9]},\
                   yticklabels=True,cmap=hmcmap,edgecolor='white', linewidth=edgelw,linecolor='grey',\
                       **kwargs)#cbar_ax=ax_cbar
cax.tick_params(left=False, bottom=False,labelleft=False, right=True, length=0, width=0,labelright=True,pad=0)
f2position=k2.get_position().get_points()
##shift f2 left to make it look like a single heatmap
shift_f2_x_by=0.14
newposition=Bbox([[f2position[0][0]-shift_f2_x_by, f1position[0][1]], [f1position[1][0]-shift_f2_x_by, f1position[1][1]]])
k2.set_position(newposition)
#attain new position of f1 as new axis has been added
f1position2=k.ax_heatmap.get_position().get_points()
reduce_x_of_f1_by=0.18
f1position3=Bbox([[f1position2[0][0], f1position2[0][1]], [f1position2[1][0]-reduce_x_of_f1_by, f1position2[1][1]]])
k.ax_heatmap.set_position(f1position3)
##also change the size of the dendrogram by the same value
coldenrooffsetextra=0.01
dendropos=k.ax_col_dendrogram.get_position().get_points()
newdendropos=Bbox([[dendropos[0][0], dendropos[0][1]+coldenrooffsetextra], [dendropos[1][0]-reduce_x_of_f1_by, dendropos[1][1]]])
k.ax_col_dendrogram.set_position(newdendropos)
###offset row dendro by tiny amount
rowdenrooffsetextra=0.005
rowdendropos=k.ax_row_dendrogram.get_position().get_points()
newrowdendropos=Bbox([[rowdendropos[0][0], rowdendropos[0][1]], [rowdendropos[1][0]-rowdenrooffsetextra, rowdendropos[1][1]]])
k.ax_row_dendrogram.set_position(newrowdendropos)

####
k2.set_xticks([0.55+x for x in list(range(0, len(cazy_genome_countdf_filt.columns)))])
k2.set_xticklabels(ordered_heatmap_columns)
k2.tick_params(left=False, bottom=False,labelleft=False, right=True, length=5, width=2,labelright=True,pad=0)
k2.tick_params(axis='y', which='both', length=55, width=2,pad=0)
plt.setp(k2.yaxis.get_majorticklabels(), rotation=0, fontsize=10)
plt.setp(k2.xaxis.get_majorticklabels(), fontsize=xticklabelsize)
###shift cmaps and legends
kcaxcoords=k.cax.get_position().get_points()
shrink_kcax_by=0.05
#[[width, height], [xmax, ymax]]
k.cax.set_position(Bbox([[0.292, 0.95],[0.4, 0.971]]))
k.cax.tick_params(left=False, bottom=True,labelleft=False, right=True, length=0, width=0,labelright=True,pad=0)
## plot the color bar for f2
colheight=0.95
boxwidth = 1
boxheight = 0.8
##plot k2 col_colors
for caz in k2.get_xticklabels():
    a_x=caz.get_position()[0]-0.5
    a_y=caz.get_position()[1]-colheight
    ax2.add_patch(Rectangle(
        xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
        linewidth=1, edgecolor='black', lw=1, facecolor=f2cazy_colsx.T.loc[caz.get_text()][0], fill=True, clip_on=False))
##plot k1 col_colors
for function in k.ax_heatmap.get_xticklabels():
    a_x=function.get_position()[0]-0.5
    a_y=function.get_position()[1]-colheight
    k.ax_heatmap.add_patch(Rectangle(
        xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
        linewidth=1, edgecolor='black', lw=1, facecolor=pha_col_colors[function.get_text()], fill=True, clip_on=False))    
colheight=0.95
boxwidth = 0.6
boxheight = 0.8
##plot k1 row_colors
for function in k.ax_heatmap.get_yticklabels():
    a_x=function.get_position()[0]-1.75
    a_y=function.get_position()[1]-0.4
    k.ax_heatmap.add_patch(Rectangle(
        xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
        linewidth=1, edgecolor='black', lw=1, facecolor=f1cm_row_cols_axn[function.get_text()], fill=True, clip_on=False))    
k2.tick_params(axis='x', which='both', pad=-5)
ax2.text(1.25,-4.3, "PHA genes (n)", fontsize=10)
ax2.text(10.5,-4.3, "CAZy genes (n)", fontsize=10)
#remove the ylabels from the first figure
k.ax_heatmap.set_yticklabels(['' for x in k.ax_heatmap.get_yticklabels()])

###
algaeleg=ax2.legend(handles=f2cazy_leg_handlesx, labels=f2cazy_leg_labels,loc=(0.78, 1.02), ncol=2,  fontsize=10, facecolor='gray',columnspacing=0.2, labelspacing=0.0,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1.1, handleheight=0.4, handletextpad=0.2)
algaeleg.set_title('Target algae group',prop={'size':10})
for _, spine in k2.spines.items():
    spine.set_visible(True)
for _, spine in k.ax_heatmap.spines.items():
    spine.set_visible(True)
###add third axis for density plot
gs3=matplotlib.gridspec.GridSpec(1,1, left=0.7, right=.8)
# create axes within this new gridspec
ax3 = k.fig.add_subplot(gs3[0])
###set axis height same as heatmap
hmpos=ax2.get_position().get_points()
ax3pos=ax3.get_position().get_points()
densitypos=Bbox([[ax3pos[0][0], ax3pos[0][1]], [hmpos[1][0], hmpos[1][1]]])
ax3.set_position(densitypos)
###format data
###red, green, blue density
uniquealgaedict={}
totalalgaedict={}
##loop over genome axn
for indx in cazy_genome_countdf_filt.reindex(dgramyticks).index:
    uniquealgaedict[indx]={'Green':0, 'Red':0, 'Brown':0}
    totalalgaedict[indx]={'Green':0, 'Red':0, 'Brown':0}
    ##loop over the columns which are cazy names
    for col in cazy_genome_countdf_filt.reindex(dgramyticks).columns:
        ###check value not== 0:
        if not cazy_genome_countdf_filt.reindex(dgramyticks).loc[indx][col]==0:
            targetalgae=CAZy_familiy_df.loc[CAZy_familiy_df['CAZy family']==col]['Algae'].iloc[0]
            if targetalgae in uniquealgaedict[indx]:
                 uniquealgaedict[indx][targetalgae]=uniquealgaedict[indx][targetalgae]+1
                 totalalgaedict[indx][targetalgae]=totalalgaedict[indx][targetalgae]+cazy_genome_countdf_filt.reindex(dgramyticks).loc[indx][col]
t=pd.DataFrame.from_dict(uniquealgaedict).T.reindex(dgramyticks)
t=pd.DataFrame.from_dict(totalalgaedict).T.reindex(dgramyticks)
#size of bar=0.75
ypos=list(range(0, len(t)))[::-1]
msmult=20
##plot scatter size [x1,x2],[y1,y2],markersize=[], marker='o', markerfacecolor=, markeredgecolor=
ax3.scatter([0.5]*len(t), ypos, s=list(t['Green'].astype(float)*msmult), marker='o', color='#11ff11')
ax3.scatter([1.1]*len(t), ypos, s=list(t['Brown'].astype(float)*msmult), marker='o', color='#d2a968')
ax3.scatter([1.7]*len(t), ypos, s=list(t['Red'].astype(float)*msmult), marker='o', color='#f03f2e')
ax3.set_ylim([0-0.5, len(t)-1+0.5])
ax3.set_xlim([0,2])
ax3.spines['right'].set_visible(False)
ax3.set_xticks([])
ax3.set_yticks([])
plt.tight_layout()
#fig=plt.gcf()
#fig.savefig("22g_PHA_CAZy_heatmaps_paper_n70_total.png", dpi=750, bbox_inches='tight')






#########################################################
############ FORMAT AND LOAD FOR THE 40K genomes ########
#########################################################
CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)

pha_gene_list=['PhaA', 'PhaC', 'PhaE', 'PhaG', 'PhaD', 'PhaF', 'FabG', 'PhaR', 'PhaZ', 'PhaC2', 'maoC', 'PhaB',\
     'PhaP', 'PhaJ', 'PhaQ', 'PhaI', 'Pct', 'Pcs', 'PrpE', 'PduP', 'Cat2', 'HadA']

master_PHA_genome_df_C_C2only_kwmatches=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_new_C_C2_E_RefSeqs_biosample_KW_TAX.csv", sep=',', index_col=0, dtype='unicode')

with open('22g_40k_UID_gb_to_axn.json', 'r') as fp:
    bac_UID_gb_to_axn = json.load(fp)

master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion']=master_PHA_genome_df_C_C2only_kwmatches.index.map({int(k): v for (k, v) in bac_UID_gb_to_axn.items()})
cazy_genome_countdf=pd.read_json('Genome_CAZy_count_dict_40k.json')  
cazy_genome_countdf=cazy_genome_countdf.T  
###get rid of columns that have no data
cazy_genome_countdf_filt=cazy_genome_countdf[cazy_genome_countdf.columns[cazy_genome_countdf.sum(axis=0)>0]]
validated_pha_count_df=master_PHA_genome_df_C_C2only_kwmatches[pha_gene_list]
validated_pha_count_df['Genome accession']=validated_pha_count_df.index.map(master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion'])
merged = pd.merge(validated_pha_count_df, cazy_genome_countdf, left_on='Genome accession', right_index=True)
#merged=merged.drop(['Genome accession'], axis=1)
merged=merged.set_index('Genome accession')
validated_pha_count_df=validated_pha_count_df.set_index('Genome accession')
###filter out columns with no data
validated_pha_count_df=validated_pha_count_df[validated_pha_count_df.columns[validated_pha_count_df.astype(float).sum(axis=0)>0]]

merged=merged.astype(float).astype(int)
###get rid of columns that have no data
merged=merged[merged.columns[merged.sum(axis=0)>0]]
###get rid of rows that have no data
merged=merged[merged.columns[merged.sum(axis=0)>1]]

#merged.to_csv("22g_2_ALL_taxa_Master_PHA_genome_MERGED_40k_axn.csv", sep=',')
#validated_pha_count_df.to_csv("Validated_pha_count_df.csv", sep=','))

#merged=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_MERGED_40k.csv", sep=',', index_col=0, dtype='unicode')
#validated_pha_count_df=pd.read_csv("Validated_pha_count_df.csv", sep=',', index_col=0, dtype='unicode')

#change all gene densities to one, so that unique counts can be determined and use to filter
cazy_genome_countdf_all_eql_1=cazy_genome_countdf_filt.copy(deep=True)
cazy_genome_countdf_all_eql_1[cazy_genome_countdf_all_eql_1>1]=1
####use this to determine the count per green, brown, or red algae types 
greenlist=[x for x in list(CAZy_familiy_df[CAZy_familiy_df['Algae'].str.contains('Green')]['CAZy family']) if x in cazy_genome_countdf_all_eql_1.columns and not x=='GH3']
redlist=[x for x in list(CAZy_familiy_df[CAZy_familiy_df['Algae'].str.contains('Red')]['CAZy family']) if x in cazy_genome_countdf_all_eql_1.columns and not x=='GH3']
brownlist=[x for x in list(CAZy_familiy_df[CAZy_familiy_df['Algae'].str.contains('Brown')]['CAZy family']) if x in cazy_genome_countdf_all_eql_1.columns and not x=='GH3']
g_unique_filter=3
r_unique_filter=3
b_unique_filter=3
greenidx=cazy_genome_countdf_all_eql_1[cazy_genome_countdf_all_eql_1[greenlist].sum(axis=1)>g_unique_filter]
redidx=cazy_genome_countdf_all_eql_1[cazy_genome_countdf_all_eql_1[redlist].sum(axis=1)>r_unique_filter]
brownidx=cazy_genome_countdf_all_eql_1[cazy_genome_countdf_all_eql_1[brownlist].sum(axis=1)>b_unique_filter]
cazy_genome_uniques_filt=merged.loc[pd.DataFrame(greenidx+redidx+brownidx).index]
###get rid of columns that have no data
cazy_genome_uniques_filt=cazy_genome_uniques_filt[cazy_genome_uniques_filt.columns[cazy_genome_uniques_filt.sum(axis=0)>0]]
###get rid of rows that have no data
cazy_genome_uniques_filt=cazy_genome_uniques_filt[cazy_genome_uniques_filt.columns[cazy_genome_uniques_filt.sum(axis=0)>1]]

filtereddf=merged[cazy_genome_countdf_all_eql_1.sum(axis=1)>3].astype(float)
validated_pha_count_df_filt=validated_pha_count_df[validated_pha_count_df.index.isin(cazy_genome_uniques_filt.index)].astype(float)
###get rid of columns that have no data
validated_pha_count_df_filt=validated_pha_count_df_filt[validated_pha_count_df_filt.columns[validated_pha_count_df_filt.sum(axis=0)>0]]
###get rid of rows that have no data
validated_pha_count_df_filt=validated_pha_count_df_filt[validated_pha_count_df_filt.columns[validated_pha_count_df_filt.sum(axis=0)>1]]

####top xx percentile for phylogenetic tree
percentile=0.15
grnpct=greenidx[greenidx.sum(axis=1).rank(ascending=False, pct=True)<=percentile]
redpct=redidx#[redidx.sum(axis=1).rank(ascending=False, pct=True)<=percentile]
brwnpct=brownidx[brownidx.sum(axis=1).rank(ascending=False, pct=True)<=percentile]
toppct_set=set(list(grnpct.index)+list(redpct.index)+list(brwnpct.index))

#sns.clustermap(filtereddf)


#########
import matplotlib.gridspec  
from matplotlib.transforms import Bbox  
import matplotlib.colors as colors
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import colors

def truncate_colormap(cmap, minval=0.1, maxval=0.9, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def make_algae_cazy_clustermap_color_map(input_df, lw):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':(0.4, 0.6, 0)}
    legend_handles=[]
    legrepeats=[]
    legend_labels=[]
    for idx in list(range(0, len(input_df))):
        rowcoldict[input_df.iloc[idx]['CAZy family']]=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()]
        if not CAZy_familiy_df.iloc[idx]['Algae'] in legrepeats:
            legend_handles.append(patches.Patch(facecolor=colrs[str(CAZy_familiy_df.iloc[idx]['Algae']).lower()],\
                                                label=str(CAZy_familiy_df.iloc[idx]['Algae']), edgecolor='black', lw=lw))
            legrepeats.append(str(CAZy_familiy_df.iloc[idx]['Algae']))
            legend_labels.append(str(CAZy_familiy_df.iloc[idx]['Algae']))
    return pd.DataFrame.from_dict([rowcoldict]), legend_handles, legend_labels

def make_clustermap_color_map(taxlevel, dataframe, name, lw):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    dataframe=dataframe.fillna('Unknown')
    num_unique_tax=len(set(list(dataframe[taxlevel])))
    #sns.palplot(sns.color_palette("husl", 12))
    colrs=list(sns.color_palette("Paired",14))+list(sns.color_palette("hls", 12)[::-1])
    c=0
    legend_handles_group=[]
    for x in set(list(dataframe[taxlevel])):
        if not x=='Unknown':
             rowcoldict[x]=(colrs[c])
             legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=rowcoldict[x], label=x, lw=lw))
             c+=1
        elif x=='Unknown':
             rowcoldict[x]=(1,1,1,1)
             legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=(1,1,1,1), label='Unknown', lw=lw))
             c+=1            
    ###use this color dict to create a list of colours from the dataframe
    col_list=[]
    col_list_white=[]
    for taxa in list(dataframe[taxlevel]):
        if not taxa=='Unknown':
             col_list.append(rowcoldict[taxa])
             col_list_white.append((1,1,1,1))
        elif taxa=='Unknown':
             col_list.append((1,1,1,1))
             col_list_white.append((1,1,1,1))                
    df_colors=pd.Series(col_list, dataframe.Genome_accesion, name=str(name))
    df_colors_white=pd.Series(col_list_white, dataframe.Genome_accesion, name=str(name))
    return df_colors_white,df_colors, legend_handles_group

linewidth=0
taxlevel='family'
min_tax_groupby_size=1
##filter to only this included above
fullfilt=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches.Genome_accesion.isin(validated_pha_count_df_filt.index)]
##filter by groupyby minimum taxa size
fullfilt2=fullfilt[fullfilt[taxlevel].groupby(fullfilt[taxlevel]).transform('size')>min_tax_groupby_size]
###filter main dataframe by new values by cross referencing genome accessions as above
validated_pha_count_df_filt2=validated_pha_count_df_filt[validated_pha_count_df_filt.index.isin(fullfilt2.Genome_accesion)]

f1cm_row_cols_axn_white,f1cm_row_cols_axn, f1leg_handles=make_clustermap_color_map(taxlevel, fullfilt2, '', linewidth)

f2cazy_colsx, f2cazy_leg_handlesx, f2cazy_leg_labels=make_algae_cazy_clustermap_color_map(CAZy_familiy_df, linewidth)

#make pha dict grouping genes by function
pha_function_dict={'PhaA':'Synthesis intermediate', 'PhaC':'Synthesis terminal', 'PhaE':'Synthesis terminal',\
                   'PhaG':'Synthesis intermediate', 'PhaD':'Phasin regulator', 'PhaF':'Phasin', 'PhaP':'Phasin',\
                   'PhaZ':'Depolymerase', 'PhaB':'Synthesis intermediate', 'PhaR':'Synthesis regulator', 'maoC':'Synthesis intermediate',\
                   'PhaC2':'Synthesis terminal','PhaI':'Phasin', 'FabG':'Synthesis intermediate',\
                   'PhaJ':'Synthesis intermediate', 'PhaQ':'Phasin regulator', 'PhaP':'Phasin',\
                   'Pcs':'Synthesis intermediate', 'PrpE': 'Synthesis intermediate'}
pha_col_colors={}
phaset=[]
for name in cazy_genome_uniques_filt.columns:
    if name in pha_function_dict and not name in phaset:
        pha_col_colors[name]=''
        phaset.append(name)
functiondict={}
c=sns.xkcd_palette(['coral', 'blue', 'lavender', 'pumpkin', 'light grey', 'dandelion'])  
count=0
for p in phaset:
    if not pha_function_dict[p] in functiondict:
        functiondict[pha_function_dict[p]]=c[count]
        count+=1
for n in pha_col_colors:
    pha_col_colors[n]=functiondict[pha_function_dict[n]]  
pha_col_colors_ser=pd.Series(pha_col_colors)
col_legends=[]
for x in functiondict:
    col_legends.append(patches.Patch(edgecolor='black', facecolor=functiondict[x], label=x, lw=1))

pha_col_lower={x.lower():pha_col_colors[x] for x in pha_col_colors}
pha_func_dict_lower={x.lower():pha_function_dict[x] for x in pha_function_dict}

########################################################################
#     CREATE TAX LEVEL GROUPBY DICT FOR CO-OCCURANCE NETWORK            #
########################################################################
def make_clustermap_color_map_from_index(index_level, desired_cbar_level, dataframe, master_dataframe, name, lw):
    ##create a dictionary to map taxa to colour e.g. taxa:red
    rowcoldict={}
    master_dataframe=master_dataframe.fillna('Unknown')
    listoftaxa=set(list(dataframe.index))
    ###get list of taxa at the higher desired tax level
    desired_tax_list=[]
    for tax in listoftaxa:
        if not tax=='Unknown':
             returned_tax=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches[index_level].\
                                                             fillna('Unknown').str.contains(tax)][desired_cbar_level].iloc[0]
        elif tax =='Unknown'
             returned_tax='Unknown'
        if returned_tax not in desired_tax_list:
            desired_tax_list.append(returned_tax)
    num_unique_tax=len(desired_tax_list)
    ##DEFINE PALETTE BASED ON THE NUMBER OF HIGHER LEVEL TAXA RETURNED
    colrs=list(sns.color_palette("Paired",14))+list(sns.color_palette("hls", 12)[::-1])
    c=0
    legend_handles_group=[]
    for x in listoftaxa:
        if not x=='Unknown':
             rowcoldict[x]=(colrs[c])
             legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=rowcoldict[x], label=x, lw=lw))
             c+=1
        elif x=='Unknown':
             rowcoldict[x]=(1,1,1,1)
             legend_handles_group.append(patches.Patch(edgecolor='black', facecolor=(1,1,1,1), label='Unknown', lw=lw))
             c+=1            
    ###use this color dict to create a list of colours from the dataframe
    col_list=[]
    col_list_white=[]
    for taxa in list(dataframe[taxlevel]):
        if not taxa=='Unknown':
             col_list.append(rowcoldict[taxa])
             col_list_white.append((1,1,1,1))
        elif taxa=='Unknown':
             col_list.append((1,1,1,1))
             col_list_white.append((1,1,1,1))                
    df_colors=pd.Series(col_list, dataframe.Genome_accesion, name=str(name))
    df_colors_white=pd.Series(col_list_white, dataframe.Genome_accesion, name=str(name))
    return df_colors_white,df_colors, legend_handles_group 

grp_by_tax_level='order'
gpby_input_df=master_PHA_genome_df_C_C2only_kwmatches
carryover_columns=[x for x in pha_gene_list if x in list(gpby_input_df.columns)]
tax_coooc_df=master_PHA_genome_df_C_C2only_kwmatches[carryover_columns].astype(float)
###pha only
tax_coooc_df[grp_by_tax_level]=master_PHA_genome_df_C_C2only_kwmatches[grp_by_tax_level]
tax_coooc_df_grouped_PHA_only=tax_coooc_df.groupby(by=grp_by_tax_level).sum()
tax_coooc_df_grouped_PHA_only_log2=np.log2(tax_coooc_df_grouped_PHA_only).replace([np.inf, -np.inf], 0)


sns.clustermap(tax_coooc_df_grouped_PHA_only_log2[tax_coooc_df_grouped_PHA_only_log2.sum(axis=1)>0])
###pha and cazy
#merged contains all pha genes for all 41,745 genomes with at least 1 identified cazyme
#map the taxonomic info at the selected level back into merged so it can be grouped
merged_tax=merged.merge(master_PHA_genome_df_C_C2only_kwmatches.set_index("Genome_accesion")[grp_by_tax_level], left_index=True,\
            right_index=True)
merged_tax_grouped=merged_tax.groupby(by=grp_by_tax_level).sum()    
merged_tax_grouped_log=np.log2(merged_tax_grouped).replace([np.inf, -np.inf], 0)

f1cm_row_cols_axn_white,f1cm_row_cols_axn, f1leg_handles=make_clustermap_color_map('class', tax_coooc_df_grouped_PHA_only_log2, '', linewidth)

########################################################################
#     CREATE TAX LEVELGROUPBY DICT FOR CO-OCCURANCE NETWORK            #
#                         --END--                                      #
########################################################################


d40k=True

sns.set(font_scale = 1.5)
sns.set(style="white")
hmcmap='GnBu'
norm = colors.TwoSlopeNorm(vmin=0, vcenter=2, vmax=9)
kwargs={'norm':norm}
clustcmap=truncate_colormap(plt.get_cmap('bone_r'), 0.05, 1)
##intensify the lower values in the heatmap to make them all visible
xticklabelsize=10
##bbox [[xmin, ymin], [xmax, ymax]]. -keep all but xmin
cbar_bbox=Bbox([[0.5, 0.5], [0.7, 0.6]])
k=sns.clustermap(validated_pha_count_df_filt2.astype(float), metric="Euclidean",\
                 row_colors=f1cm_row_cols_axn_white, edgecolor='black', linewidth=linewidth,\
                     dendrogram_ratio=(0.16,0.06),yticklabels=True, cmap=clustcmap,\
                     cbar_kws = {'use_gridspec': False, 'orientation': 'horizontal',"ticks":[0, 1]}) #col_colors=pha_col_colors_ser
#change rowcolor tick labels
k.ax_row_colors.tick_params(pad=-5)
plt.setp(k.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=8)
plt.setp(k.ax_heatmap.xaxis.get_majorticklabels(), fontsize=xticklabelsize, rotation=90)
legend2=k.ax_heatmap.legend(loc=(1.525,-0.25),handles=f1leg_handles, ncol=2, fontsize=10, facecolor='gray',columnspacing=1, labelspacing=0.25,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1, title='Taxonomy')#title='Taxonomy'
legend2col=k.ax_heatmap.legend(loc=(-0.25,-0.19),handles=col_legends, ncol=1, fontsize=10, facecolor='gray',columnspacing=1, labelspacing=0.25,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1.5, title='PHA function')
k.ax_heatmap.add_artist(legend2) 
k.ax_heatmap.add_artist(legend2col)     
for a in k.ax_row_dendrogram.collections:
    a.set_linewidth(2)
for a in k.ax_col_dendrogram.collections:
    a.set_linewidth(2)  
k.ax_heatmap.set_xticklabels(k.ax_heatmap.get_xmajorticklabels())
k.ax_heatmap.tick_params(right=False, bottom=False, pad=-5) 
k.ax_heatmap.set_ylabel('')
#########k.ax_heatmap.set_ylabel(labels=labelpad=0)
##extract new order of ylabels to transfer over to the heatmap
dgramyticks = [lbl.get_text() for lbl in k.ax_heatmap.get_yticklabels()]
##extract the clustermap heatmap position, to ensure fig 2 is the same size
##get position returns [[xmin, ymin], [xmax, ymax]]. -keep all but xmin
f1position=k.ax_heatmap.get_position().get_points()
# set the gridspec to only cover half of the figure
k.gs.update(left=0.05, right=0.45)  
gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.42)
# create axes within this new gridspec
ax2 = k.fig.add_subplot(gs2[0])  
#make an axis for heatmap cbar
cax = inset_axes(ax2,
                 width="20%",  # width: 40% of parent_bbox width
                 height="3%",  # height: 10% of parent_bbox height
                 loc='lower left',
                 bbox_to_anchor=(0.41, 1.030, 1.5, .82),#([x, y, width, height])
                 bbox_transform=ax2.transAxes,
                 borderpad=0)
##define the order of the heatmap by category of algae
ordercathm=['Brown', 'Brown/Green', 'Green', 'Red']
ordered_heatmap_columns=[]
for category in ordercathm:
    for n in list(cazy_genome_countdf_filt.columns):
        if not n in ordered_heatmap_columns:
            if n in list(CAZy_familiy_df.loc[CAZy_familiy_df['Algae'] == category]['CAZy family']):
                ordered_heatmap_columns.append(n)
k2=sns.heatmap(cazy_genome_countdf_filt.reindex(index=dgramyticks, columns=ordered_heatmap_columns), ax=ax2, cbar_ax=cax,\
               cbar_kws = {'use_gridspec': False, 'orientation': 'horizontal', "ticks":[0, 3, 6, 9]},\
                   yticklabels=False,cmap=hmcmap,edgecolor='white', linewidth=linewidth,linecolor='grey',\
                       **kwargs)#cbar_ax=ax_cbar
cax.tick_params(left=False, bottom=False,labelleft=False, right=True, length=0, width=0,labelright=True,pad=0)
f2position=k2.get_position().get_points()
##shift f2 left to make it look like a single heatmap
shift_f2_x_by=0.14
newposition=Bbox([[f2position[0][0]-shift_f2_x_by, f1position[0][1]], [f1position[1][0]-shift_f2_x_by, f1position[1][1]]])
k2.set_position(newposition)
#attain new position of f1 as new axis has been added
f1position2=k.ax_heatmap.get_position().get_points()
reduce_x_of_f1_by=0.18
f1position3=Bbox([[f1position2[0][0], f1position2[0][1]], [f1position2[1][0]-reduce_x_of_f1_by, f1position2[1][1]]])
k.ax_heatmap.set_position(f1position3)
##also change the size of the dendrogram by the same value
coldenrooffsetextra=0.01
dendropos=k.ax_col_dendrogram.get_position().get_points()
newdendropos=Bbox([[dendropos[0][0], dendropos[0][1]+coldenrooffsetextra], [dendropos[1][0]-reduce_x_of_f1_by, dendropos[1][1]]])
k.ax_col_dendrogram.set_position(newdendropos)
###offset row dendro by tiny amount
rowdenrooffsetextra=0.005
#rowdenrooffsetextra=0.000
rowdendropos=k.ax_row_dendrogram.get_position().get_points()
newrowdendropos=Bbox([[rowdendropos[0][0], rowdendropos[0][1]], [rowdendropos[1][0]-rowdenrooffsetextra, rowdendropos[1][1]]])
k.ax_row_dendrogram.set_position(newrowdendropos)

####
k2.set_xticks([0.55+x for x in list(range(0, len(cazy_genome_countdf_filt.columns)))])
k2.set_xticklabels(ordered_heatmap_columns)
k2.tick_params(left=False, bottom=False,labelleft=False, right=True, length=5, width=2,labelright=True,pad=0)
k2.tick_params(axis='y', which='both', length=55, width=2,pad=0)
plt.setp(k2.yaxis.get_majorticklabels(), rotation=0, fontsize=10)
plt.setp(k2.xaxis.get_majorticklabels(), fontsize=xticklabelsize)
###shift cmaps and legends
kcaxcoords=k.cax.get_position().get_points()
shrink_kcax_by=0.05
#[[width, height], [xmax, ymax]]
k.cax.set_position(Bbox([[0.292, 0.947],[0.4, 0.969]]))
k.cax.tick_params(left=False, bottom=True,labelleft=False, right=True, length=0, width=0,labelright=True,pad=0)
## plot the color bar for f2
colheight=0.95
boxwidth = 1
boxheight = 0.8
d40kextra=2.2
if d40k==True:
    boxheight=boxheight+d40kextra
##plot k2 col_colors
for caz in k2.get_xticklabels():
    a_x=caz.get_position()[0]-0.5
    a_y=caz.get_position()[1]-colheight
    if d40k==True:
         ax2.add_patch(Rectangle(
             xy=(a_x, a_y-d40kextra) ,width=boxwidth, height=boxheight,
             linewidth=linewidth, edgecolor='black', lw=linewidth, facecolor=f2cazy_colsx.T.loc[caz.get_text()][0], fill=True, clip_on=False))
    else:
         ax2.add_patch(Rectangle(
             xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
             linewidth=linewidth, edgecolor='black', lw=linewidth, facecolor=f2cazy_colsx.T.loc[caz.get_text()][0], fill=True, clip_on=False))
##plot k1 col_colors
for function in k.ax_heatmap.get_xticklabels():
    a_x=function.get_position()[0]-0.5
    a_y=function.get_position()[1]-colheight
    if d40k==True:
         k.ax_heatmap.add_patch(Rectangle(
             xy=(a_x, a_y-d40kextra), width=boxwidth, height=boxheight,
             linewidth=linewidth, edgecolor='black', lw=linewidth, facecolor=pha_col_colors[function.get_text()], fill=True, clip_on=False))    
    else:
         k.ax_heatmap.add_patch(Rectangle(
             xy=(a_x, a_y), width=boxwidth, height=boxheight,
             linewidth=linewidth, edgecolor='black', lw=linewidth, facecolor=pha_col_colors[function.get_text()], fill=True, clip_on=False))    

colheight=0.95
boxwidth = 0.6
boxheight = 1
d40kextrawidth=0.75
moveleft=0.1
if d40k==True:
    boxwidth=boxwidth+d40kextrawidth
##plot k1 row_colors
for function in k.ax_heatmap.get_yticklabels():
    a_x=function.get_position()[0]-1.75
    a_y=function.get_position()[1]-0.4
    k.ax_heatmap.add_patch(Rectangle(
        xy=(a_x-(boxwidth/2)+moveleft, a_y) ,width=boxwidth, height=boxheight,
        linewidth=linewidth, edgecolor='black', lw=linewidth, facecolor=f1cm_row_cols_axn[function.get_text()], fill=True, clip_on=False))    
k2.tick_params(axis='x', which='both', pad=-5)
ax2.text(1.25,-20, "PHA genes (n)", fontsize=10)
ax2.text(12.5,-20, "CAZy genes (n)", fontsize=10)
#remove the ylabels from the first figure
k.ax_heatmap.set_yticklabels(['' for x in k.ax_heatmap.get_yticklabels()])

###
algaeleg=ax2.legend(handles=f2cazy_leg_handlesx, labels=f2cazy_leg_labels,loc=(0.78, 1.02), ncol=2,  fontsize=10, facecolor='gray',columnspacing=0.2, labelspacing=0.0,framealpha=0.1, frameon=True, fancybox=True, edgecolor='black', handlelength=1.1, handleheight=0.4, handletextpad=0.2)
algaeleg.set_title('Target algae group',prop={'size':10})
for _, spine in k2.spines.items():
    spine.set_visible(True)
for _, spine in k.ax_heatmap.spines.items():
    spine.set_visible(True)
###add third axis for density plot
gs3=matplotlib.gridspec.GridSpec(1,1, left=0.7, right=.8)
# create axes within this new gridspec
ax3 = k.fig.add_subplot(gs3[0])
###set axis height same as heatmap
hmpos=ax2.get_position().get_points()
ax3pos=ax3.get_position().get_points()
densitypos=Bbox([[ax3pos[0][0], ax3pos[0][1]], [hmpos[1][0], hmpos[1][1]]])
ax3.set_position(densitypos)
###format data
###red, green, blue density
uniquealgaedict={}
totalalgaedict={}
##loop over genome axn
for indx in cazy_genome_countdf_filt.reindex(dgramyticks).index:
    uniquealgaedict[indx]={'Green':0, 'Red':0, 'Brown':0}
    totalalgaedict[indx]={'Green':0, 'Red':0, 'Brown':0}
    ##loop over the columns which are cazy names
    for col in cazy_genome_countdf_filt.reindex(dgramyticks).columns:
        ###check value not== 0:
        if not cazy_genome_countdf_filt.reindex(dgramyticks).loc[indx][col]==0:
            targetalgae=CAZy_familiy_df.loc[CAZy_familiy_df['CAZy family']==col]['Algae'].iloc[0]
            if targetalgae in uniquealgaedict[indx]:
                 uniquealgaedict[indx][targetalgae]=uniquealgaedict[indx][targetalgae]+1
                 totalalgaedict[indx][targetalgae]=totalalgaedict[indx][targetalgae]+cazy_genome_countdf_filt.reindex(dgramyticks).loc[indx][col]
t=pd.DataFrame.from_dict(uniquealgaedict).T.reindex(dgramyticks)
#t=pd.DataFrame.from_dict(totalalgaedict).T.reindex(dgramyticks)
#size of bar=0.75
ypos=list(range(0, len(t)))[::-1]
msmult=10
##plot scatter size [x1,x2],[y1,y2],markersize=[], marker='o', markerfacecolor=, markeredgecolor=
ax3.scatter([0.5]*len(t), ypos, s=list(t['Green'].astype(float)*msmult), marker='o', alpha=0.5, color='#11ff11')
ax3.scatter([1.1]*len(t), ypos, s=list(t['Brown'].astype(float)*msmult), marker='o', alpha=0.5, color='#d2a968')
ax3.scatter([1.7]*len(t), ypos, s=list(t['Red'].astype(float)*msmult), marker='o', alpha=0.5, color='#f03f2e')
ax3.set_ylim([0-0.5, len(t)-1+0.5])
ax3.set_xlim([0,2])
ax3.spines['right'].set_visible(False)
ax3.set_xticks([0.5,1.1,1.7])
ax3.set_xticklabels(['G', 'B', 'R'])
ax3.set_yticks([])
plt.tight_layout()
#fig=plt.gcf()
#fig.savefig("40k_22g_PHA_CAZy_heatmaps_paper_n70_total.png", dpi=750, bbox_inches='tight')











################################################################
###  for ALL NON-HIGH QUALITY GENOMES### all 801 genomes    ####
################################################################
cazy_genome_countdf_all=pd.read_json('Genome_CAZy_count_dict_small.json')  
cazy_genome_countdf_all=cazy_genome_countdf_all.T  

CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)
###use F6 as it is not filtered
master_PHA_genome_df_F6=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M.csv", sep=',', index_col=0, dtype='unicode')
GOI=['phaa', 'phac', 'phae', 'phag', 'phad', 'phaf', 'fabg', 'phar', 'phaz', 'phac2', 'maoc', 'phab',\
     'phap', 'phaj', 'phaq', 'phai', 'pct', 'pcs', 'prpe', 'pdup', 'cat2', 'hada']
validated_pha_count_df_all=master_PHA_genome_df_F6[pha_gene_list]
##map genome accession onto pha count dataframe  
validated_pha_count_df_all['Genome accession']=validated_pha_count_df_all.index.map(master_PHA_genome_df_F6['Genome_accesion'])
###merge the dataframes using the genomaccession from pha df to the index of the cazy df  
###merge is pd.merge(left df, right df, how left, how right)
merged_all = pd.merge(validated_pha_count_df_all, cazy_genome_countdf_all, left_on='Genome accession', right_index=True)
merged_all=merged_all.drop(['Genome accession'], axis=1)    
validated_pha_count_df=validated_pha_count_df_all.set_index('Genome accession')
merged_all=merged_all.astype(float).astype(int)
###get rid of columns that have no data
merged_all=merged_all[merged_all.columns[merged_all.sum(axis=0)>0]]
###get rid of rows that have no data
merged_all=merged_all[merged_all.columns[merged_all.sum(axis=0)>1]]


################################################################
###                  for ALL 40K genomes                    ####
################################################################
cazy_genome_countdf_all=pd.read_json('Genome_CAZy_count_dict_40k.json')  
cazy_genome_countdf_all=cazy_genome_countdf_all.T  

CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)
###use F6 as it is not filtered
master_PHA_genome_df_F6=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_MERGED_40k_axn.csv", sep=',', index_col=0, dtype='unicode')
GOI=['phaa', 'phac', 'phae', 'phag', 'phad', 'phaf', 'fabg', 'phar', 'phaz', 'phac2', 'maoc', 'phab',\
     'phap', 'phaj', 'phaq', 'phai', 'pct', 'pcs', 'prpe', 'pdup', 'cat2', 'hada']

validated_pha_count_df_all=master_PHA_genome_df_F6[[gene for gene in pha_gene_list if gene in master_PHA_genome_df_F6.columns]]
##map genome accession onto pha count dataframe  
validated_pha_count_df_all['Genome accession']=validated_pha_count_df_all.index.map(master_PHA_genome_df_F6['Genome accession'])
###merge the dataframes using the genomaccession from pha df to the index of the cazy df  
###merge is pd.merge(left df, right df, how left, how right)
merged_all = pd.merge(validated_pha_count_df_all, cazy_genome_countdf_all, left_on='Genome accession', right_index=True)
merged_all=merged_all.drop(['Genome accession'], axis=1)    
validated_pha_count_df=validated_pha_count_df_all.set_index('Genome accession')
merged_all=merged_all.astype(float).astype(int)
###get rid of columns that have no data
merged_all=merged_all[merged_all.columns[merged_all.sum(axis=0)>0]]
###get rid of rows that have no data
merged_all=merged_all[merged_all.columns[merged_all.sum(axis=0)>1]]

################################################################
###                    NETWORKX                             ####
################################################################
import numpy as np
import networkx as nx

#### Use this for ALL analysis ###

mergedint=merged_all.astype(float).astype(int)
#mergedint=merged.astype(float).astype(int)
###remove data less than x, i.e. filter edges by minimum value
minimum_edges=2
##filter nodes that are not present in any genome
mergedint=mergedint[mergedint.columns[mergedint.sum(axis=0)>0]]
##construct co-occurance matrix
coocc_count = mergedint.T.dot(mergedint)

coocc_corr = mergedint.corr(method='pearson')

# create a graph from your dataframe
fig, ax = plt.subplots(figsize=(8,8))
plt.axis('off')
G = nx.from_pandas_adjacency(coocc_count)
#G = nx.from_pandas_adjacency(coocc_corr)
layout = nx.spring_layout(G,iterations=10,k=1)

nodecolor=[]
colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':(0.4, 0.6, 0)}
algaecoldict={}
for cazyfam in CAZy_familiy_df['CAZy family']:
    if cazyfam not in algaecoldict:
         algaecoldict[cazyfam]=colrs[str(CAZy_familiy_df.loc[CAZy_familiy_df['CAZy family'] == cazyfam]['Algae'].iloc[0]).lower()]
for label in G.nodes:
    if label in algaecoldict:
         nodecolor.append(algaecoldict[label])
    if not label in algaecoldict:
        nodecolor.append('skyblue')
###edge colors
edgelist=G.edges(data=True)
edgeweights=[]
edgeweights_mult=[]
edgemultiplier=5
for x in edgelist:
    weight=list(x[2].values())[0]
    edgeweights.append(weight)
    edgeweights_mult.append(weight*edgemultiplier)
##generate a normalised distribution for the colors based on max and min
norm = mpl.colors.Normalize(vmin=min(edgeweights), vmax=max(edgeweights)) 
cmapname='coolwarm'     
orig_cmap = mpl.cm.ScalarMappable(norm=norm, cmap=cmapname)
edgecolors=[orig_cmap.to_rgba(norm(x)) for x in edgeweights]
cbar_ax = fig.add_axes([0.09, 0.86, 0.25, 0.02])
cb=fig.colorbar(orig_cmap, cax=cbar_ax, orientation="horizontal")
cb.ax.set_xticklabels([-0.2, 0, 0.5, 1], weight='bold')
cb.ax.tick_params(labelsize=12) 
cb.outline.set_edgecolor('black')
cb.outline.set_linewidth(2)
nodeborders=['black']*len(G.nodes())
nodeborderlinewidths=[2]*len(G.nodes())
nx.draw(G, layout, node_color=nodecolor,edgecolors=nodeborders,linewidths=nodeborderlinewidths, edge_color=edgecolors, width=edgeweights_mult, with_labels=True, ax=ax)
#fig.savefig("example_network.png", dpi=500)

######################################
####    ARC STANKEY DIAGRAM       ####
#### NB/ requires edges from above####
####        using RAW counts      ####
######################################

from scipy import interpolate
fig, ax = plt.subplots(figsize=(14,8))
nodes=list(coocc_count.columns)
ordercat=['Green','Brown/Green', 'Brown', 'pha', 'Red']
nodealpha=1
ordercoldict={'Green':(0.2, 0.8, 0.2, nodealpha), 'Brown':(0.45, .25, 0, nodealpha), 'Brown/Green':(0.4, 0.6, 0, nodealpha), 'pha':(0.57, 0.57, 0.57, nodealpha), 'Red':(0.8, 0, 0, nodealpha)}
orderednodes_raw=[]
orderednodecolors_raw=[]
ordered_sum_raw=[]
colordict={}
###create a dict to store params
param_dict={}
params=['Nodes', 'Colors', 'Sizes']
for cat in ordercat:
    param_dict[cat]={}
    for para in params:
        param_dict[cat][para]=[]
for category in ordercat:
    if not category == 'pha':
        for n in nodes:
            if n in list(CAZy_familiy_df.loc[CAZy_familiy_df['Algae'] == category]['CAZy family']):
                param_dict[category]['Nodes'].append(n)
                param_dict[category]['Colors'].append(ordercoldict[category])
                param_dict[category]['Sizes'].append(coocc_count[n].sum())
                colordict[n]=ordercoldict[category]
                orderednodes_raw.append(n)
                orderednodecolors_raw.append(ordercoldict[category])
                ordered_sum_raw.append(coocc_count[n].sum())                
    if category == 'pha':
        for n in nodes:
            if n.lower().startswith(category) or n.lower().startswith('fab') or n.lower().startswith('pcs') or \
                n.lower().startswith('prp') or n.lower().startswith('mao') or n.lower().startswith('had'):
                param_dict[category]['Nodes'].append(n)
                param_dict[category]['Colors'].append(ordercoldict[category])
                param_dict[category]['Sizes'].append(coocc_count[n].sum())
                orderednodes_raw.append(n)
                orderednodecolors_raw.append(ordercoldict[category])  
                ordered_sum_raw.append(coocc_count[n].sum())  
                colordict[n]=ordercoldict[category]
###use the data in params_dict to order the data from large to small, then join all the data lists together
orderednodes=[]
orderednodecolors=[]
ordered_sum=[]
for cat in ordercat:
    ###sort by size, color, node and append each from sorted zip, to keep indexes matching
    ordered_sum=ordered_sum+[x for x, _, _, in sorted(zip(param_dict[cat]['Sizes'], param_dict[cat]['Colors'], param_dict[cat]['Nodes']), reverse=True)]                  
    orderednodecolors=orderednodecolors+[x for _, x, _, in sorted(zip(param_dict[cat]['Sizes'], param_dict[cat]['Colors'], param_dict[cat]['Nodes']), reverse=True)]    
    orderednodes=orderednodes+[x for _, _, x, in sorted(zip(param_dict[cat]['Sizes'], param_dict[cat]['Colors'], param_dict[cat]['Nodes']), reverse=True)]    

count=list(range(0, len(nodes)))
##y=y position of x axis nodes
y=[0,3,0]
plt.figure(figsize=(10,10))
separation=5
x=[1, 1+separation*1, 1+separation*2]
#x = [1,4,7]
start=x[-1]
node_position_dict={}
xtickpositions=[]
###normalise marker sizes
maxms=400
minms=20
markersizes_norm=[maxms if x >= maxms else x+minms if x <= minms else x for x in ordered_sum]
ms_multiplier=0.1
markersizes_norm=[x*ms_multiplier for x in markersizes_norm]
###plot the nodes in the loop below
for i in list(range(0,len(orderednodes))):
   #plot the x axis nodes with equal spacing e.g. 3 pts apart
   if i==0:
      x = x
   else:
     x[0]=start
     x[1]=x[0]+separation
     x[2]=x[1]+separation
   start=x[-1]
   new_x=[x[0],x[-1]]
   new_y=[y[0],y[-1]]
   if not orderednodes[i] in node_position_dict:
        node_position_dict[orderednodes[i]]={}
        node_position_dict[orderednodes[i]]['x_coords']=new_x
        node_position_dict[orderednodes[i]]['y_coords']=new_y
        node_position_dict[orderednodes[i]]['nodecolor']=orderednodecolors[i]
   ###plot nodes
   ax.plot(new_x[0],new_y[0],"o",color=orderednodecolors[i],markeredgewidth=2,ms=markersizes_norm[i], markeredgecolor='grey')
   xtickpositions.append(new_x[0])
   print(orderednodes[i], new_x, new_y)
   ##check the node is not last and plot lines
   line='no'
   if line == 'yes':
        if not i == len(orderednodes)-1:
            ax.plot(new_x,[0,0],color='grey',linewidth=5, zorder=-1)
   #xticklabs.append(nodes[i])
##plot between
plot_between=[('phaA','GH3'),('GH29','GH129')]
##get the max value of all arcs to set the y axis limit for the plot
ymaxvalue=0
##requires edges, weights
connection_count=0
##set max value for edgeweights
maxarc=50
arcedgeweights=[maxarc if x >= maxarc else x for x in edgeweights]
linewidth_multiplier=0.2
arcedgeweights=[x*linewidth_multiplier for x in arcedgeweights]
for node_pair in G.edges():
#for node_pair in plot_between:    
     arcstartxy=node_position_dict[node_pair[0]]
     arcendxy=node_position_dict[node_pair[1]]
     ##source
     x1, y1 = node_position_dict[node_pair[0]]['x_coords'][0],node_position_dict[node_pair[0]]['y_coords'][0]
     ##target
     x2, y2 = node_position_dict[node_pair[1]]['x_coords'][0],node_position_dict[node_pair[1]]['y_coords'][0]
     #swap values if the first x position is larger, or a whole circle will be drawn
     if x2<x1:
          x1, x2 = x2, x1
     # calculate the arc
     mxmy = mx, my = [(x1 + x2) / 2, (y1 + y2) / 2]
     r = np.sqrt((x1 - mx)**2 + (y1 - my)**2)
     width = 2 * r
     height = 2 * r
     if height/2 >= ymaxvalue:
         ymaxvalue=height/2
     start_angle = np.arctan2(y1 - my, x1 - mx) * 180 / np.pi
     end_angle = np.arctan2(my - y2, mx - x2) * 180 / np.pi
     # draw, use the first node in the pair to designate color
     if node_pair[0].lower().startswith(('pha', 'fab')) and not node_pair[1].lower().startswith(('pha', 'fab')):
         arcfacecolor=node_position_dict[node_pair[1]]['nodecolor']
     else:
         arcfacecolor=node_position_dict[node_pair[0]]['nodecolor']
     arc = patches.Arc(mxmy, width, height, start_angle, end_angle, color=arcfacecolor, linewidth=arcedgeweights[connection_count], alpha=0.6)
     ax.add_patch(arc)
     connection_count+=1
     #print(x1, y1, x2, y2)
ax.set_xticks(xtickpositions)
ax.set_xticklabels(orderednodes, rotation=90)
ax.tick_params(labelleft=False, axis='both', which='major', pad=-10)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#ax.set_xlim(-0.1, 1.1) # you need to set the appropriate limits explicitly!
ax.set_ylim(-1*(max(markersizes_norm)/2)*.5, ymaxvalue*1.05)
plt.show()
#fig.savefig('Arc digram network co-occurance from 675.png', dpi=750)

######################################################################
#####          CREATE EDGELIST FROM CO-OCCURANCE                 #####
######################################################################
##create a dict containing indexes for all strings and node e.g. GH3=0, GH2=1 etc
index_dict={}
index_dict_r={}
count=0
for c in orderednodes:
    index_dict[c]=count
    index_dict_r[count]=c
    count+=1

pha_gene_list=['PhaA', 'PhaC', 'PhaE', 'PhaG', 'PhaD', 'PhaF', 'FabG', 'PhaR', 'PhaZ', 'PhaC2', 'maoC', 'PhaB',\
     'PhaP', 'PhaJ', 'PhaQ', 'PhaI', 'Pct', 'Pcs', 'PrpE', 'PduP', 'Cat2', 'HadA']

coocc_count_pha_col=coocc_count[np.array(coocc_count.columns)[coocc_count.columns.isin(pha_gene_list)]]
coocc_count_pha_row=coocc_count_pha_col[coocc_count_pha_col.index.isin(pha_gene_list)]
phaorderednodes=[n for n in orderednodes if n in pha_gene_list]

############## OTHER DATA ########################
# the diagram also requires:
    #1. verticies/node order
    #2. vertices/node color
    #3. vertcices size
vertex_data={}
##start from zero if filtering
oncolors=[]
osum=[]
index_dict_new_r={}
newindex=[]
idx=0
for o in phaorderednodes:
    oldindex=index_dict[o]
    newindex.append(idx)
    oncolors.append(orderednodecolors[oldindex])
    osum.append(ordered_sum[oldindex])
    index_dict_new_r[idx]=index_dict_r[oldindex]
    idx+=1
    
vertex_data['order']=newindex  
vertex_data['orderstr']=phaorderednodes
vertex_data['vertexcolor']=oncolors
vertex_data['vertexsize']=[float(n) for n in osum]
index_dict_new=dict((v, k) for k, v in index_dict_new_r.items()) 

with open('40k_pha_only_22g_gt_index.json', 'w') as fp:
    json.dump(index_dict_new, fp)  

with open('40k_pha_only_22g_gt_index_r.json', 'w') as fp:
    json.dump(index_dict_new_r, fp)  

with open('40k_pha_only_22g_gt_vertex_data.json', 'w') as fp:
    json.dump(vertex_data, fp)  

import itertools
# get all possible pairs of (skillid1, skillid2)
edges = list(itertools.combinations(coocc_count_pha_row.columns, 2))  
# find associated weights in the original df
edges_with_weights = [(node1, node2, coocc_count_pha_row.loc[node1][node2]) for (node1, node2) in edges]
# put it all in a new dataframe
edge_df = pd.DataFrame(edges_with_weights, columns=["source", "target", "count"]) 
#filter non zero values
edge_df=edge_df[edge_df['count']>0]
#CONVERT all strings into integer index values from index_dict
edge_df_idx=pd.DataFrame(index=edge_df.index)
edge_df_idx['source']=edge_df['source'].map(index_dict_new)
edge_df_idx['target']=edge_df['target'].map(index_dict_new)
edge_df_idx['count']=edge_df['count']

edge_df_idx.to_csv("40k_22g_Arc_diagram_edgelist_pha_only.csv")



























######################################################################
#####          CORRELATION PLOT WITH DOTS                        #####
######################################################################


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import math

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))

###normalise marker sizes
maxms=3
minms=0.1
coocc_count_clipped=coocc_count.clip(upper=maxms, lower=minms)
ms_multiplier=5
#coocc_count_clipped=coocc_count
normalised_weights=(coocc_count_clipped-coocc_count_clipped.to_numpy().min())/(coocc_count_clipped.to_numpy().max()-coocc_count_clipped.to_numpy().min())
coocc_corr_raw = mergedint.corr(method='pearson')
coocc_corr_trim=coocc_corr_raw.dropna(axis=1, how='all')
coocc_corr_trim=coocc_corr_trim.dropna(axis=0, how='all')
###use this to get rid of the redundant repeats from the original data, and make diagonal plot
mask = np.zeros_like(coocc_corr_trim, dtype=bool)
mask[np.triu_indices_from(mask)] = True
###set max and min values in count, to scale with size of the markers
coocc_corr_trim[mask] = np.nan
coocc_corr_trim=coocc_corr_trim[::-1] # reverse the data frame for a nice diagonal
normalised_weights=normalised_weights[::-1]
normalised_weights_np=normalised_weights.to_numpy() #turn df into numpy array

ylabels = coocc_corr_trim.index
xlabels = coocc_corr_trim.columns
fig, ax = plt.subplots()

###draw grid lines manually to cover only the required area
axcount=len(coocc_corr_trim.columns)
ycount=0
for num in list(range(0, len(coocc_corr_trim.columns))):
    xaxis_x0x1=(num,num)
    xaxis_y0y1=(0,axcount-0.5)
    y_axis_x0x1=(0,axcount-0.5)
    y_axis_y0y1=(ycount,ycount)
    ycount+=1
    axcount+=-1
    ##columns
    plt.plot(xaxis_x0x1, xaxis_y0y1, '-', color=(0.5,0.5,0.5,0.5), zorder=-1,lw=.5)
    ###rows
    plt.plot(y_axis_x0x1, y_axis_y0y1, '-', color=(0.5,0.5,0.5,0.5),zorder=-1,lw=0.5)
    print(y_axis_x0x1, y_axis_y0y1)
   
cooc_corr_cmap=[value for value in coocc_corr_trim.to_numpy().flatten() if not math.isnan(value)]
colormap = cm.coolwarm_r
norm = MidpointNormalize(vmin=np.min(cooc_corr_cmap), vmax=np.max(cooc_corr_cmap), midpoint=0)
#cmap = 'RdBu_r'
s_map = cm.ScalarMappable(norm=norm, cmap=colormap)
##loop over the rows
ypos=0 #row
for row in coocc_corr_trim.index:
    xpos=0 # columns
    for col in coocc_corr_trim.columns:
        if not math.isnan(coocc_corr_trim.loc[row][col])==True:
             #plt.scatter(xpos, ypos, s=normalised_weights.loc[row][col]*3, marker='k')
             ax.plot(xpos, ypos, markersize=normalised_weights.loc[row][col]*ms_multiplier, marker='o',\
                      markerfacecolor=s_map.to_rgba(coocc_corr_trim.loc[row][col]), markeredgecolor=s_map.to_rgba(coocc_corr_trim.loc[row][col]))
             xpos+=1
    ypos+=1
Xsize=len(xlabels)
Ysize=len(ylabels)
ax.set(xticks=np.arange(Xsize), yticks=np.arange(Ysize),
       xticklabels=xlabels, yticklabels=ylabels)
ax.set_xticks(np.arange(Xsize+1)-0.5, minor=True)
ax.set_yticks(np.arange(Ysize+1)-0.5, minor=True)
ax.tick_params(axis='x', which='both', rotation=90, pad=1, labelsize=8)
ax.tick_params(axis='y', which='both', pad=2, labelsize=8)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(-0.5, len(coocc_corr_trim.columns))
ax.set_ylim(-0.5, len(coocc_corr_trim.columns))
cbax=fig.colorbar(s_map, location='left')
cbax.ax.set_ylabel('Pearson correlation', rotation=90, labelpad=-5)
cbax.set_ticks([-0.8,-0.4,0,0.4,0.8])
cbax.set_ticklabels([-0.8, -0.4,0,0.4,0.8])
cbax.ax.tick_params(labelsize=10, pad=0)
cbax.ax.tick_params(size=2)
### add color boxes to delineate what cazy families attack which algae
colheight=0.5
boxwidth = 0.7
boxheight = 0.7
offsetx=boxwidth/2
##plot k2 col_colors
for caz in ax.get_xticklabels():
    a_x=caz.get_position()[0]-offsetx
    a_y=caz.get_position()[1]-colheight-1
    ax.add_patch(Rectangle(
        xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
        linewidth=1, edgecolor='black', lw=1, facecolor=colordict[caz.get_text()], fill=True, clip_on=False))   
### add color boxes to delineate what cazy families attack which algae
colheight=0.5
boxwidth = 0.7
boxheight = 0.7
offsety=boxheight/2
##plot k2 col_colors
for caz in ax.get_yticklabels():
    a_x=caz.get_position()[0]-colheight-1
    a_y=caz.get_position()[1]-offsety
    ax.add_patch(Rectangle(
        xy=(a_x, a_y) ,width=boxwidth, height=boxheight,
        linewidth=1, edgecolor='black', lw=1, facecolor=colordict[caz.get_text()], fill=True, clip_on=False))   
fig.savefig("22g_Pearson corelaton_801_genomes_rowcolcolors.png", dpi=750)
plt.show()






################################################################
###                    Archittecture                        ####
################################################################
from ete3 import NCBITaxa, Tree, TreeStyle, TextFace, SeqMotifFace, CircleFace, RectFace, NodeStyle, faces, AttrFace, ProfileFace, add_face_to_node
from Bio import Phylo
import matplotlib.gridspec 
#import ete3
ncbi = NCBITaxa()


###load data dict with gene positions in
with open('CAZy_data_dict_small.json', 'r') as fp:
    gene_data_dict = json.load(fp)

CAZy_file='Algal_polysaccharide_degrading_CAZy_families.xlsx'  
CAZy_familiy_df=pd.read_excel(CAZy_file, sheet_name='Sheet1', header=0)
###use F6 as it is not filtered
master_PHA_genome_df_F6=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_df_from_nuccore_RefSeqs_biosample_KW_filtered_count_assembly_tax_axn_amp2_n50_M_filtered.csv", sep=',', index_col=0, dtype='unicode')

genomes_of_interest=['NZ_CP030092.1', 'NZ_CP060494.1', 'NZ_MAJC01000002.1', 'NZ_CP012023.1',\
                     'NZ_MIND01000018.1', 'NZ_LT855380.1', 'NZ_CP017009.1', 'NZ_JABWTA010000001.1',\
                     'NZ_CP047025.1', 'NZ_CP031155.1', 'CP011805.1']    
genomes_of_interest=[]
number_uniques=2
for axn in uniquealgaedict:
    for algaetype in uniquealgaedict[axn]:
        if uniquealgaedict[axn][algaetype]>=number_uniques and not axn in genomes_of_interest:
            genomes_of_interest.append(axn)
GOIdf=master_PHA_genome_df_F6[master_PHA_genome_df_F6['Genome_accesion'].isin(genomes_of_interest)]            
####################################
### 40k only genomes of interest ###
####################################
genomes_of_interest=toppct_set            
master_PHA_genome_df_F6=pd.read_csv("22g_2_ALL_taxa_Master_PHA_genome_MERGED_40k_axn.csv", sep=',', index_col=0, dtype='unicode')
GOIdf=master_PHA_genome_df_F6[master_PHA_genome_df_F6['Genome accession'].isin(genomes_of_interest)]
with open('CAZy_data_dict_40k.json', 'r') as fp:
    gene_data_dict = json.load(fp)
####################################
####################################
####################################    
    
##get accessions from gi numbers
gi_to_axn_dict={}
gi_to_axn_dict_r={}
taxidlist=[]
taxid_axn_list={}
#querylist=list(GOIdf['Genome_accesion'])
###40k####
querylist=list(GOIdf['Genome accession'])
handlefetch=Entrez.efetch(db="nuccore", id=querylist, rettype='docsum', retmode='xml')
records=Entrez.read(handlefetch)
axnlist=[]
for idx in list(range(0, len(records))):
    gi_to_axn_dict[records[idx]['Caption']]=records[idx]['Id']
    gi_to_axn_dict_r[records[idx]['Id']]=records[idx]['AccessionVersion']##or caption
    axnlist.append(records[idx]['Id'])
    taxidlist.append(int(records[idx]['TaxId']))
    taxid_axn_list[int(records[idx]['TaxId'])]=records[idx]['AccessionVersion']##or caption
taxid_sciname=ncbi.get_taxid_translator(taxidlist)
taxid_sciname_r=dict((v, k) for k, v in taxid_sciname.items())  
colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':'#669901'}
####get genome sizes
gensize={} 
genomefilepath='22g_High_quality_genome_fastas'
if os.path.exists(genomefilepath)==False:
    os.mkdir(genomefilepath)
maxgensize=0
for axn in axnlist:
    filenm=str(axn)+'.fasta'
    if not filenm in os.listdir(genomefilepath):       
        #dowload it
        print('Fetching genome: ' + str(axn))
        handlefetch=Entrez.efetch(db="nuccore", id=axn, rettype='fasta', retmode='text')   
        print("Fetched... reading...")
        t=SeqIO.read(handlefetch, 'fasta')
        print("Read... writing...")
        SeqIO.write(t,os.path.join(genomefilepath,filenm), 'fasta')
        print('Written... parsing data...')
    genobj=SeqIO.read(os.path.join(genomefilepath, filenm), 'fasta')
    gensize[gi_to_axn_dict_r[axn]]=len(genobj.seq)
    if len(genobj.seq)>maxgensize:
            maxgensize=len(genobj.seq)        
####use accessions to get the taxids to build the phylogenetic tree
architecture_tree = ncbi.get_topology(taxidlist, intermediate_nodes=False, rank_limit=None)

def convert_taxid_to_sciname(intree):
    for node in intree.traverse():
        node.name=str(node.sci_name).replace('. Incertae Sedis','').replace('Family','').replace('=', '_')       
    return intree

taxid_sciname={x:taxid_sciname[x].replace('. Incertae Sedis','').replace('Family','').replace('=', '_') for x in taxid_sciname}
taxid_sciname_r=dict((v, k) for k, v in taxid_sciname.items())  

architecture_tree_sciname=convert_taxid_to_sciname(architecture_tree)
#print(architecture_tree_sciname.get_ascii(attributes=["sci_name", "rank"]))


###save tree in newick format
#architecture_tree_sciname.write(format=3, features=["sci_name", "rank"], outfile="22g_architecture_tree.nw")
###load tree in Phylo 

#tree = Phylo.read("22g_architecture_tree.nw", "newick")  
tree = Phylo.read("22g_40k_architecture_tree.nw", "newick")  
#######################################################################
####get the pha cluster genes start and stop positions for the plot####
#######################################################################
GOI=['phaa', 'phac', 'phae', 'phag', 'phad', 'phaf', 'fabg', 'phar', 'phaz', 'phac2', 'maoc', 'phab',\
     'phap', 'phaj', 'phaq', 'phai', 'pct', 'pcs', 'prpe', 'pdup', 'cat2', 'hada']
final_folder_name='PHA_genomes_for_tree_gbwithparts'
PHA_gene_position_dict_full={}
PHA_gene_position_dict_min={}
cwd=os.getcwd()
#check folders exist to in case script needs re-running
if not os.path.exists(os.path.join(cwd,final_folder_name)):
        os.makedirs(final_folder_name)
for axn in genomes_of_interest:
    axn_handle=str(axn)
    filetype='.gb'
    fn=axn_handle+filetype
    ###if file not in folder, download it
    if not str(axn)+'.gb' in os.listdir(os.path.join(cwd, final_folder_name)):
        rettypestr='gbwithparts'
        retmodestr='text'
        filetype='.gb'
        readas='gb'
        try:
                 print('Fetching genome: ' +str(axn) + str('. UID: ' + str(axn_handle)))
                 handlefetch=Entrez.efetch(db="nuccore", id=axn, rettype=rettypestr, retmode=retmodestr)  
                 print("Fetched... reading...")
                 t=SeqIO.read(handlefetch, readas)
                 print("Read... writing...")
                 SeqIO.write(t,os.path.join(cwd, final_folder_name, fn), readas)
                 print('Written... parsing data...')
        except Exception:
                 print('Error downloadng ' + str(axn))
    elif str(axn)+'.gb' in os.listdir(os.path.join(cwd, final_folder_name)):
        print(str(axn)+'.gb already downloaded and present: extracting info...')
    ###open file in final folder and extract relevant info
    if axn not in PHA_gene_position_dict_full:
        PHA_gene_position_dict_min[axn]={}
        PHA_gene_position_dict_full[axn]={}
    for seq in SeqIO.parse(os.path.join(cwd, final_folder_name, fn), 'gb'):
         for f in seq.features:
         #loop over features in feature table for only genes
              if f.type == "CDS":
              #if a "gene", check for gene name throught the gene key in the subsequent dict
                   if "gene" in f.qualifiers:
                         if f.qualifiers['gene'][0].lower() in GOI:
                             # print(f.qualifiers, f.location)
                              gene=f.qualifiers['gene'][0]
                              locus=f.qualifiers['locus_tag'][0]
                              description=f.qualifiers['product'][0]
                              startstop=f.location
                              start=int(f.location.start)
                              end=int(f.location.end)
                             # print(startstop)
                              #print(gene, description)
                              if gene+'_count' not in PHA_gene_position_dict_full[axn]:
                                   PHA_gene_position_dict_full[axn][gene+'_count']=0
                                   PHA_gene_position_dict_full[axn][gene]=1
                              #give the gene a name like PhaC_0, the next will be PhaC_1 etc
                              gene_num=gene+'_'+str(PHA_gene_position_dict_full[axn][gene+'_count'])
                              PHA_gene_position_dict_full[axn][gene_num]={}
                              PHA_gene_position_dict_full[axn][gene_num]['locus']=locus
                              PHA_gene_position_dict_full[axn][gene_num]['definition']=description
                              PHA_gene_position_dict_full[axn][gene_num]['position']=startstop
                              PHA_gene_position_dict_full[axn][gene+'_count']+=1    
                              PHA_gene_position_dict_min[axn][gene_num]={}
                              PHA_gene_position_dict_min[axn][gene_num]['start']=start
                              PHA_gene_position_dict_min[axn][gene_num]['end']=end
 
                              
pha_col_lower={x.lower():pha_col_colors[x] for x in pha_col_colors}
pha_func_dict_lower={x.lower():pha_function_dict[x] for x in pha_function_dict}

    

t= Tree("22g_architecture_tree.nw", format=3)
t= Tree("22g_40k_architecture_tree.nw", format=3)

##IDS to prune
prunelist=['Vibrio astriarenae', 'Aestuariibacter sp. A3R04'] 
prune=False
if prune==True:
    for name in prunelist:
        n=t.search_nodes(name=name)[0]
        n.delete()

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_scale=False
scaledownfactor=25000
arrowwidthpercent=0.01
maxgenomescaled=int(maxgensize/scaledownfactor)
arrowwidthvalue=int(maxgenomescaled*arrowwidthpercent)
fig_active_phas={}
for node in t.get_leaves():
    nodename=node.name
    ###add genome backbone by utilising the size
    taxid=taxid_sciname_r[nodename]
    axn=taxid_axn_list[taxid]
    genomesize=int(float(gensize[axn])/scaledownfactor)
    #create empy list to store all motifs in
    motifs=[]
    ### draw backbones with width as genome length
    #motifs.append([0,  genomesize, "[]", None, 0, "black", "black", None])
    ###draw cazy families
    # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
    ccount=0
    maxcount=0
    shift=5
    countup=0.9
    ###add the cazy genes to the motifs
    for CDS in gene_data_dict[axn]:    
        overallshift=shift*len(gene_data_dict[axn])
        gene=gene_data_dict[axn][CDS]['Genes'][0]
        startpos, endpos = int(gene_data_dict[axn][CDS]['Positions']['start']), int(gene_data_dict[axn][CDS]['Positions']['stop'])
        startpos = int(startpos/scaledownfactor)
        endpos=int(endpos/scaledownfactor)
        endposmanip=startpos+arrowwidthvalue
        height=12
        tricolor=colrs[CAZy_familiy_df.loc[CAZy_familiy_df['CAZy family'] == gene]['Algae'].iloc[0].lower()]
        motifs.append([int(startpos+(shift*ccount)), int(endpos+(shift*ccount)), "[]", 500, 12, tricolor, tricolor, None])
        if not ccount>=countup*12:
            ccount+=0.6 
    ###add dashed lines every 1mb to the motifs
    one_megabase_scale=1000000/scaledownfactor
    number_dashes=int(genomesize/one_megabase_scale)
    for dash in list(range (0, number_dashes+2)):
        motifs.append([int(dash*one_megabase_scale), int(dash*one_megabase_scale), "[]", 3, 25, '#d8dde6','#d8dde6', None])
    ###add pha genes to the motifs
    phacount=0
    phashift=5
    phawidth=2
    for phagene in PHA_gene_position_dict_min[axn]:
        phagenelower=phagene.split('_')[0].lower()
        phageneupper=phagene.split('_')[0]
        phastart=int(PHA_gene_position_dict_min[axn][phagene]['start']/scaledownfactor)+int(phashift*phacount)
        phaend=int(PHA_gene_position_dict_min[axn][phagene]['end']/scaledownfactor)+int(phashift*phacount)+phawidth
        motifs.append([phastart, phaend, 'o', 0, 20, mpl.colors.to_hex(pha_col_lower[phagenelower]), mpl.colors.to_hex(pha_col_lower[phagenelower]), None])
        phacount+=0.5
        if phagenelower not in fig_active_phas:
            fig_active_phas[phagenelower]=phageneupper
    #motifs.append([0,  genomesize, "[]", None, 2, "grey", "grey", None])
    ### add all motifs to the node
    seqFace = SeqMotifFace(seq=None, motifs=motifs)
    (t & node.name).add_face(seqFace, 0, "aligned")
seqFace = SeqMotifFace(seq=None, motifs=[[0,  100, "[]", None, 0, "grey", "grey", None]])

# Set bold red branch to the root node
style = NodeStyle()
style["fgcolor"] = "#000000"
style["size"] = 0
style["vt_line_color"] = "#000000"
style["hz_line_color"] = "#000000"
style["vt_line_width"] = 2
style["hz_line_width"] = 2
style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
style["hz_line_type"] = 0

t.set_style(style)
for n in t.get_descendants():
    n.img_style=style

t.add_face(seqFace, 0, "aligned")
t.convert_to_ultrametric()

#####write legend
col=0
maxcol=9
##add green/red/brown legends
for y in colrs:
    if col>=maxcol:
        col=0
    ts.legend.add_face(RectFace(8,20,  fgcolor=colrs[y], bgcolor=colrs[y]), column=col)
    col+=1
    ts.legend.add_face(TextFace(y.capitalize()+" acting CAZy"), column=col)
    col+=1 
    ts.legend.add_face(CircleFace(10, "white"), column=col)    
    col+=1
ts.legend_position=3  

plot_funcs_in_legend=True
funcs=[]
for gene in fig_active_phas:
    if plot_funcs_in_legend==True:
         phafuncdict={x.lower():pha_function_dict[x] for x in pha_function_dict}
         func=phafuncdict[gene]
         if not func in funcs:
              if col>=maxcol:
                  col=0
              ts.legend.add_face(CircleFace(8, mpl.colors.to_hex(functiondict[func])), column=col)
              col+=1
              ts.legend.add_face(TextFace(func), column=col)
              col+=1
              #blank col
              ts.legend.add_face(CircleFace(10, "white"), column=col)
              col+=1   
              funcs.append(func)
    
    else:
         if col>=maxcol:
            col=0
         ts.legend.add_face(CircleFace(8, mpl.colors.to_hex(pha_col_lower[gene])), column=col)
         col+=1
         ts.legend.add_face(TextFace(fig_active_phas[gene]), column=col)
         col+=1
         #blank col
         ts.legend.add_face(CircleFace(10, "white"), column=col)
         col+=1
    
t.show(tree_style=ts)
#t.render('40K_phylo_tree_top_0.15pctnile_algae_degraders_func.png', tree_style=ts, w=2400, units="mm")
#t.render('22g_Phylo_tree_top_algae_degraders_2_uniques.png', tree_style=ts, w=2400, units="mm")





# ####start plotting###
# ###gridspec = num rows, num columns, height=[integer heigh for each row], width=[int for each col]
# gs = matplotlib.gridspec.GridSpec(1, 2, height_ratios=[1],
#                              width_ratios=[1, 1], hspace=0, wspace=0) 
# phyl_ax=plt.subplot(gs[0])   
# t=Phylo.draw(tree, axes=phyl_ax, do_show=False)
# arc_ax=plt.subplot(gs[1])  
# #phyl_ax.set_ylim([0-0.5, len(gensize)-0.5])
# #phyl_ax.set_xlim([0, maxgensize*1.025])  
# ###get a list of the leaf terminals to plot arictecture
# terminal_list=tree.get_terminals()
# ###get dictionary of taxid:sci name
# taxid_sciname=ncbi.get_taxid_translator(taxidlist)
# taxid_sciname_r=dict((v, k) for k, v in taxid_sciname.items())
# colrs={'red':'#f03f2e', 'green':'#11ff11', 'brown':'#d2a968', 'brown/green':#669901}
# arrowthickness=0.4
# arrowwidthpercent=0.025 # as a percentage of the width of the total plot, i.e. of maxgensize
# arrowwidthvalue=maxgensize*arrowwidthpercent
# ypos=0.1
# ###use the taxid list, that's ordered from the phylo tree, to get the genome axn in order
# for taxid in taxid_sciname:
#     genaxn=taxid_axn_list[taxid]
#     arc_ax.plot((ypos, gensize[genaxn]), (ypos, ypos), color='black', lw=1) 
#     for CDS in gene_data_dict[genaxn+'.1']:
#         start=float(gene_data_dict[genaxn+'.1'][CDS]['Positions']['start'])
#         gene=gene_data_dict[genaxn+'.1'][CDS]['Genes'][0]
#         if len(gene_data_dict[genaxn+'.1'][CDS]['Genes'])>1:
#             print(genaxn, CDS, gene, '- greater than 1 gene in CDS, handle this exception!\n',gene_data_dict[genaxn+'.1'][CDS]['Genes'])
#         tricolor=colrs[CAZy_familiy_df.loc[CAZy_familiy_df['CAZy family'] == gene]['Algae'].iloc[0].lower()]
#         trianglepos=[[start, ypos-arrowthickness],[start+arrowwidthvalue, ypos],[start, ypos+arrowthickness]] # top x y, middle x y, bottom x y
#         tri = plt.Polygon(trianglepos, color=tricolor)
#         arc_ax.add_patch(tri)
#     ypos+=1
# arc_ax.set_ylim([0-0.5, len(gensize)-0.5])
# arc_ax.set_xlim([0, maxgensize*1.025])    
# arc_ax.tick_params(axis='x', which='both', bottom='on', length=5, width=2, direction='out', pad=2, color='black')    
# arc_ax.set_xlabel("Size (million base pairs)")
# for _, spine in arc_ax.spines.items():
#     spine.set_visible(False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
###Sun et al genomes
filekey='Sun et al supplementary material\Table_4_Phylogenetic Distribution of Polysaccharide-Degrading Enzymes in Marine Bacteria.XLSX'    
SunT4df=pd.read_excel(filekey, sheet_name=0, header=0,engine="openpyxl")    
    
SunIDS=set(['_'.join(p.split('_',2)[:2]) for p in list(SunT4df['Annotated_PDE'])])
print(len(SunIDS))


eg='GCF_000468575.1'
handle=Entrez.elink(dbfrom="refseq", id=eg, db='nuccore', rettype='acc')
record=Entrez.read(handle)

handle=Entrez.esearch(db="assembly", term=eg, retmode="xml")
record=Entrez.read(handle)
print(record)

18310


PhaArecordgen=esearch_database_for_querystring('PhaA', 'AND Bacteria[Orgn] AND genome[Description]', 'nuccore', 100000)
PhaBrecordgen=esearch_database_for_querystring('PhaB', 'AND Bacteria[Orgn] AND genome[Description]', 'nuccore', 100000)
PhaCrecordgen=esearch_database_for_querystring('PhaC', 'AND Bacteria[Orgn] AND genome[Description]', 'nuccore', 100000)

Alist=PhaArecordgen['IdList']
Blist=PhaBrecordgen['IdList']
Clist=PhaCrecordgen['IdList']

for x in Clist:
    if x in Blist:
        print('b')
    if x in Alist:
        print('a')
        


record[0]['GBSeq_definition']
##NB 5 is the len(record[0]['GBSeq_feature-table'])-1 index of that dict....
record[0]['GBSeq_feature-table'][5]['GBFeature_quals'][2]['GBQualifier_value']

handle=Entrez.efetch(db="protein", id=fetchquerystr[0:190], rettype='gb', retmode='xml')
bytehandle = BytesIO(bytes(handle.read(), 'utf-8'))
#record=Entrez.read(handle)
record = Entrez.read(bytehandle)




testterm=' '.join([' '  + str(x) +'[Orgn]' for x in target_taxa_families_gt50])

handlex=Entrez.esearch(db="genome", term="Phac[Gene] AND Acetobacteraceae[Orgn]", retmode="xml", retmax=1500)
recordx=Entrez.read(handlex)
print(len(recordx['IdList']))



with open('outfile.xml', 'w') as outfile:
    outfile.write(handle.read())
with open('outfile.xml', 'rb') as sourcefile:
     record = Entrez.read(sourcefile)




genus='Xanthomonas '
lvl='genus'
a=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches[lvl].fillna('Unknown').str.contains(genus)]
b=a[a['Genome_accesion'].isin(validated_pha_count_df_filt.index)]
axnlist=master_PHA_genome_df_C_C2only_kwmatches.loc[b.index]['Genome_accesion']
c=merged.loc[axnlist]
c=c[c.columns[c.sum(axis=0)>0]]


list_of_tax=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion'].isin(validated_pha_count_df_filt.index)]['genus']

#lg=greenlist
#print('G Mean:' + str(c[[x for x in lg if x in c.columns]].sum(axis=1).mean()) + ' SEM: ' + str(c[[x for x in lg if x in c.columns]].sum(axis=1).sem()) + ' N: ' + str(len(c)))
#lb=brownlist
#print('B Mean:' + str(c[[x for x in lb if x in c.columns]].sum(axis=1).mean()) + ' SEM: ' + str(c[[x for x in lb if x in c.columns]].sum(axis=1).sem()) + ' N: ' + str(len(c)))

lr=redlist
print('R Mean:' + str(c[[x for x in lr if x in c.columns]].sum(axis=1).mean()) + ' SEM: ' + str(c[[x for x in lr if x in c.columns]].sum(axis=1).sem()) + ' N: ' + str(len(c)))


level='genus'
list_of_tax=set(master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion'].isin(validated_pha_count_df_filt.index)][level])
list_of_tax=[g for g in list_of_tax if (type(g)==str)==True]
v=[]
err=[]
tdata={}
for t in list_of_tax:
     if not (type(t)==str) == True:
          if np.isnan(t)==True:
              continue
     genus=t
     a=master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches[level].fillna('Unknown').str.contains(genus)]
     b=a[a['Genome_accesion'].isin(validated_pha_count_df_filt.index)]
     axnlist=master_PHA_genome_df_C_C2only_kwmatches.loc[b.index]['Genome_accesion']
     c=merged.loc[axnlist]
     c=c[c.columns[c.sum(axis=0)>0]]
     tdata[t]=list(c[[x for x in lr if x in c.columns]].columns)
     v.append(c[[x for x in lr if x in c.columns]].sum(axis=1).mean())
     err.append(c[[x for x in lr if x in c.columns]].sum(axis=1).sem())

f, ax = plt.subplots(figsize=(12, 8))    
ax.bar(list(range(0,len(v))), v, yerr=err)
ax.set_xticks(list(range(0,len(v))))
ax.set_xticklabels(list_of_tax, rotation=90, fontsize=8)


master_PHA_genome_df_C_C2only_kwmatches[master_PHA_genome_df_C_C2only_kwmatches['Genome_accesion'].str.contains('NZ_CP005973.1')]['species']
        
        
#cols=sourceaccession, sourceid, linkname protein_nuccor, linkname protein_genome, linkname protein_assembly, taxonomy 
