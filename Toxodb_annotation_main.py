#!/usr/bin/env python3


#Name: Pablo Angulo Lara
#Program: Toxodb_annotationV3

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Program libraries

import requests  
import re        
import os       
import pickle
import _pickle as cPickle
import bz2
import sys
import collections
import pandas as pd
import statistics as stats
from multiprocessing import Pool

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#PARAMETERS

feature = "gene"
ref_genome = "ME49_53"
Ref_genome_isolate = f'Ref {ref_genome}'

create_newdb = "FALSE"
region_file = f'{ref_genome}_bins_1KB.bed'
gene_file = f'Toxodb_{ref_genome}.txt'
chr_file = f'{ref_genome}.bed'

genedb_file = f'{ref_genome}_genedb_gff3_NSsites.pbz2'
regiondb_file = f'{ref_genome}_regiondb_1KB.pbz2'

genes_to_regiondb = "FALSE"

gff3_info = "FALSE"
gff3_file = f'{ref_genome}_gff3.txt'

N_S_sites_info = "FALSE"
n_s_sites_file = "N_S_sites.txt"

snp_info = "FALSE"
snp_file = "Type_II_filter.table"

cnv_info = "TRUE"
cnv_dir = "Gene_CNV/"

newdb_path = f'{ref_genome}_{feature}db_gff3_NSsites_CNV.pbz2'
save_database = "TRUE"
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#CLASSES

class Gene:
    def __init__(self,params):
        self.id = params.get('id',"__Unknown__")
        self.name = params.get('name',"__Unknown__")
        self.organism = params.get('organism',"__Unknown__")
        self.product = params.get('product',"__Unknown__")
        self.chromosome = params.get('chromosome',"__Unknown__")
        self.chromosome_length = params.get('chromosome_length',"__Unknown__")
        self.start = params.get('start',"__Unknown__")
        self.end = params.get('end',"__Unknown__")
        self.strand = params.get('strand',"__Unknown__")
        self.type = params.get('type',"__Unknown__")
        self.transcript_id = params.get('transcript_id',"__Unknown__")
        self.n_transcripts = params.get('n_transcripts',"__Unknown__")
        self.n_exons = params.get('n_exons',"__Unknown__")
        self.utr_length = params.get('utr_length',dict())
        self.transcript_length = params.get('transcript_length',"__Unknown__")
        self.cds_length = params.get('cds_length',"__Unknown__")
        self.apollo_product_description = params.get('apollo_product_description',"__Unknown__")          
        self.protein = params.get('protein',dict())
        self.GO = params.get('GO',dict())
        self.ortholog = params.get('ortholog',dict())
        self.paralog = params.get('paralog',dict())
        self.N_sites = params.get('N_sites',"__Unknown__")
        self.S_sites = params.get('S_sites',"__Unknown__")
        self.snp = params.get('snp',dict())
        self.cnv = params.get('cnv',dict())
        self.depth = params.get('depth',dict())
        self.gff3 = params.get('gff3',dict())
        self.N_sites = params.get('N_sites',"__Unknown__")
        self.S_sites = params.get('S_sites',"__Unknown__")
        self.pmids = params.get('pmids',collections.defaultdict(list))
    def import_cnv(self,df,isolate):
        df.columns = ['Gene_id','Chromosome','Estimated_ploidy','TPM','Gene_dose','Haploid_number','Ortholog_id',
                      'Copies_in_reference','Cluster_gene_dose','Cluster_haploid_number','Gene_list']
        df = df[df.Gene_id == self.id]
        df.reset_index(inplace=True,drop=True)
        if not df.empty:
            df = df.iloc[0,:]
            self.cnv[isolate] = df.to_dict()
        else:
            self.cnv[isolate] = "No CNV info"
    def import_gff3(self,gff3_df):
        df = gff3_df[gff3_df.ID.str.contains(self.id)]
        for row in zip(df['seqid'],df['type'],df['start'],df['end'],df['strand'],df['phase'],df['ID'],
                       df['description'],df['Parent'],df['protein_source_id'],df['Name'],df['Note']):
            if row[0] == self.chromosome and row[2] >= self.start and row[3] <= self.end:
                self.gff3[row[6]] = {'ID':row[6],'Type':row[1],'Chromosome':row[0],'Start':int(row[2]),'End':int(row[3]),'Gene_id':self.id,'Parent':row[8],'Strand':row[4],
                                         'Phase':row[5],'Description':[7],'Protein_source_id':row[9],'Name':row[10],'Note':row[11]}
    def import_N_S_sites(self,df):
        df = df[df['product'] == self.id]
        df.reset_index(inplace=True,drop=True)
        if not df.empty:
            self.S_sites = float(df.loc[0,'S_sites'])
            self.N_sites = float(df.loc[0,'N_sites'])
        else:
            self.S_sites = 'None'
            self.N_sites = 'None'     
    
class Region:
    def __init__(self,params):
        self.id = params.get('id',"__Unknown__")
        self.genes = params.get('genes',dict())
        self.chromosome = params.get('chromosome',"__Unknown__")
        self.chromosome_length = params.get('chromosome_length',"__Unknown__")
        self.start = params.get('start',"__Unknown__")
        self.end = params.get('end',"__Unknown__")
        self.snp = params.get('snp',dict())
        self.depth = params.get('depth',dict())
        self.cnv = params.get('cnv',collections.defaultdict(list))
        self.gff3 = params.get('gff3',dict())
    def import_genes(self):
        for gene_id,gene in genedb.items():
            if gene.chromosome == self.chromosome and gene.start >= self.start and gene.start <= self.end:
                self.genes[gene_id] = gene
            if gene.chromosome == self.chromosome and gene.end >= self.start and gene.start <= self.end:
                self.genes[gene_id] = gene
    def import_gff3(self,df):
        df = gff3_df[(gff3_df.seqid == self.chromosome) & (((gff3_df.start >= self.start) & (gff3_df.start <= self.end)) | ((gff3_df.end >= self.start) & (gff3_df.end <= self.end)))]
        for row in zip(df['seqid'],df['type'],df['start'],df['end'],df['strand'],df['phase'],df['ID'],
                       df['description'],df['Parent'],df['protein_source_id'],df['Name'],df['Note']):
            self.gff3[row[6]] = {'ID':row[6],'Type':row[1],'Chromosome':row[0],'Start':int(row[2]),'End':int(row[3]),'Gene_id':self.id,
                                 'Parent':row[8],'Strand':row[4],'Phase':row[5],'Description':[7],'Protein_source_id':row[9],'Name':row[10],'Note':row[11]}
            
                   
def import_snp(id):
    snp_dict = {}
    df = snp_df[(snp_df.CHROM == db[id].chromosome) & (snp_df.POS >= db[id].start) & (snp_df.POS <= db[id].end)]
    df.reset_index(inplace=True,drop=True)
    for isolate in sample_list:
        snp_dict[isolate] = dict()
        gff3_features = ['exon','CDS','five_prime_UTR','three_prime_UTR',feature]
        features = list(set([re.sub('\w*UTR','UTR',i['Type']) for i in db[id].gff3.values() if i['Type'] in gff3_features]))
        features.append(feature)
        for feat in gff3_features:
            feat = re.sub('\w*UTR','UTR',feat)
            if feat in features:
                snp_dict[isolate][feat] = {}
                snp_dict[isolate][feat]['count'] = 0
                snp_dict[isolate][feat]['vector'] = []
                snp_dict[isolate][feat]['sequence'] = []
                snp_dict[isolate][feat]['general_info'] = {}
                if isolate == Ref_genome_isolate:
                    snp_dict[isolate][feat]['total_count'] = 0
            else:
                snp_dict[isolate][feat] = 'No Features'
        if isolate == Ref_genome_isolate:   
            for row in zip(df.index,df['ID'],df['CHROM'],df['POS'],df['TYPE'],df['REF'],df['ALT'],df['ANN']):
                allele_list = df.loc[row[0],genotype_cols].to_list()
                allele_list = [x for x in allele_list if x != '.']
                allele_freqs = collections.Counter(allele_list).most_common()
                ref_allele = row[5]
                if row[5] == '*':
                    ref_allele = '-'
                max_allele_length = max(df.loc[row[0],genotype_cols].str.len())
                if len(ref_allele) != max_allele_length:
                    ref_allele += '-' * (max_allele_length - len(ref_allele))
                alt_alleles = row[6].split(',')
                num_alleles = len(set(list(allele_list)))
                snpeff = re.findall(f'\|(\w*variant)',row[7])
                snp_dict[isolate][feature]['general_info'][row[1]] = {"Chromosome":row[2],"Position":row[3],"Type":row[4],"Ref":row[5],"Alt":row[6],
                                                                    "Allele_count":num_alleles,'Feature_info':[],'Allele_mod':ref_allele,
                                                                    'Major_allele':allele_freqs[0],'Minor_allele':allele_freqs[-1],'Allele_freqs':allele_freqs,'SnpEff':snpeff}                 
                snp_dict[isolate][feature]['total_count'] += 1
                snp_dict[isolate][feature]['sequence'].append(ref_allele)
                snp_dict[isolate][feature]['vector'].append(0)
                for gff3_feat in db[id].gff3.values():
                    if row[3] >= gff3_feat['Start'] and row[3] <= gff3_feat['End'] and gff3_feat['Type'] in gff3_features:
                        gff3_feat_type = re.sub('\w*UTR','UTR',gff3_feat['Type'])
                        snp_dict[isolate][gff3_feat_type]['sequence'].append(ref_allele)
                        snp_dict[isolate][gff3_feat_type]['vector'].append(0)
                        snp_dict[isolate][gff3_feat_type]['total_count'] += 1
                        snp_dict[isolate][feature]['general_info'][row[1]]['Feature_info'].append(gff3_feat)
                        snp_dict[isolate][gff3_feat_type]['general_info'][row[1]] = snp_dict[isolate][feature]['general_info'][row[1]]
        else:
            isolate_colname = f'{isolate}.GT'
            for row in zip(df.index,df['ID'],df['CHROM'],df['POS'],df['TYPE'],df['REF'],df[isolate_colname],df['ALT'],df['ANN']):
                allele_list = df.loc[row[0],genotype_cols].to_list()
                allele_list = [x for x in allele_list if x != '.']
                allele_freqs = collections.Counter(allele_list).most_common()
                major_allele = allele_freqs[0][0]
                allele_dict = {row[5]:0}
                allele = row[6]
                snpeff_allele = row[6]
                if row[6] == '.':
                    allele = major_allele
                    snpeff_allele = 'X'
                if allele == '*':
                    allele = '-'
                    snpeff_allele = 'X'
                max_allele_length = max(df.loc[row[0],genotype_cols].str.len())
                if len(allele) != max_allele_length:
                    allele += '-' * (max_allele_length - len(allele))
                alt_alleles = row[7].split(',')
                count = 1
                for a in alt_alleles:
                    allele_dict[a] = count
                    count += 1
                allele_dict['.'] = allele_dict[major_allele]
                num_alleles = len(set(list(allele_list)))
                snpeff = re.findall(f'{snpeff_allele}\|(\w*variant)',row[8])
                snp_dict[isolate][feature]['sequence'].append(allele)
                snp_dict[isolate][feature]['vector'].append(allele_dict[row[6]])
                snp_dict[isolate][feature]['general_info'][row[1]] = {"Chromosome":row[2],"Position":row[3],"Type":row[4],'Variant':'REF',"Ref":row[5],"Allele":row[6],"Alt":row[7],
                                                                    "Allele_count":num_alleles,'Feature_info':[],'Allele_mod':allele,
                                                                    'Allele_vector':allele_dict[row[6]],'Major_allele':allele_freqs[0],'Minor_allele':allele_freqs[-1],
                                                                    'Allele_freqs':allele_freqs,'SnpEff':snpeff}                                                            
                if row[5] != row[6] and row[6] != '.':
                    snp_dict[isolate][feature]['count'] += 1
                    snp_dict[isolate][feature]['general_info'][row[1]]['Variant'] = 'ALT'
                if row[6] == '.':
                    snp_dict[isolate][feature]['general_info'][row[1]]['Variant'] = 'Missing'
                for gff3_feat in db[id].gff3.values():
                    if row[3] >= gff3_feat['Start'] and row[3] <= gff3_feat['End'] and gff3_feat['Type'] in gff3_features:
                        gff3_feat_type = re.sub('\w*UTR','UTR',gff3_feat['Type'])
                        snp_dict[isolate][gff3_feat_type]['sequence'].append(allele)
                        snp_dict[isolate][gff3_feat_type]['vector'].append(allele_dict[row[6]])
                        snp_dict[isolate][feature]['general_info'][row[1]]['Feature_info'].append(gff3_feat)
                        snp_dict[isolate][gff3_feat_type]['general_info'][row[1]] = snp_dict[isolate][feature]['general_info'][row[1]]
                        if row[5] != row[6] and row[6] != '.':
                            snp_dict[isolate][gff3_feat_type]['count'] += 1
    print(id)
    new_snp = [id,snp_dict]
    return(new_snp)
                                      

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#FUNCTIONS

def run_multiprocessing(func,input_list):
    p = Pool(processes=8)
    out = p.imap_unordered(func,input_list)
    return out
  
def save_db(db,db_path):
    db_file = bz2.BZ2File(db_path,'w')
    cPickle.dump(db,db_file)
    db_file.close()

def load_db(db_path):
    db_file = bz2.BZ2File(db_path,'rb')
    db = cPickle.load(db_file)
    db_file.close()
    return db
    
def cnv_processing(file):
    isolate = file.rsplit('_geneCNVs.tsv')[0]
    isolate = isolate.rsplit('/')[-1]
    df = pd.read_table(file)
    for id in db:
        db[id].import_cnv(df,isolate)
    cnv_out = [isolate,db]
    print(f'ISOLATE {isolate} PROCESSED')
    return cnv_out
    
def gene_attributes(gene_id):
    df = toxodb_df[toxodb_df['Gene ID'] == gene_id]
    gene = df.to_dict(orient='records')[0]
    for key, value in gene.items():
        value = str(value)
        gene[key] = value.split(';')
    pos = gene['Genomic Location (Gene)'][0]
    start = int(re.findall(':([0-9,]*)',pos)[0].replace(',',''))
    end = int(re.findall('\.\.([0-9,]*)',pos)[0].replace(',',''))
    gene['start'] = start
    gene['end'] = end
    chr = gene['Genomic Sequence ID'][0]
    chr_length = int(chr_df.loc[chr_df.Chromosome == chr,'End'].values[0])
    try:
        gene_dict = {'id':gene['Gene ID'][0],'organism':gene['Organism'][0],
             'product':gene['Product Description'][0],'n_exons':int(gene['# Exons in Gene'][0]),
             'n_transcripts':int(gene['# Transcripts'][0]),'start':int(gene['start']),'end':int(gene['end']),
             'utr_length':{'3_utr_length':float(gene["Annotated 3' UTR length"][0]),
                           '5_utr_length':float(gene["Annotated 5' UTR length"][0])},
             'strand':gene['Gene Strand'][0],'type':gene['Gene Type'][0],
             'transcript_length':int(gene['Transcript Length'][0]),
             'apollo_product_description':gene['Apollo Product Description'][0],
             'name':gene['Gene Name or Symbol'][0],'chromosome':chr,'chromosome_length':chr_length,
             'ortholog':{'count':int(gene['Ortholog count'][0]),'group':gene['Ortholog Group'][0]},
             'paralog':{'count':int(gene['Paralog count'][0])},'cds_length':float(gene['CDS Length'][0]),
             'protein':{'Length':float(gene['Protein Length'][0]),
                        'Molecular weight':float(gene['Molecular Weight'][0]),
                        'Isoelectric point':float(gene['Isoelectric Point'][0]),
                        'Interpro desc':gene['Interpro Description'],'PFam desc':gene['PFam Description'],
                        'PirSF desc':gene['PirSF Description'],
                        'Prositefamilies desc':gene['Prositefamilies Description'],
                        'Smart desc':gene['Smart Description'],
                        'Superfamily desc':gene['Superfamily Description'],
                        'TigrFam desc':gene['TigrFam Description'],
                        'Predicted Location (TAGM-MAP)':gene['Predicted Location (TAGM-MAP)'][0],
                        'Top Predicted Location (TAGM-MCMC)':gene['Top Predicted Location (TAGM-MCMC)'][0],
                        'TM Domains':float(gene['# TM Domains'][0]),'EC numbers':gene['EC numbers'],
                        'EC numbers from OrthoMCL':gene['EC numbers from OrthoMCL']},
             'GO':{'Computed Component':{gene['Computed GO Component IDs'][i]:gene['Computed GO Components'][i] for i in range(len(gene['Computed GO Component IDs'])) if i < len(gene['Computed GO Components'])},
                   'Computed Function':{gene['Computed GO Function IDs'][i]:gene['Computed GO Functions'][i] for i in range(len(gene['Computed GO Function IDs'])) if i < len(gene['Computed GO Functions'])},
                   'Computed Process':{gene['Computed GO Process IDs'][i]:gene['Computed GO Processes'][i] for i in range(len(gene['Computed GO Process IDs'])) if i < len(gene['Computed GO Processes'])},
                   'Curated Component':{gene['Curated GO Component IDs'][i]:gene['Curated GO Components'][i] for i in range(len(gene['Curated GO Component IDs'])) if i < len(gene['Curated GO Components'])},
                   'Curated Function':{gene['Curated GO Function IDs'][i]:gene['Curated GO Functions'][i] for i in range(len(gene['Curated GO Function IDs'])) if i < len(gene['Curated GO Functions'])},
                   'Curated Process':{gene['Curated GO Process IDs'][i]:gene['Curated GO Processes'][i] for i in range(len(gene['Curated GO Process IDs'])) if i < len(gene['Curated GO Processes'])}}} 
    except:
        print(gene_id)
    new_gene = [gene_id,gene_dict]
    return new_gene

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#MAIN CODE

if create_newdb == "TRUE" and feature == "gene":
    print(f"CREATING {ref_genome} GENE DATABASE")
    chr_df = pd.read_table(chr_file,header=None)
    chr_df.columns = ['Chromosome','Start','End']
    gene_list = []
    db = dict()
    toxodb_df = pd.read_table(gene_file)
    gene_list = toxodb_df['Gene ID'].tolist()
    res = run_multiprocessing(gene_attributes,gene_list)
    for new_gene in res:
        db[new_gene[0]] = Gene(new_gene[1])

if create_newdb == "TRUE" and feature == "region":
    print(f"CREATING {ref_genome} REGION DATABASE")
    chr_df = pd.read_table(chr_file,header=None)
    chr_df.columns = ['Chromosome','Start','End']
    region_df = pd.read_table(region_file,header=None)
    region_df.columns = ['Chromosome','Start','End']
    db = dict()
    count = 0 
    for row in zip(region_df['Chromosome'],region_df['Start'],region_df['End']):
        region_id = f'{row[0]}_{row[1]}_{row[2]}'
        chr_length = int(chr_df.loc[chr_df.Chromosome == row[0],'End'].values[0])
        if int(row[2]) < chr_length:
            end = int(row[2]) - 1
        if int(row[2]) == chr_length:
            end = int(row[2])
        new_region = {'id':region_id,'chromosome':row[0],'start':int(row[1]),'end':end,'chromosome_length':chr_length}
        db[region_id] = Region(new_region)
        count += 1
        print(count)
      
if create_newdb == "FALSE":
    if feature == 'region':
        db = load_db(regiondb_file)                
    if feature == 'gene':
        db = load_db(genedb_file)

if gff3_info == "TRUE":
    if feature == 'gene':
        print(f'SAVING GFF3 INFO INTO GENE DATABASE:\n')
    if feature == 'region':
        print(f'SAVING GFF3 INFO INTO REGION DATABASE:\n')
    gff3_df = pd.read_table(gff3_file)  
    count = 0
    for id in db:
        db[id].import_gff3(gff3_df)
        count += 1
        print(f'{count}: {id}')

if N_S_sites_info == "TRUE":
    print(f"SAVING N_S sites INFO INTO GENE DATABASE ")
    df = pd.read_table(n_s_sites_file)
    count = 0
    for id in db:
        db[id].import_N_S_sites(df)
        count += 1
        print(f'{count}: {id}')      

if genes_to_regiondb == "TRUE":
    print(f'SAVING GENES INTO REGION DATABASE:\n')
    count= 0
    genedb = load_db(genedb_file)
    for id in db:
        db[id].import_genes()
        count += 1
        print(f'{count}: {id}')
     
if snp_info == "TRUE":
    if feature == 'gene':
        print(f'SAVING SNP VARIANTS INTO GENE DATABASE:\n')
    if feature == 'region':
        print(f'SAVING SNP VARIANTS INTO REGION DATABASE:\n')
    snp_df = pd.read_table(snp_file)
    genotype_cols = snp_df.filter(regex='.GT',axis=1).columns
    genotype_cols = genotype_cols.to_list()
    sample_list = [i.replace('.GT','') for i in genotype_cols]
    #genotype_cols.append('REF')
    sample_list.append(Ref_genome_isolate)
    id_list = db.keys()
    res = run_multiprocessing(import_snp,id_list)
    for new_snp in res:
        db[new_snp[0]].snp = new_snp[1]
    
            
if cnv_info == "TRUE":
    if feature == 'gene':
        print(f'SAVING CNV INFO INTO GENE DATABASE:\n')
    if feature == 'region':
        print(f'SAVING CNV INFO INTO REGION DATABASE:\n')
    file_list = []
    for file in os.listdir(cnv_dir):
        if file.endswith("geneCNVs.tsv"):
            file_list.append((cnv_dir+file))
    res = run_multiprocessing(cnv_processing,file_list)
    for out in res:
        for id in db:
            db[id].cnv[out[0]] = out[1][id].cnv[out[0]]

if save_database == "TRUE":
    save_db(db,newdb_path)
  
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
