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
            self.cnv[isolate] = {}
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
            

