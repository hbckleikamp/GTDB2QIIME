# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:19:05 2021

@author: hbckleikamp
"""


#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()

#%% standard Modules
import pandas as pd
import numpy as np
import Bio
from Bio import SeqIO

#%% Parameters

#ssu reps
reps_input_folder=Path(basedir,"GTDB_ssu_reps")
reps_input_files=[str(Path(reps_input_folder,d)) for d in ["ar53_ssu_reps_r207.fna","bac120_ssu_reps_r207.fna"]]
reps_output_base_filename="GTDB_ssu_reps"

#ssu all
all_input_folder=Path(basedir,"GTDB_ssu_all")
all_input_file=str(Path(all_input_folder,"ssu_all_r207.fna"))
all_output_base_filename="GTDB_ssu_all"

#metadata files used for taxonomic annotation
metadata_filepaths=[str(Path(basedir,"GTDB-metadata",i)) for i in ["bac120_taxonomy.tsv","ar53_taxonomy.tsv"]]

minimum_length=1200 #length trimming of fragmented GTDB ssu sequences 

ranks=["superkingdom","phylum","class","order","family","genus","species"]
#%% Reps

output_folder="GTDB_reps_emu"
path=str(Path(basedir,output_folder))
if not os.path.exists(path): os.mkdir(path)

out_fa=reps_output_base_filename+".fa"
out_txt=reps_output_base_filename+".tsv"

#get taxonomies
tdf=[]
for t in metadata_filepaths:
    tdf.append(pd.read_csv(t,header=None,sep="\t"))#["Accession"]+ranks)
tdf=pd.concat(tdf)

tdf.columns=["accessions","lineage"]
tdf[ranks]=tdf["lineage"].str.rsplit(";",expand=True)

records=[]
for file in reps_input_files:
    records.append(pd.DataFrame([[str(i.id),str(i.seq)] for i in Bio.SeqIO.parse(file,"fasta")],columns=["accessions","seqs"]))
records=pd.concat(records)
fl_rec=records[records["seqs"].apply(len)>minimum_length].drop_duplicates()
fl_rec["tax_id"]=np.arange(1,len(fl_rec)+1)

fasta=("\n".join(">"+fl_rec["tax_id"].astype(str)+":"+fl_rec["accessions"]+"\n"+fl_rec["seqs"])+"\n").replace("\x00","")

#write to output
if out_fa  in os.listdir(reps_input_folder): os.remove(str(Path(reps_input_folder,out_fa)))
with open(str(Path(basedir,output_folder,"species_taxid.fa")),"w") as fa: 
    fa.write(fasta)

tdf=tdf.merge(fl_rec,on="accessions")
tdf=tdf[["tax_id"]+ranks[::-1]]
tdf.to_csv(str(Path(basedir,output_folder,"taxonomy.tsv")),index=False,sep="\t")

#%% All

output_folder="GTDB_all_emu"
path=str(Path(basedir,output_folder))
if not os.path.exists(path): os.mkdir(path)

out_fa=all_output_base_filename+".fa"
out_txt=all_output_base_filename+".tsv"

#get taxonomies
lines=[]
with open(all_input_file ,"r") as f:
    for line in f.readlines():
        line=line.replace("\x00","")
        
        if line.startswith(">"):
            lines.append(line.split("[")[0].strip().split(" ",1))
tdf=pd.DataFrame(lines,columns=["accessions","lineage"])
tdf["accessions"]=tdf["accessions"].str.replace(">","")
tdf[ranks]=tdf["lineage"].str.rsplit(";",expand=True)


records=pd.DataFrame([[str(i.id),str(i.seq)] for i in Bio.SeqIO.parse(all_input_file,"fasta")],columns=["accessions","seqs"])
fl_rec=records[records["seqs"].apply(len)>minimum_length]
fl_rec["tax_id"]=np.arange(1,len(fl_rec)+1)


fasta=("\n".join(">"+fl_rec["tax_id"].astype(str)+":"+fl_rec["accessions"]+"\n"+fl_rec["seqs"])+"\n").replace("\x00","")

#write to output
if out_fa  in os.listdir(reps_input_folder): os.remove(str(Path(reps_input_folder,out_fa)))
with open(str(Path(basedir,output_folder,"species_taxid.fa")),"w") as fa: 
    fa.write(fasta)

tdf=tdf.merge(fl_rec,on="accessions")
tdf=tdf[["tax_id"]+ranks[::-1]]
tdf.to_csv(str(Path(basedir,output_folder,"taxonomy.tsv")),index=False,sep="\t")
