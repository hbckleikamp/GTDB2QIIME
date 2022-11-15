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

#%% Reps

out_fa=reps_output_base_filename+".fa"
out_txt=reps_output_base_filename+".tsv"

#get taxonomies
ranks=["superkingdom","phylum","class","order","family","genus","species"]
tdf=[]
for t in metadata_filepaths:
    tdf.append(pd.read_csv(t,header=None,sep="\t"))#["Accession"]+ranks)
tdf=pd.concat(tdf)

tdf.columns=["tax_id","lineage"]
tdf[ranks]=tdf["lineage"].str.rsplit(";",expand=True)
tdf=tdf[["tax_id"]+ranks[::-1]]
tdf.to_csv(out_txt,index=False,sep="\t")

#%%


records=[]
for file in reps_input_files:
    records.append(pd.DataFrame([[str(i.id),str(i.seq)] for i in Bio.SeqIO.parse(file,"fasta")],columns=["headers","seqs"]))
records=pd.concat(records)
fl_rec=records[records["seqs"].apply(len)>minimum_length]
fl_tax=tdf[tdf.iloc[:,0].isin(fl_rec["headers"])].apply("\t".join,axis=1).tolist()
fl_rec=fl_rec[["headers","seqs"]].apply("\n".join,axis=1).tolist()

#write to output

if out_fa  in os.listdir(reps_input_folder): os.remove(str(Path(reps_input_folder,out_fa)))
# if out_txt in os.listdir(reps_input_folder): os.remove(str(Path(reps_input_folder,out_txt)))
with open(out_fa,"a+") as fa: #, open(str(Path(reps_input_folder,out_txt)),"a+") as txt:

    for ix in range(len(fl_rec)):
        print(ix)
        fa.write ((">"+fl_rec[ix]+"\n").replace("\x00",""))
        # txt.write((fl_tax[ix]+"\n").replace("\x00","")) #remove bytestring


#%% All
out_fa=all_output_base_filename+".fa"
out_txt=all_output_base_filename+".tsv"

#get taxonomies
lines=[]
with open(file ,"r") as f:
    for line in f.readlines():
        line=line.replace("\x00","")
        
        if line.startswith(">"):
            lines.append(line.split("[")[0].strip().split(" ",1))
tdf=pd.DataFrame(lines,columns=["short_heads","lineage"])




records=pd.DataFrame([[str(i.id),str(i.seq)] for i in Bio.SeqIO.parse(all_input_file,"fasta")],columns=["headers","seqs"])
records["headers"]=">"+records["headers"]
records["short_heads"]=records["headers"].apply(lambda x: x.split("~")[0])


records=records.merge(tdf,on="short_heads",how="inner")

fl_rec=records[records["seqs"].apply(len)>minimum_length]
fl_tax=fl_rec[["headers","lineage"]].apply("\t".join,axis=1).tolist()
fl_rec=fl_rec[["headers","seqs"]].apply("\n".join,axis=1).tolist()

#write to output

if out_fa  in os.listdir(all_input_folder): os.remove(str(Path(all_input_folder,out_fa)))
# if out_txt in os.listdir(all_input_folder): os.remove(str(Path(all_input_folder,out_txt)))
with open(out_fa,"a+") as fa:#, open(str(Path(all_input_folder,out_txt)),"a+") as txt:

    for ix in range(len(fl_rec)):
        print(ix)
        fa.write ((fl_rec[ix]+"\n").replace("\x00",""))
        # txt.write((fl_tax[ix][1:]+"\n").replace("\x00","")) #remove bytestring
        
tdf.columns=["tax_id","lineage"]
tdf["tax_id"]=tdf["tax_id"].str.replace(">","")
tdf[ranks]=tdf["lineage"].str.rsplit(";",expand=True)
tdf=tdf[["tax_id"]+ranks[::-1]]
tdf.to_csv(out_txt,index=False,sep="\t")

