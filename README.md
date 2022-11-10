# GTDB2QIIME

Set of simple auxiliary python scripts to help create GTDB databases for annotation with QIIME

<br>
This collection of scripts is designed to facillitate QIIME annotations with GTDB reps ssu, or GTDB all ssu 
It consists of 2 separate small scripts, that can be run in Spyder, or reused for other purposes.
Since QIIME varies in dependancies and requirements on different operating systems, automated Qiime installation is not included,
and should be done following: https://docs.qiime2.org/2021.8/
<br><br>
Running a pipeline would consist of:<br>
-1. GTDB_ssu_download.py: to download recent GTDB ssu fasta files and taxonomy metadata<br>
-2. GTDB_ssu_prep.py: to format GTDB ssu into seperate taxonomy and fasta files, with option to trim fragmented sequences below certain basepair length <br>



#### Licensing

The pipeline is licensed with standard MIT-license. <br>
If you would like to use this pipeline in your research, please cite the following papers: 
      
-Bolyen E, et al. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852–857. https://doi.org/10.1038/s41587-019-0209-9
<br>
-Parks, D.H., et al. 2020. A complete domain-to-species taxonomy for Bacteria and Archaea. Nature Biotechnology, https://doi.org/10.1038/s41587-020-0501-8.


#### Contact:
-Hugo Kleimamp (Developer): hbckl@bio.aau.dk <br> 
-Martin Pabst: M.Pabst@tudelft.nl<br>


#### Recommended links to other repositories:
https://github.com/qiime2/qiime2
