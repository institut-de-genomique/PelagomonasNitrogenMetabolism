# Functional annotation of *P. calceolata* RCC100  
This file indicates the procedure used to generate the functionnal annotation version 2 of *P. calceolata* RCC100 genes.  
- The version 1 is available here : https://doi.org/10.1038/s42003-022-03939-z (Supplementary Data 2).  
- The gene models and proteins are accessible here : https://www.ncbi.nlm.nih.gov/Traces/wgs/CAKKNE01?display=download
- Version 2 of the functional annotation is here : ZENODO.  

## Methods and tools  
For this study, we updated the functional annotation *P. calceolata* genes. We used InterProScan v5.61.93.0 (Quevillon et al., 2005) and the following protein domain databases : Pfam, Gene3D, TIGRfam, SMART, CDD and Gene Onthology. All matches with a pvalue below 1x10-5 were retained. A protein alignment against the NR database (24-08-2023 version) was performed with diamond v2.1.2 (Buchfink et al., 2021). The best match is retained if the evalue is below 1x10-5. We used the HMM search tool KofamKoala v1.3.0 to identify KEGG Orthologues (version of November 2023) (Aramaki et al., 2020). Annotations with an e-value < 1x10−5 and a score above the HMM threshold were retained.  We identified homologies with protein clusters of the Eggnog database with eggnog-mapper tool version 2.1.12 using the very-sensitive mode and diamond aligner (Cantalapiedra et al., 2021; Huerta-Cepas et al., 2019). Finally, for the prediction of protein localization, we used DeepLoc version 2.0 (Thumuluri et al., 2022) and TargetP version 2 in eukaryote mode (Almagro Armenteros et al., 2019). We applied this methodology on the 16,667 gene models, translated in the 6 frames in full and on the protein predicted by Gmove (Dubarry et al., 2016). If functional annotations are identified on several frames, the frame with the best score is retained. This new version of the functional annotation of P. calceolata genes is available on Zenodo.

## Preparation of the protein file to annotate
Unique name ended by _P if it's a predicted protein _1 _2,...or _6 if it is a translation (in full with X instead of "*").

```
perl -ne 'chomp;if($_=~/^>/){$_.="_P"}; 's/\\*//g' ; print "$_\n"' Pelagomonas_calceolata.pep.fa > Pelagomonas_calceolata.pep.rename.fa
```
Translation in the 6 frames in full for all gene models with emboss/6.6.0 (https://emboss.sourceforge.net/)
```
transeq -frame 6 -trim -clean -sequence Pelagomonas_calceolata_gene-models.fa -outseq Pelagomonas_calceolata_gene-models_6frametranslation.fa
cat Pelagomonas_calceolata_gene-models_6frametranslation.fa Pelagomonas_calceolata.pep.rename.fa > Pelago_AllProt.fa
```

## interproscan version 5.61.93.0 
Approximately 4h for 2000 protein sequences with 36 CPUs  

```
mkdir Pelago_interpro
cd Pelago_interpro
interproscan -cpu 36 -t p --tempdir temp-interpro -dp -goterms -i ../Pelago_AllProt.fa -f TSV -o Pelago_AllProt.interpro.tsv
cd ../
```


## Diamond blastp v 2.1.2 with NR database (24_08_2023)
Approximately 2h for 2000 sequences with 24 CPUs

```
mkdir Pelago_diamond/
cd Pelago_diamond/
diamond blastp --query ../Pelago_AllProt.fa --db nr.dmnd --outfmt 6 --evalue 0.00001 --unal 0 --out Pelago_AllProt.diamondnr.tsv --threads 12
#best match selection
sort -uk1,1 Pelago_AllProt.diamondnr.tsv ;done > Pelago_AllProt.diamondnr.bestmatch.tsv
```

## Kofamscan (version 1.3.0) with kegg database (November 2023) 
Approximatively 1h for 100,000 proteins with 12 cpu)  
kegg database : ftp://ftp.genome.jp:21/pub/db/kofam/ko_list.gz  
kegg profiles : ftp://ftp.genome.jp:21/pub/db/kofam/profiles.tar.gz  

```
tar -xvzf profiles.tar.gz
gzip ko_list.gz
mkdir Pelago_kofamscan/
cd Pelago_kofamscan/
mkdir temp-kofam/
exec_annotation -k ko_list -p profiles --cpu 36 --tmp-dir temp-kofam/ --e-value 0.00001 -o Pelago_AllProt.kofamscan.tsv ../Pelago_AllProt.fa
cd ../
```
## EggNog-mapper version 2.1.12
Approximatively 1h for 110,000 proteins with 36 cpu and 100G of memory  

```
mkdir PelagoV4_EggNog
cd PelagoV4_EggNog
emapper.py --dbmem --cpu 36 -m diamond  --sensmode very-sensitive --no_file_comments -i ../Pelago_AllProt.fa -o Pelago_AllProt_eggnog
cd ../
```

## DeepLoc version 2  for protein subcellular localization
Approximatively 1 minute for 10 proteins with 1 CPU and 15G of memory  

```
mkdir PelagoV4_DeepLoc2
cd PelagoV4_DeepLoc2
deeploc2 -m Accurate -f ../Pelago_AllProt.fa -o Pelago_AllProt_DeepLocAccurate.csv;done
cd ../
```
## TargetP version 2 for protein subcellular localization based on signal peptides

```
mkdir PelagoV4_TargetP
cd PelagoV4_TargetP
targetp -fasta Pelago_AllProt.fa -org pl -format short -prefix Pelago_AllProt -batch 1000
cd ../
```
## Concatenation of all results 
perl script in this directory. If functional annotations are identified on several frames, the frame with the best score is retained.
kegg Brite hierarchy is available here https://www.kegg.jp/kegg/download/ 

```
perl FuncAnnotConcatenation.pl -interpro Pelago_interpro/Pelago_AllProt.interpro_full.tsv -nr Pelago_diamond/Pelago_AllProt.diamondnr.bestmatch.tsv -nrDB nr -kofam Pelago_kofamscan/Pelago_AllProt.kofamscan.tsv -eggnog Pelago_EggNog/Pelago_AllProt_eggnog.emapper.annotations -keggHierarchy KO_Network.keg -deeploc Pelago_DeepLoc2/Pelago_AllProt_DeepLocAccurate.csv -targetp Pelago_TargetP/Pelago_AllProt_summary.targetp2 -out Pelago_AllProt.annotConcat
```
