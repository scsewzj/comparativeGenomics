#### Check the different types of sequences in the gene assembly
	
```bash
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 1 -d ':' | sort | uniq
# chromosome
# superscaffold

# Check for isoforms in the proteome
```bash
grep -i "variant" filtered_pep.fasta | head
grep -i "isoform" filtered_pep.fasta
```

# Each isoform comes from a different gene that doesn't have another occurence in the proteome file 
```bash
grep "gene:OB04G16260" Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa
grep "gene:OB03G41510" Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa
grep "gene:OB03G16230" Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa

--------
grep "^>" filtered_pep.fasta | grep "mitochond"
# OB12G24960.1 pep chromosome:Oryza_brachyantha.v1.4b:12:13826905:13832854:1 gene:OB12G24960 transcript:OB12G24960.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:methylcrotonyl-CoA carboxylase alpha chain, mitochondrial / 3-methylcrotonyl-CoA carboxylase 1 (MCCA) [Source:Projected from Arabidopsis thaliana (AT1G03090) TAIR;Acc:AT1G03090]
grep "^>" filtered_pep.fasta | grep "mitochond" | wc -l
# 36
grep "^>" filtered_pep.fasta | grep "Source:Projected" | wc -l
# 3601
```

# Created the filtered fasta (named "filtered_pep.fasta")
```bash
./filter_fasta.sh
```

## 2. BLASTP
# Prepare the Protein DB
```bash
makeblastdb -in filtered_pep.fasta -dbtype prot -out filtered_proteins_db

# Building a new DB, current time: 12/03/2025 16:21:47
# New DB name:   /Project_Comparative_Genomics/Brachyantha/filtered_proteins_db
# New DB title:  filtered_pep.fasta
# Sequence type: Protein
# Keep MBits: T
# Maximum file size: 3000000000B
# Adding sequences from FASTA; added 31355 sequences in 2.80157 seconds.
```

# Run the BLASTP
```bash
blastp -query filtered_pep.fasta -db filtered_proteins_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 10 -out all_vs_all.tsv
```

# Move db parameters into a directory (filtered_proteins_db)
```bash
mkdir filtered_proteins_db; mv -i filtered_proteins_db.* filtered_proteins_db/
```

# There are multiple hits:
```bash
awk ' {if ($1>$2){print $2,$1}; if ($2>$1){print $1,$2} } ' all_vs_all.tsv | sort -d | uniq -c | sort -nr -k1 | head -25
```

# Among them, some self hits. We remove them:
```bash
awk '$1 != $2' all_vs_all.tsv > no_self_hits.tsv
```

We then add the length of each sequence:
Let's get the lengths:
```bash
awk '/^>/ {if (seq != "") print id, length(seq);
                split($0, a, " ");      # split header by spaces
                id = substr(a[1], 2);   # remove ">" and keep AccessionID
                seq = ""; next}{seq = seq $0} END {print id, length(seq)
        }' filtered_pep.fasta > lengths_pep.tsv
```

- First, we add the length:
```bash
awk 'NR==FNR {len[$1]=$2; next}
        {print $0, len[$1], len[$2] }' lengths_pep.tsv no_self_hits.tsv > blast_with_lengths.tsv
```

- We keep the longest alignment.:
```bash
awk '{ if ($1 < $2) pair=$1 FS $2; else pair=$2 FS $1;
        if (pair in maxlen) {if ($4 > maxlen[pair]) { maxlen[pair]=$4; best[pair]=$0 }
            } else {maxlen[pair]=$4; best[pair]=$0}
        } END {for (p in best) print best[p]
    }' blast_with_lengths.tsv > treated_blast.tsv
```

# Filterings
## Strict
* % identity ≥ 40
* Query coverage ≥ 70%
* Subject coverage ≥ 70%
```bash
awk '{ qcov = ($8-$7+1)/$13*100; scov = ($10-$9+1)/$14*100; if ($3>=40 && qcov>=70 && scov>=70) print }' treated_blast.tsv > strict.tsv
```

## Moderate homologs (more permissive)
* % identity ≥ 30
* Query coverage ≥ 50%
* *Subject coverage ≥ 50%
```bash
awk '{ qcov = ($8-$7+1)/$13*100; scov = ($10-$9+1)/$14*100; if ($3>=30 && qcov>=50 && scov>=50) print }' treated_blast.tsv > moderate.tsv
```

## Very relaxed homologs (detect remote similarities)
* % identity ≥ 20
* Query coverage ≥ 30% OR Subject coverage ≥ 30%
* E-value ≤ 1e-5
```bash
awk '{ qcov = ($8-$7+1)/$13*100; scov = ($10-$9+1)/$14*100; if ($3>=20 && (qcov>=30 || scov>=30) && $11<=1e-5) print }' treated_blast.tsv > relaxed.tsv
```

# Prepare for Clustering:
```bash
awk '{print $1,"\t",$2,"\t",$12}' treated_blast.tsv > clustering_O_brachyantha_raw.tsv # BEGIN{print "geneIDA\tgeneIDB\tbitScore"}

awk '{print $1,"\t",$2,"\t",$12}' moderate.tsv > clustering_O_brachyantha_moderate.tsv
awk '{print $1,"\t",$2,"\t",$12}' strict.tsv > clustering_O_brachyantha_strict.tsv
awk '{print $1,"\t",$2,"\t",$12}' relaxed.tsv > clustering_O_brachyantha_relaxed.tsv
```

The clustering was made with the use on [*Galaxy server*](https://usegalaxy.eu/root?tool_id=mcl)
```bash
mv -i ../../../Galaxy23-\[MCL\ on\ data\ 19\].tabular clusters/MCL_OryzaB_strict.tabular
mv -i ../../../Galaxy24-\[MCL\ on\ data\ 20\].tabular clusters/MCL_OryzaB_moderate.tabular
mv -i ../../../Galaxy29-\[MCL\ on\ data\ 27\].tabular clusters/MCL_OryzaB_raw.tabular
mv -i ../../../Galaxy30-\[MCL\ on\ data\ 28\].tabular clusters/MCL_OryzaB_relaxed.tabula
```

# TBA
Worked on `display_cluster.R`

# Make a proteome file without the filtered out files 
```bash
./filter_fasta.sh Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa filtered_cds.fasta

python prepare_kaks_input.py clusters/MCL_OryzaB_strict.tabular filtered_cds.fasta filtered_pep.fasta ks_filtered_input

for d in ks_filtered_input/*; do for f in $d/*_prot.fasta; do a="${f%_prot.fasta}"; clustalw2 -quiet -align -infile=$f -outfile=$a.ali.aln; done; echo "$d - DONE"; done;
```

```bash
cp -ir ks_filtered_input/OB01G10020.1_OB01G43890.1 .
./execute_yn00.sh_generalize.sh
rm -r ks_filtered_input/OB01G10020.1_OB01G43890.1
mv -i OB01G10020.1_OB01G43890.1 ks_filtered_input/
```

```bash
for n in ks_filtered_input/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ./pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy;   done; done
```

```bash
./execute_yn00.sh_generalize.sh
```

```bash
python collect_KS.py

cut -f 1,2 strict.tsv | sort -k1 > duplicates_strict.tsv
cut -f 1,2 treated_blast.tsv | sort -k1 > duplicates_raw.tsv
```

```bash
python detect_sigletons.py -raw duplicates_raw.tsv -reference duplicates_strict.tsv -out_duplicates duplicates_strict.txt -out_singletons singletons_strict.txt -delimiter " "
# Found 19543 duplicates and 72597 singletons.
# Duplicates written to: duplicates_strict.txt
# Singletons written to: singletons_strict.txt

python detect_sigletons.py -raw duplicates_raw.tsv -reference duplicates_strict.tsv -out_duplicates duplicates_strict_unique.txt -out_singletons singletons_strict_unique.txt
# Found 110176 duplicates and 74104 singletons.
# Duplicates written to: duplicates_strict_unique.txt
# Singletons written to: singletons_strict_unique.txt

mkdir -p singletons duplicates
mv -i duplicates_strict.txt duplicates_strict_unique.txt duplicates
mv -i singletons_strict.txt singletons_strict_unique.txt singletons
```

# <TBA>
Worked on "detect_Ks_peaks.R"


# 4 Dec. 2025 ---
---

# 21 Dec. 2025

---
# 27 dec. 2025

## Pre-analysis of sequence

#### Number of sequences
```bash
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | sed 's/>//g' | cut -f 1 -d ' '| sort | uniq | wc -l
# 32037
```

#### Number of *pep*
```bash
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | sed 's/>//g' | cut -f 2 -d ' '| sort | uniq | wc -l
# 1

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | sed 's/>//g' | cut -f 2 -d ' '| sort | uniq
# pep
```

#### Check the different types of sequences in the gene assembly
	
```bash

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 2 -d ':' | sort | uniq
# Oryza_brachyantha.v1.4b

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 3 -d ':' | sort | uniq |wc -l
# 210

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 4 -d ':' | sort | uniq |wc -l
# 31979

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 5 -d ':' | sort | uniq |wc -l
# 31994

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 3 -d ' '| cut -f 6 -d ':' | sort | uniq
# -1
# 1
```

```bash
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 4 -d ' '| cut -f 1 -d ':' | sort | uniq
# gene
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 4 -d ' '| cut -f 2 -d ':' | sort | uniq | wc -l
# 32037
```

```bash
grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 5 -d ' '| cut -f 1 -d ':' | sort | uniq
# transcript

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 5 -d ' '| cut -f 2 -d ':' | sort | uniq | wc -l
# 32037

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 6 -d ' '| sort | uniq
# gene_biotype:protein_coding

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 7 -d ' '| sort | uniq
# transcript_biotype:protein_coding
```

#### Look at the different descriptions
```bash
grep "^>" Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | grep "description" | wc -l
# 5881

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 8-25 -d ' '| sort | uniq | head -3

# description:1,2-alpha-L-fucosidases [Source:Projected from Arabidopsis thaliana (AT4G34260) TAIR;Acc:AT4G34260]
# description:1,2-dihydroxy-3-keto-5-methylthiopentene dioxygenase [Source:UniProtKB/TrEMBL;Acc:J3LK94]

grep '^>' Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | cut -f 8-25 -d ' '| grep 'Projected from Arabidopsis thaliana' | sort | uniq | wc -l
# 3493
```

```bash
mkdir -p otherBlast; cd otherBlast

blastp -query ../filtered_pep.fasta -db ../filtered_proteins_db/filtered_proteins_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 50 -out brachyanthaBlast_1e5_max50.tsv
```

#### Remove *scaffold* from the fasta list
```bash
../filter_fasta.sh ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa ./no_scaffold.fa
```

#### Verification
```bash
grep -v "^>" ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa | wc -l
# 204282
grep -v "^>" no_scaffold.fa | wc -l
# 200541
```

#### Removal of sequences projected from A. thaliana
```bash
grep "^>" no_scaffold.fa | grep "Projected from Arabidopsis thaliana" | wc -l
# 3601
```

---
# 28th Dec. 2025

In `filter_fasta.sh`, replace `superscaffold` with `Source:Projected from Arabidopsis thaliana`

```bash
../filter_fasta.sh ./no_scaffold.fa ./o_branchyantha_only.fa

grep "^>" no_scaffold.fa | wc -l
# 31355
grep "^>" o_branchyantha_only.fa | wc -l
# 27754
```

```bash
grep "^>" no_scaffold.fa | grep "mitochond" | wc -l
# 36
grep "^>" no_scaffold.fa | grep "chloro" | wc -l
# 70
grep "^>" o_branchyantha_only.fa | grep "chloro" | wc -l
# 25
grep "^>" o_branchyantha_only.fa | grep "mitochond" | wc -l
# 17
```

### BLAST & Cluster
#### Prepare the BLAST
```bash
makeblastdb -in no_scaffold.fa -dbtype prot -out no_scaffold_db


# Building a new DB, current time: 12/28/2025 10:27:02
# New DB name:   /otherBlast/no_scaffold_db
# New DB title:  no_scaffold.fa
# Sequence type: Protein
# Keep MBits: T
# Maximum file size: 3000000000B
# Adding sequences from FASTA; added 31355 sequences in 2.73171 seconds.

makeblastdb -in o_branchyantha_only.fa -dbtype prot -out obrach_only_db


# Building a new DB, current time: 12/28/2025 10:27:37
# New DB name:   /otherBlast/obrach_only_db
# New DB title:  o_branchyantha_only.fa
# Sequence type: Protein
# Keep MBits: T
# Maximum file size: 3000000000B
# Adding sequences from FASTA; added 27754 sequences in 2.5455 seconds.
```

#### Run BLAST
```bash
blastp -query no_scaffold.fa -db no_scaffold_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 10 -out all_noScaffold_1e5_max10.tsv

blastp -query no_scaffold.fa -db no_scaffold_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 50 -out all_noScaffold_1e5_max50.tsv

blastp -query o_branchyantha_only.fa -db obrach_only_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 50 -out all_OBranch_1e5_max50.tsv

blastp -query o_branchyantha_only.fa -db obrach_only_db -outfmt 6 -evalue 1e-5 -num_threads 4 -max_target_seqs 10 -out all_OBranch_1e5_max10.tsv
```

#### Move the databases
```bash
mkdir -p obrach_only_db; mv -i obrach_only_db.* obrach_only_db
mkdir -p no_scaffold_db; mv -i no_scaffold_db.* no_scaffold_db
```

#### Prepare the cluster
```bash
../Scripts/BLAST_appendLenght.sh -i all_OBranch_1e5_max10.tsv -f ./o_branchyantha_only.fa -o append_OBranch_1e5_max10.tsv
../Scripts/BLAST_appendLenght.sh -i all_OBranch_1e5_max50.tsv -f ./o_branchyantha_only.fa -o append_OBranch_1e5_max50.tsv

../Scripts/BLAST_appendLenght.sh -i all_noScaffold_1e5_max10.tsv -f no_scaffold.fa -o append_noScaffold_1e5_max10.tsv
../Scripts/BLAST_appendLenght.sh -i all_noScaffold_1e5_max50.tsv -f no_scaffold.fa -o append_noScaffold_1e5_max50.tsv

../Scripts/BLAST_appendLenght.sh -i brachyanthaBlast_1e5_max50.tsv -f ../filtered_pep.fasta -o ./append_brachyanthaBlast_1e5_max50.tsv
```

Verification: Check the number of fields
```bash
awk '{NF}' append_noScaffold_1e5_max50.tsv | sort -n | uniq
awk 'NF != 14 {print NR, NF, $0}' append_noScaffold_1e5_max50.tsv
```

```bash
../Scripts/filter_cluster.sh -i ./append_OBranch_1e5_max10.tsv -s ./strict.tsv -m ./moderate.tsv -r ./relaxed.tsv
mkdir -p append_OBranch_1e5_max10; mv -i moderate* strict* relaxed* append_OBranch_1e5_max10

../Scripts/filter_cluster.sh -i ./append_OBranch_1e5_max50.tsv -s ./strict.tsv -m ./moderate.tsv -r ./relaxed.tsv ;
mkdir -p append_OBranch_1e5_max50; mv -i moderate* strict* relaxed* append_OBranch_1e5_max50

../Scripts/filter_cluster.sh -i ./append_noScaffold_1e5_max10.tsv -s ./strict.tsv -m ./moderate.tsv -r ./relaxed.tsv ;
mkdir -p append_noScaffold_1e5_max10; mv -i moderate* strict* relaxed* append_noScaffold_1e5_max10

../Scripts/filter_cluster.sh -i ./append_noScaffold_1e5_max50.tsv -s ./strict.tsv -m ./moderate.tsv -r ./relaxed.tsv;
mkdir -p append_noScaffold_1e5_max50; mv -i moderate* strict* relaxed* append_noScaffold_1e5_max50

../Scripts/filter_cluster.sh -i ./append_brachyanthaBlast_1e5_max50.tsv -s ./strict.tsv -m ./moderate.tsv -r ./relaxed.tsv ;
mkdir -p brachyanthaBlast_1e5_max50; mv -i moderate* strict* relaxed* brachyanthaBlast_1e5_max50
```

#### Upload the files on [Galaxy](https://usegalaxy.eu/root?tool_id=mcl)
Selected files:
||brachyanthaBlast_1e5_max50|append_noScaffold_1e5_max50|append_noScaffold_1e5_max10|append_OBranch_1e5_max50|append_OBranch_1e5_max10|
|----|----|----|----|----|----|
relaxed_clustering.tsv|relaxed_bBlast50.tsv|r_Scaff50.tsv|r_Scaff_10.tsv|r_Ob50.tsv|r_Ob10.tsv
moderate_clustering.tsv|moderate_bBlast50.tsv|m_Scaff50.tsv|m_Scaff_10.tsv|m_Ob50.tsv|<span style='color: red;'>**m_Ob10.tsv**</span>
strict_clustering.tsv|strict_bBlast50.tsv|s_Scaff50.tsv|s_Scaff_10.tsv|s_Ob50.tsv|<span style='color: red;'>**s_Ob10.tsv**</span>

#### Visualize the clusters
```bash
mkdir -p clusters/ ; mv -i /Downloads/*.tabular ./clusters/
```

### Ka/ Ks
#### Prepare the pairs
```bash
python ../prepare_kaks_input.py clusters/m_Ob10.tabular ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa ./ks_input
# Prepared 55201 gene pair directories with CDS and protein FASTA files in ./ks_input.

python ../prepare_kaks_input.py clusters/s_Ob10.tabular ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa ./ks_input_strict
# Prepared 22900 gene pair directories with CDS and protein FASTA files in ./ks_input_strict.
```

```bash
for d in ks_input/*; do for f in $d/*_prot.fasta; do a="${f%_prot.fasta}"; clustalw2 -quiet -align -infile=$f -outfile=$a.ali.aln; done; echo "$d - DONE"; done;

for d in ks_input_strict/*; do for f in $d/*_prot.fasta; do a="${f%_prot.fasta}"; clustalw2 -quiet -align -infile=$f -outfile=$a.ali.aln; done; echo "$d - DONE"; done;
```


## Updates
#### Filtering fasta file
New Usage:
```bash 
../filter_fasta.sh --input ../Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa --output ./o_branchyantha_only.fa --skip-superscaffold --skip-arabidopsis
```
Use `./filter_fasta.sh --help` for details.

#### Move Bash scripts into 1 repertory
```bash
mkdir -p ../Scripts; mv -i ../*.sh ../Scripts/
```
Replace `execute_yn00.sh_generalize.sh` with `execute_yn00_generalize.sh`

```bash
mkdir -p Data; mv -i Oryza_brachyantha.Oryza_brachyantha.v1.4b.* Data/
```

---
# 29th Dec 2025
### Ka/ Ks -- next
#### Prepare the pairs -- next
```bash
for n in ks_input/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ../pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy;   done;
done

for n in ks_input_strict/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ../pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy; done; done
```

#### Compute Ka, Ks, Omega
```bash
../Scripts/execute_yn00_generalize.sh ks_input/ ../yn00.ctl_master

../Scripts/execute_yn00_generalize.sh ks_input_strict/ ../yn00.ctl_master
```

#### Concatenate Ka, Ks, Omega in a tabular file
```bash
python ../collect_KS.py ks_input/ --out_file ./KaKs_all_ks_input.tsv
python ../collect_KS.py ks_input_strict/ --out_file ./KaKs_all_ks_input_strict.tsv
```
---
# 30th Dec 2025

## Visualization
### Clusters
#### Distribution
The distribution of the clusters was made/ shown on `display_cluster.R` .

#### GO Annotations
Among the biggest clusters, select a few genes and express them on PANTHER, DAVID.

### Ka/Ks
#### Distribution
`detect_Ks_peaks.R`

## -- To be Updated
```bash
cut -f 1,2 strict.tsv | sort -k1 > duplicates_strict.tsv
cut -f 1,2 treated_blast.tsv | sort -k1 > duplicates_raw.tsv
```

```bash
python detect_sigletons.py -raw duplicates_raw.tsv -reference duplicates_strict.tsv -out_duplicates duplicates_strict.txt -out_singletons singletons_strict.txt -delimiter " "
# Found 19543 duplicates and 72597 singletons.
# Duplicates written to: duplicates_strict.txt
# Singletons written to: singletons_strict.txt

python detect_sigletons.py -raw duplicates_raw.tsv -reference duplicates_strict.tsv -out_duplicates duplicates_strict_unique.txt -out_singletons singletons_strict_unique.txt
# Found 110176 duplicates and 74104 singletons.
# Duplicates written to: duplicates_strict_unique.txt
# Singletons written to: singletons_strict_unique.txt

mkdir -p singletons duplicates
mv -i duplicates_strict.txt duplicates_strict_unique.txt duplicates
mv -i singletons_strict.txt singletons_strict_unique.txt singletons
```
---
# Jan. 5th
```bash
for f in ../Data/*.fa ;  do echo $f; grep '^>' $f | cut -f 3 -d ' '| cut -
f 1 -d ':' | sort | uniq; done
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa
# chromosome
# superscaffold
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.dna.toplevel.fa
# chromosome
# superscaffold
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa
# chromosome
# superscaffold
```

```bash
for f in ../Data/*.fa ;  do echo $f; grep '^>' $f | grep -i "isoform"; echo "----------"; done
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa
# >OB03G16230.1 cds chromosome:Oryza_brachyantha.v1.4b:3:3783875:3788960:-1 gene:OB03G16230 gene_biotype:protein_coding transcript_biotype:protein_coding description:autoinhibited H(+)-ATPase isoform 10 [Source:Projected from Arabidopsis thaliana (AT1G17260) TAIR;Acc:AT1G17260]
# >OB03G41510.1 cds chromosome:Oryza_brachyantha.v1.4b:3:25245496:25246730:1 gene:OB03G41510 gene_biotype:protein_coding transcript_biotype:protein_coding description:cytochrome B5 isoform A [Source:Projected from Arabidopsis thaliana (AT1G26340) TAIR;Acc:AT1G26340]
# >OB04G16260.1 cds chromosome:Oryza_brachyantha.v1.4b:4:6650453:6658949:-1 gene:OB04G16260 gene_biotype:protein_coding transcript_biotype:protein_coding description:P4H isoform 1 [Source:Projected from Arabidopsis thaliana (AT2G43080) TAIR;Acc:AT2G43080]
# ----------
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.dna.toplevel.fa
# ----------
# ../Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa
# >OB03G16230.1 pep chromosome:Oryza_brachyantha.v1.4b:3:3783875:3788960:-1 gene:OB03G16230 transcript:OB03G16230.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:autoinhibited H(+)-ATPase isoform 10 [Source:Projected from Arabidopsis thaliana (AT1G17260) TAIR;Acc:AT1G17260]
# >OB03G41510.1 pep chromosome:Oryza_brachyantha.v1.4b:3:25245496:25246730:1 gene:OB03G41510 transcript:OB03G41510.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:cytochrome B5 isoform A [Source:Projected from Arabidopsis thaliana (AT1G26340) TAIR;Acc:AT1G26340]
# >OB04G16260.1 pep chromosome:Oryza_brachyantha.v1.4b:4:6650453:6658949:-1 gene:OB04G16260 transcript:OB04G16260.1 gene_biotype:protein_coding transcript_biotype:protein_coding description:P4H isoform 1 [Source:Projected from Arabidopsis thaliana (AT2G43080) TAIR;Acc:AT2G43080]
----------
```

```bash
for f in ../Data/*.fa ;  do echo $f; grep '^>' $f | grep "chromosome" | cu
t -f 3 | cut -f 3 -d ':'| uniq -c; echo "----------"; done
```
**Oryza_brachyantha.Oryza_brachyantha.v1.4b -- release-41**
| chromosome|cds|dna|pep|
| ----------|----------|----------|----------|
| 1|   4518| Not concerned | 4518
| 2|   3557| Not concerned |    3557
| 3|   3901| Not concerned |    3901
| 4|   2801| Not concerned |    2801
| 5|   2580| Not concerned |    2580
| 6|   2613| Not concerned |    2613
| 7|   2302| Not concerned |    2302
| 8|   2101| Not concerned |    2101
| 9|   1730| Not concerned |    1730
| 10|   1712| Not concerned |    1712 
| 11|   1828| Not concerned |    1828 
| 12|   1712| Not concerned |    1712 
---
# Jan. 6th

## Filter
```bash
./Scripts/filter_fasta.sh --input Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa --output pep_filtered.fa --skip-arabidopsis --skip-chlorop --skip-mitochond --skip-superscaffold

grep "^>" pep_filtered.fa | wc -l
27715
```

## BLASTP + Clustering
```bash
makeblastdb -in pep_filtered.fa -dbtype prot -out pep_db

# Building a new DB, current time: 01/06/2026 12:45:50
# New DB name:   /mnt/c/Users/Emilie-jeanne/Downloads/Project_Comparative_Genomics/Brachyantha/ComparativeGenomics/pep_db
# New DB title:  pep_filtered.fa
# Sequence type: Protein
# Keep MBits: T
# Maximum file size: 3000000000B
# Adding sequences from FASTA; added 27715 sequences in 2.46917 seconds

blastp -query pep_filtered.fa -db pep_db -outfmt 6 -evalue 1e-3 -num_threads 4 -out all_vs_all.tsv
```

```bash
mkdir -p pep_db; mv -i pep_db.* pep_db/
```


```bash
./Scripts/BLAST_appendLenght.sh -i all_vs_all.tsv -f ./pep_filtered.fa -o append_all_vs_all.tsv

awk '{print NF}' append_all_vs_all.tsv | sort -n | uniq
# 14

wc -l append_all_vs_all.tsv
# 659182 append_all_vs_all.tsv

./Scripts/filter_cluster.sh -i append_all_vs_all.tsv
mkdir -p clusters/ ; mv -i ./*.tabular ./*_pairs.tsv ./clusters/
mv -i relaxed.tsv strict.tsv moderate.tsv clusters/

# To check there are no multiple hits
cut -f 1,2 clusters/strict.tsv | sort -k 1 |uniq -c| sort -k1 -n | head -15
cut -f 1,2 clusters/strict.tsv | sort -k 1 |uniq -c| sort -k1 -n | tail -15
```

Run `blastDistribution.R`.
Run `blastCount.R`.
Run `display_cluster.R`.

```r
File: clusters/strict_mcl.tabular 
Number of clusters: 2937 
Min size: 2 
Max size: 28 
Mean size: 2.861764 
Median size: 2 

File: clusters/moderate_mcl.tabular 
Number of clusters: 3363 
Min size: 1 
Max size: 199 
Mean size: 4.25394 
Median size: 2 
```


## KaKs
```bash
python ./prepare_kaks_input.py clusters/strict_mcl.tabular Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa ./ks_pair_strict
# Prepared 13426 gene pair directories with CDS and protein FASTA files in ./ks_pair_strict.

for d in ks_pair_strict/*; do for f in $d/*_prot.fasta; do a="${f%_prot.fasta}"; clustalw2 -quiet -align -infile=$f -outfile=$a.ali.aln; done; echo "$d - DONE"; done;

for n in ks_pair_strict/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ./pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy; done; done

python ./prepare_kaks_input.py clusters/moderate_mcl.tabular Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.cds.all.fa Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.pep.all.fa ./ks_pair_moderate
# Prepared 115739 gene pair directories with CDS and protein FASTA files in ./ks_pair_moderate.

for d in ks_pair_moderate/*; do for f in $d/*_prot.fasta; do a="${f%_prot.fasta}"; clustalw2 -quiet -align -infile=$f -outfile=$a.ali.aln; done; echo "$d - DONE"; done;

for n in ks_pair_moderate/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ./pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy; done; done

# --------------------- LAST OF 6th JAN.
./Scripts/execute_yn00_generalize.sh ks_pair_strict/ ./yn00.ctl_master
python ./collect_KS.py ks_pair_strict/ --out_file ./KaKs_all_ks_input_strict.tsv

wc -l KaKs_all_ks_input_strict.tsv 
# 13427 KaKs_all_ks_input_strict.tsv

./Scripts/execute_yn00_generalize.sh ks_pair_moderate/ ./yn00.ctl_master
python ./collect_KS.py ks_pair_moderate/ --out_file ./KaKs_all_ks_input_moderate.tsv


conda install bioconda::mcscanx
conda install bioconda::orthofinder

python detect_sigletons.py -filtered clusters/moderate.tsv -reference append_all_vs_all.tsv
# Found 14306 duplicates and 5883 singletons.
# Duplicates written to: moderate.duplicate.tsv
# Singletons written to: moderate.singleton.tsv
python detect_sigletons.py -filtered clusters/strict.tsv -reference append_all_vs_all.tsv
# Found 8405 duplicates and 11784 singletons.
# Duplicates written to: strict.duplicate.tsv
# Singletons written to: strict.singleton.tsv
```

# Comparing genomes reveals features
## Conserved regions
## Distinct regions
## Evolution of genome across the time
### Gain/ Loss of function

# Jan. 8-10th

Worked on:
* **Jan 9th**
* detect_Ks_peaks.R
* **Jan 8th**
* Ks_SEdS.R
* display_cluster.R
* **Jan. 7th**
* substitutionLvls.R
* **Jan. 6th**
* blastCount.R
* blastDistribution.R

```bash
cd workdir/
python ./collect_KS.py ks_pair_strict/ --out_file ./KaKs_all_ks_input_strict.tsv
ls -l ks_pair_strict/ | wc -l

wc -l KaKs_all_ks_input_strict.tsv 
head -5 KaKs_all_ks_input_strict.tsv 
ls -l ks_pair_strict/ | tail -4

conda install bioconda::orthofinder
python detect_sigletons.py -filtered clusters/moderate.tsv -reference append_all_vs_all.tsv

sort moderate.duplicate.tsv|head -5

python detect_sigletons.py -filtered clusters/moderate.tsv -reference append_all_vs_all.tsv
python detect_sigletons.py -filtered clusters/strict.tsv -reference append_all_vs_all.tsv

for n in ks_pair_moderate/*; do for file in "$n"/*_cds.fasta; do a=$(basename "$file"); pair="${a%_cds.fasta}"; ./pal2nal.pl $n/$pair.ali.aln $file -output paml > $n/$pair.ali.phy; done; done

./Scripts/execute_yn00_generalize.sh ks_pair_moderate/ ./yn00.ctl_master
python ./collect_KS.py ks_pair_moderate/ --out_file ./KaKs_all_ks_input_moderate.tsv

mcscanx # Still dooesn't work

grep "^>" Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.41.chr.gff3 | grep "isoform" | wc -l
grep "^>" Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.41.chr.gff3 | grep "iso" | wc -l
grep "isoform" Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.41.chr.gff3 | wc -l
grep "isoform" Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.41.chr.gff3


cut -f 1,2 append_all_vs_all.tsv | sort -d -k 1 | uniq -c | sort -k 1 -n | head -5
      # 1 OB01G10010.1    OB04G17220.1
      # 1 OB01G10010.1    OB06G31350.1
      # 1 OB01G10010.1    OB08G22190.1
      # 1 OB01G10010.1    OB08G29380.1
      # 1 OB01G10010.1    OB09G23820.1
cut -f 1,2 append_all_vs_all.tsv | sort -d -k 1 | uniq -c | sort -k 1 -n | tail -5
      # 1 OB12G27120.1    OB05G17990.1
      # 1 OB12G27120.1    OB05G18310.1
      # 1 OB12G27120.1    OB07G18120.1
      # 1 OB12G27120.1    OB10G11220.1
      # 1 OB12G27120.1    OB12G15400.1

ls Data/Oryza_sativa.IRGSP-1.0.*
# Data/Oryza_sativa.IRGSP-1.0.41.chr.gff3.gz
# Data/Oryza_sativa.IRGSP-1.0.cds.all.fa.gz
# Data/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz

for file in Data/*; do gunzip $file; done

ls Data/Oryza_sativa.IRGSP-1.0.*
# Data/Oryza_sativa.IRGSP-1.0.41.chr.gff3  Data/Oryza_sativa.IRGSP-1.0.pep.all.fa
# Data/Oryza_sativa.IRGSP-1.0.cds.all.fa
for file in Data/Oryza_sativa.IRGSP-1.0.*; do head -5 $file; done
# ##gff-version 3
# ##sequence-region   1 1 43270923
# ##sequence-region   10 1 23207287
# ##sequence-region   11 1 29021106
# ##sequence-region   12 1 27531856
# >BAC19866 cds chromosome:IRGSP-1.0:Mt:84886:90873:-1 gene:BAC19866 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:nad7 description:NADH dehydrogenase subunit 7
# ATGACGACTAGGAACGGGCAAATCAAGAATTTCACTTCGAATTCCGGACCTCAACATCCT
# GCTGCTCATGGTGTTTCACGATCAGTATTGGAAATGAACGGAGAAGTGGTGGAACGTGCG
# GAACCACATATTGGATCACTCCAACTAGAGGGACTGAGAAATTAATCGAGTACAAAACTT
# ATCTTCAAGCTTTACCTTATTTTGATCGTTCAGACTATGTTTCTACGATGGCCCAAGAAC
# >BAC19866 pep chromosome:IRGSP-1.0:Mt:84886:90873:-1 gene:BAC19866 transcript:BAC19866 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:nad7 description:NADH dehydrogenase subunit 7
# MTTRNGQIKNFTSNSGPQHPAAHGVSRSVLEMNGEVVERAEPHIGSLQLEGLRN*SSTKL
# IFKLYLILIVQTMFLRWPKNTLILQP*RDF*IVRYHYELNIYECYSVK*LEFQIIHLLQL
# LMLWMWEHQLRSFGLLRSGRNCWNSMKESREPGCMPVSYDLVEWHKICLLAYVEILIPPH
# NNLLLVSTN*KRCQPATVSGNND*WILVLSLHSKQRIGDSVV*C*EVLGYAGIREEQHLT
for file in Data/Oryza_sativa.IRGSP-1.0.41.chr.gff3; do grep -v "#"  $file|head -5; done
# 1       RAP-DB  chromosome      1       43270923        .       .       .    ID=chromosome:1;Alias=Chr1
# 1       ensembl gene    2983    10815   .       +       .       ID=gene:Os01g0100100;biotype=protein_coding;description=RabGAP/TBC domain containing protein. (Os01t0100100-01);gene_id=Os01g0100100;logic_name=irgspv1.0-20170804-genes
# 1       ensembl mRNA    2983    10815   .       +       .       ID=transcript:Os01t0100100-01;Parent=gene:Os01g0100100;biotype=protein_coding;transcript_id=Os01t0100100-01
# 1       ensembl exon    2983    3268    .       +       .       Parent=transcript:Os01t0100100-01;Name=Os01t0100100-01.exon1;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Os01t0100100-01.exon1;rank=1;version=1
# 1       ensembl five_prime_UTR  2983    3268    .       +       .       Parent=transcript:Os01t0100100-01
for file in Data/Oryza_sativa.IRGSP-1.0.41.chr.gff3; do grep -v "^#"  $file|head -5; done
# 1       RAP-DB  chromosome      1       43270923        .       .       .    ID=chromosome:1;Alias=Chr1
# 1       ensembl gene    2983    10815   .       +       .       ID=gene:Os01g0100100;biotype=protein_coding;description=RabGAP/TBC domain containing protein. (Os01t0100100-01);gene_id=Os01g0100100;logic_name=irgspv1.0-20170804-genes
# 1       ensembl mRNA    2983    10815   .       +       .       ID=transcript:Os01t0100100-01;Parent=gene:Os01g0100100;biotype=protein_coding;transcript_id=Os01t0100100-01
# 1       ensembl exon    2983    3268    .       +       .       Parent=transcript:Os01t0100100-01;Name=Os01t0100100-01.exon1;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Os01t0100100-01.exon1;rank=1;version=1
# 1       ensembl five_prime_UTR  2983    3268    .       +       .       Parent=transcript:Os01t0100100-01
for file in Data/Oryza_sativa.*.fa ; do grep "^>" $file|head -5; echo ""; done
# >BAC19866 cds chromosome:IRGSP-1.0:Mt:84886:90873:-1 gene:BAC19866 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:nad7 description:NADH dehydrogenase subunit 7
# >BAC19897 cds chromosome:IRGSP-1.0:Mt:338572:340147:1 gene:BAC19897 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:cox1 description:Cytochrome c oxidase subunit1
# >BAC19886 cds chromosome:IRGSP-1.0:Mt:294692:297681:1 gene:BAC19886 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS description:Ribosomal protein L2
# >BAC19869 cds chromosome:IRGSP-1.0:Mt:170836:174298:1 gene:BAC19869 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:RPS3 description:Ribosomal protein S3
# >Os08t0254300-00 cds chromosome:IRGSP-1.0:8:9399521:9402233:1 gene:Os08g0254300 gene_biotype:protein_coding transcript_biotype:protein_coding description:Similar to OSIGBa0112M24.2 protein. (Os08t0254300-00)
# 
# >BAC19866 pep chromosome:IRGSP-1.0:Mt:84886:90873:-1 gene:BAC19866 transcript:BAC19866 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:nad7 description:NADH dehydrogenase subunit 7
# >BAC19897 pep chromosome:IRGSP-1.0:Mt:338572:340147:1 gene:BAC19897 transcript:BAC19897 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:cox1 description:Cytochrome c oxidase subunit1
# >BAC19886 pep chromosome:IRGSP-1.0:Mt:294692:297681:1 gene:BAC19886 transcript:BAC19886 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS description:Ribosomal protein L2
# >BAC19869 pep chromosome:IRGSP-1.0:Mt:170836:174298:1 gene:BAC19869 transcript:BAC19869 gene_biotype:nontranslating_CDS transcript_biotype:nontranslating_CDS gene_symbol:RPS3 description:Ribosomal protein S3
# >Os08t0254300-00 pep chromosome:IRGSP-1.0:8:9399521:9402233:1 gene:Os08g0254300 transcript:Os08t0254300-00 gene_biotype:protein_coding transcript_biotype:protein_coding description:Similar to OSIGBa0112M24.2 protein. (Os08t0254300-00)
# blastp -query Data/Oryza_sativa.IRGSP-1.0.pep.all.fa -db pep_db -outfmt 6 -evalue 1e-3 -num_threads 4 -out interspecie.blast
# # BLAST Database error: No alias or index file found for protein database [pep_db] in search path [/home/epic_robinson/workdir/compGEN::]. Please verify the spelling of the BLAST database and its molecule type.
# blastp -query Data/Oryza_sativa.IRGSP-1.0.pep.all.fa -db pep_db/pep_db -outfmt 6 -evalue 1e-3 -num_threads 4 -out interspecie.blast

makeblastdb -in Data/Oryza_sativa.IRGSP-1.0.pep.all.fa -dbtype prot -out os_db
# Building a new DB, current time: 01/09/2026 18:53:50
# New DB name:   /home/epic_robinson/workdir/compGEN/os_db
# New DB title:  Data/Oryza_sativa.IRGSP-1.0.pep.all.fa
# Sequence type: Protein
# Keep MBits: T
# Maximum file size: 3000000000B
# Adding sequences from FASTA; added 42377 sequences in 0.839024 seconds.
blastp -query pep_filtered.fa -db os_db -outfmt 6 -evalue 1e-3 -num_threads 4 -out ob_vs_os_db.blast
mkdir -p os_db; mv -i os_db.* os_db

# ------------------------------------------------------

python makeBedFile.py Data/Oryza_brachyantha.Oryza_brachyantha.v1.4b.41.chr.gff3 pep_filtered.fa

python PyScripts/getTAG.py -o strict_TAGs.tsv
python PyScripts/getTAG.py -b clusters/moderate.tsv -m clusters/moderate_mcl.tabular -o moderate_TAGs.tsv

python PyScripts/TAGAnalysis.py
Maximum array size: 16

TAG arrays by size bins:
  2-2: 624
  3-5: 213
  6-9: 20
  10<=: 4

TAG pairs by spacer range:
  0-0: 1300
  1-5: 885
  6-10: 115

TAG arrays by max spacer range:
  0-0: 861
  1-5: 0
  6-10: 0

(compGen) epic_robinson@698b0d05b5e5:~/workdir/compGEN$ python PyScripts/getTAG.py -b clusters/moderate.tsv -m clusters/moderate_mcl.tabular -o moderate_TAGs.tsv
TAG arrays written to moderate_TAGs.tsv
(compGen) epic_robinson@698b0d05b5e5:~/workdir/compGEN$ python PyScripts/TAGAnalysis.py 
Maximum array size: 21

TAG arrays by size bins:
  2-2: 854
  3-5: 364
  6-9: 43
  10<=: 12

TAG pairs by spacer range:
  0-0: 2163
  1-5: 1997
  6-10: 352

TAG arrays by max spacer range:
  0-0: 1273
  1-5: 0
  6-10: 0
```

# Jan. 11th
```bash
# Final - right way
python PyScripts/TAGAnalysis.py -b clusters/strict.tsv -m clusters/strit_mcl.tabular -g gene_pos.bed
```
Spacer_range|TAG_pairs|   TAG_arrays|  Max_array|   Arrays_size2|Arrays_3_5|  Arrays_6_9|  Arrays_10plus
|----|----|----|----|----|----|----|----|
(0, 0) | 1635|898| 9|   706| 183| 9|   0
(1, 5) | 3297|733| 17|  445| 242| 36|  10
(6, 10) |1725|181| 21|  99|  55|  18|  9


```bash
python PyScripts/TAGAnalysis.py -b clusters/strict.tsv -m clusters/strict_mcl.tabular -g gene_pos.bed
```
Spacer_range|TAG_pairs|   TAG_arrays|  Max_array|   Arrays_size2|Arrays_3_5|  Arrays_6_9|  Arrays_10plus
|----|----|----|----|----|----|----|----|
(0, 0)|903| 557| 9|   459| 94|  4|   0
(1, 5)|1657|479| 16|  318| 142| 15|  4
(6, 10) |623| 96|  16|  58|  30|  5|   3


Keep only the top hits
```bash
# Sort by query (col1) and descending bit score (col12)
sort -k1,1 -k12,12nr OB_vs_Os.blastp > OB_vs_Os.sorted.blastp
# Keep the first occurrence for each query (the best hit)
awk '!seen[$1]++' OB_vs_Os.sorted.blastp > OB_to_Os_besthit.txt

wc -l OB_to_Os_besthit.txt
# 21857 OB_to_Os_besthit.txt
```

Find  *O. Sativa* homologs genes to *O. Brachyantha*
```bash
cut -f 1,2 OB_to_Os_besthit.txt > OB_to_Os_IDs.txt
# OB01G10010.1    Os01t0100100-01
# OB01G10030.1    Os01t0100500-01
# OB01G10040.1    Os01t0100600-01

awk '{print $1, NF}' clusters/strict_mcl.tabular > clusters/strict_mcl_IDs.txt
awk '{print $1, NF}' clusters/moderate_mcl.tabular > clusters/moderate_mcl_IDs.txt

sort -k 1 clusters/moderate_mcl_IDs.txt > clusters/moderate_mcl_IDs.sorted.txt; rm clusters/moderate_mcl_IDs.txt
sort -k 1 clusters/strict_mcl_IDs.txt > clusters/strict_mcl_IDs.sorted.txt; rm clusters/strict_mcl_IDs.txt

awk '{if(NR>1){print $1 "\t" $2}}' ../Oryza_sativa/data/protein_info.txt > Data/Os_Prot2Gen.txt
# Os12t0640700-01 Os12g0640700
# Os12t0640800-01 Os12g0640800

join -1 1 -2 1 OB_to_Os_IDs.txt clusters/strict_mcl_IDs.sorted.txt > strict_OB_OS_joint.tsv
join -1 1 -2 1 OB_to_Os_IDs.txt clusters/moderate_mcl_IDs.sorted.txt > moderate_OB_OS_joint.tsv
# OB01G10050.1 Os01t0100800-01 3
# OB01G10200.1 Os01t0102700-01 2

sort -k 2 strict_OB_OS_joint.tsv | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > strict_OB2GO_for_PANTHER.txt
# Os01t0101150-00 OB01G10070.1 2 Os01g0101150
# Os01t0103000-01 OB01G10230.1 2 Os01g0103000
sort -k 2 moderate_OB_OS_joint.tsv | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > moderate_OB2GO_for_PANTHER.txt

awk '{print $4 "\t" $3}' strict_OB2GO_for_PANTHER.txt > strict_OB2GO_for_PANTHER.clean.txt
awk '{print $4 "\t" $3}' moderate_OB2GO_for_PANTHER.txt > moderate_OB2GO_for_PANTHER.clean.txt
# Os01g0100800    3
# Os01g0102700    2

sort -k 1 moderate_OB2GO_for_PANTHER.clean.txt | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}'| sort -k 2 -nr > moderate_OB2GO_for_PANTHER.clean.uniq.txt
# Os05g0319700    199
# Os07g0134200    126
# Os04g0101400    118
sort -k 1 strict_OB2GO_for_PANTHER.clean.txt | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}'| sort -k 2 -nr > strict_OB2GO_for_PANTHER.clean.uniq.txt

awk '{if($2 >= 10){print $1}}' strict_OB2GO_for_PANTHER.clean.uniq.txt > strict_OB2GO_for_PANTHER.Overrep.txt
awk '{if($2 >= 10){print $1}}' moderate_OB2GO_for_PANTHER.clean.uniq.txt > moderate_OB2GO_for_PANTHER.Overrep.txt
# Os12g0640700
# Os12g0636800
# Os12g0484600
```

## Analyze Overrepresentation on TAGs

```bash
python PyScripts/TAGAnalysis.py -b clusters/strict.tsv -m clusters/strict_mcl.tabular -g gene_pos.bed
Spacer_range    TAG_pairs       TAG_arrays      Max_array       Arrays_size2    Arrrays_6_9       Arrays_10plus
# TAGs-0 saved to: TAGS/TAGs_0_spacers.txt
# (0, 0)  903     557     9       459     94      4       0
# TAGs-0 saved to: TAGS/TAGs_0_spacers.txt
# (1, 5)  1657    479     16      318     142     15      4
# TAGs-0 saved to: TAGS/TAGs_0_spacers.txt
# (6, 10) 623     96      16      58      30      5       3
# Non-TAGs saved to: TAGS/non_TAGs.txt

awk '{print $1, NF}' TAGS/TAGs_0_spacers.txt | sort > TAGS/TAGs_0_spacers_IDs.txt;
join -1 1 -2 1 OB_to_Os_IDs.txt TAGS/TAGs_0_spacers_IDs.txt > TAGS/temp0.txt;
sort -k 2 TAGS/temp0.txt | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > TAGS/OB2GO_temp.txt;
rm TAGS/temp0.txt;
awk '{print $4 "\t" $3}' TAGS/OB2GO_temp.txt > TAGS/OB2GO_temp.clean.txt;
rm TAGS/OB2GO_temp.txt;
sort -k 1 TAGS/OB2GO_temp.clean.txt | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}'| sort -k 2 -nr > TAGS/OB2GO_temp.uniq.txt;
rm TAGS/OB2GO_temp.clean.txt;
awk '{print $1}' TAGS/OB2GO_temp.uniq.txt > TAGS/TAGs_OB2GO_Overrep.txt;
rm TAGS/OB2GO_temp.uniq.txt;
```

GO analyses on [PANTHER](https://pantherdb.org/) : Select "Statistical overrepresentation test"
Visualization with `EnrichmentPlot.R`

# Jan 12th

## Rerun for statistics & plots
```bash
python PyScripts/mcl_cluster_size_plot.py -i clusters/strict_mcl.tabular
# Cluster and gene counts per category:
#   size_2: 1906 clusters, 3812 unique genes
#   size_3_9: 987 clusters, 3976 unique genes
#   size_10_19: 39 clusters, 500 unique genes
#   size_20_plus: 5 clusters, 117 unique genes
# Plot saved to: plots/cluster_distribution_strict_mcl.png

python PyScripts/mcl_cluster_size_plot.py -i clusters/moderate_mcl.tabular
# Cluster and gene counts per category:
#   size_2: 1722 clusters, 3444 unique genes
#   size_3_9: 1428 clusters, 6141 unique genes
#   size_10_19: 143 clusters, 1894 unique genes
#   size_20_plus: 69 clusters, 2826 unique genes
# Plot saved to: plots/cluster_distribution_moderate_mcl.png
```

```bash
# === For the non-TAGs & singletons
## non-TAGs
awk '{print $1, NF}' TAGS/non_TAGs.txt | sort > TAGS/non_TAGs_0_spacers_IDs.txt
join -1 1 -2 1 OB_to_Os_IDs.txt TAGS/TAGs_0_spacers_IDs.txt > TAGS/temp0.txt
sort  -k 2 TAGS/temp0.txt | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > TAGS/OB2GO_temp.txt
rm TAGS/temp0.txt
awk '{print $4 "\t" $3}' TAGS/OB2GO_temp.txt > TAGS/OB2GO_temp.clean.txt
rm TAGS/OB2GO_temp.txt
sort -k 1 TAGS/OB2GO_temp.clean.txt | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}'| sort -k 2 -nr > TAGS/OB2GO_temp.uniq.txt
rm TAGS/OB2GO_temp.uniq.txt
awk '{print $1}' TAGS/OB2GO_temp.uniq.txt > TAGS/non-TAGs_OB2GO_Overrep.txt

## singletons
awk '{print $1, NF}' strict.singleton.tsv | sort > TAGS/strict.singleton_IDs.txt;
join -1 1 -2 1 OB_to_Os_IDs.txt TAGS/strict.singleton_IDs.txt > TAGS/temp0.txt;
sort -k 2 TAGS/temp0.txt | join -1 2 -2 1 - Data/Os_Prot2Gen.txt > TAGS/OB2GO_temp.txt;
rm TAGS/temp0.txt;
awk '{print $4 "\t" $3}' TAGS/OB2GO_temp.txt > TAGS/OB2GO_temp.clean.txt;
rm TAGS/OB2GO_temp.txt;
sort -k 1 TAGS/OB2GO_temp.clean.txt | uniq -c | sort -nr -k 1 | awk '{print $2 "\t" $3}'| sort -k 2 -nr > TAGS/OB2GO_temp.uniq.txt;
rm TAGS/OB2GO_temp.clean.txt;
awk '{print $1}' TAGS/OB2GO_temp.uniq.txt > TAGS/strict.singleton_OB2GO_Overrep.txt;
rm TAGS/OB2GO_temp.uniq.txt;
```