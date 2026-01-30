

# Euterpe precatoria Mitochondrial Genome Assembly

This repository provides the code and methodology used to assemble the mitochondrial genome of the Amazonian palm Euterpe precatoria. Due to the high structural complexity and recombination activity typical of the Arecaceae family, this assembly involved a hybrid approach of automated toolsets and manual graph resolution.
The assembly represents a complex, recombination-driven architecture spanning approximately 1.27 Mb.

## Methods: Assembly & Curation

### 1. Reference-Guided Read Filtering
To isolate mitochondrial signal from the total genomic PacBio HiFi data, we mapped raw reads to known Arecaceae mitogenomes and a preliminary OATK contig.

* **References:**
    * *Cocos nucifera* ([NC_031696](https://www.ncbi.nlm.nih.gov/nuccore/NC_031696))
    * *Phoenix dactylifera* ([NC_016740](https://www.ncbi.nlm.nih.gov/nuccore/NC_016740))

**Example Mapping Command:**

```bash
minimap2 -ax map-hifi cocos_mito_NC_031696.fasta ../all_reads_palm.fastq.gz > cocos_all_reads.sam
```

### 2. Read Processing & Deduplication

We used SAMtools to filter for mapped reads and remove potential duplicates to ensure a high-quality input for the assembler.

```Bash
# Convert to BAM and filter for mapped reads (-F 4)```
samtools view -b -F 4 cocos_all_reads.sam > mapped_reads.bam 

# Prepare for duplicate marking (Collate > Fixmate > Sort)
samtools collate -o collated.bam mapped_reads.bam
samtools fixmate -m collated.bam fixmate.bam
samtools sort -o positionsort.bam fixmate.bam

# Remove duplicates and export to FASTA
samtools markdup -r positionsort.bam markdup.bam
samtools fasta markdup.bam > input_reads_flye.fa
```

### 3. De Novo Assembly & Graph Resolution

The cleaned reads were assembled using Flye. Despite high coverage, the graph remained unresolved due to high-frequency recombination.

```Bash
flye --pacbio-hifi input_reads_flye.fa --o flye_output
```

Manual Path Curation:The assembly graph was visualized in Bandage. We identified a Nuclear Mitochondrial DNA (NUMT) segment (edge_4) and excluded it. The final master sequence was resolved by manually tracing the following recombination path:
```
edge_5+, edge_12+, edge_1+, edge_14+, edge_16+, edge_7+, edge_10+, edge_18+, edge_19+, edge_18-, edge_10-, edge_7-, edge_16-, edge_17+, edge_11+, edge_17-, edge_16+, edge_7+, edge_10+, edge_18+, edge_19+, edge_18-, edge_10-, edge_7-, edge_16-, edge_17+, edge_11+, edge_17-, edge_16+, edge_7+, edge_10+, edge_18+, edge_19+, edge_18-, edge_10-, edge_7-, edge_16-, edge_17+, edge_11+, edge_17-, edge_16+, edge_7+, edge_10+, edge_13-, edge_1+, edge_2+, edge_6-, edge_15-, edge_3+
```
### 4. Technical Validation & Quality Control

To validate the final 1,268,308 bp assembly, we mapped the original reads back to the curated master sequence to calculate depth and consistency.

```Bash
# Map reads to the final assembly
minimap2 -ax map-hifi E_precatoria_MT_final.fasta all_reads_palm.fastq.gz > validation.sam

# Sort and index
samtools view -b validation.sam | samtools sort -o validation_sorted.bam
samtools index validation_sorted.bam

# Calculate mean coverage depth
samtools depth validation_sorted.bam > mt_depth.txt
awk '{sum+=$3; count++} END {print "Mean Depth: ",sum/count}' mt_depth.txt
```

## Genome Annotation


### 1. Method

The resolved sequence was annotated using the PMGA (Plant Mitochondrial Genome Annotation) pipeline.

**Homology Search**: Homology-based identification against 319 plant mitogenomes.

**tRNA/rRNA Prediction**: Identified using tRNAscan-SE 2.0 and ARAGORN.RNA 


### 2. Results

#### Genome Features

The assembly reveals a significantly expanded mitochondrial architecture compared to other palms.


```
Metric,Value
Total Length,"1,268,308 bp"
Mean GC Content,42.9%
Mean Read Coverage,746X
Mapping Rate,1.20%
Total Gene Features,177
Unique Genes,76
```


#### Summary

The 76 unique genes identified are distributed as follows:
```
Protein-coding: 41 unique genes (56 total features).
tRNA Genes: 32 unique genes (118 features).
rRNA Genes: 3 unique genes (rrn5, rrn18, rrn26).
```
This profile reflects the characteristic gene duplication and fragmentation documented in the Arecaceae family.


## Software Requirements

Minimap2 (v2.28)  
Flye (v2.9.5)  
Samtools (v1.21)  
Bandage (v0.8.1)  
PMGA Pipeline  

