# Euterpe precatoria Mitochondrial Genome Assembly

This repository provides the code and methodology used to assemble the mitochondrial genome of the Amazonian palm Euterpe precatoria. Due to the high structural complexity and recombination activity typical of the Arecaceae family, this assembly involved a hybrid approach of automated toolsets and manual graph resolution.

## Bioinformatic Implementation

### 1.Read Filtering and Preparation

We used minimap2 to pull mitochondrial-related reads and samtools to deduplicate them.


Mapping against reference
```
minimap2 -ax map-hifi [reference].fasta all_reads_palm.fastq.gz > [reference]_all_reads.sam
```
Deduplication and cleaning
```
samtools view -b -F 4 [reference].sam | samtools sort -n -o sorted_name.bam
samtools fixmate -m sorted_name.bam fixmate.bam
samtools sort fixmate.bam -o sorted_pos.bam
samtools markdup -r sorted_pos.bam cleaned.bam
samtools fasta cleaned.bam > final_input_reads.fa
```

### 2. De Novo Assembly
The assembly was performed using Flye specifically tuned for HiFi data.
Bash
```
flye --pacbio-hifi final_input_reads.fa --o flye_assembly
```
### 3. Manual Graph Curation
The assembly graph was manually resolved in Bandage to navigate the recombination topology.Master Content Sequence: 1,268,308 bpExcluded Segments: edge_4 (identified as a NUMT).

### 4. Annotation
Annotation was executed via the PMGA pipeline, incorporating tRNAscan-SE, ARAGORN, and Deepred-MT.

📊 Summary StatisticsFeatureCountUnique Genes76Protein-coding41tRNAs32rRNAs3Coverage746X
