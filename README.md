## Mapping reads with STAR 

STAR requires two steps:
1. Index creation.
2. Read mapping

### 1. Index creation
To create an index of a reference genome, we need the fasta file of a genome and GTF/GFF file with genomic feature coordinates.
Indexing options:
    --runThreadN: number of threads
    --runMode *genomeGenerate*
    --genomeDir */path/to/index_dir* : directory where index files will be stored. The directory needs to be created beforehand
    --genomeFastaFiles */path/to/fasta* : path to one or more genome fasta files
    --sjbGTFfile */path/to/gtf* : path to GTF file with annotations
    --sjdbOverhang *read length - 1* : set this maximum read length of the sequencing reads minus one base. Defaults to 99. Rebuild index if the read length in the project significantly deviates from the one used to build the index.
    --sjdbGTFtagExonParent : use this option with GFF3 files!  

Create genome index:
    ```
    STAR --runThreadN num_threads --runMode genomeGenerate --genomeFastaFiles path/to/genome.fasta \ 
    --sjdbGTFfile path/to/annotation.gtf --sjdbOverhang max(readLength)-1 --genomeDir /path/to/index_directory
    
    ```
