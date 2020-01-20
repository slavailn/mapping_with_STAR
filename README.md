## Mapping reads with STAR 

STAR requires two steps:
1. Index creation.
2. Read mapping

### 1. Index creation
To create an index of a reference genome, we need the fasta file of a genome and GTF/GFF file with genomic feature coordinates.
Indexing options:
   * --runThreadN: number of threads
   * --runMode *genomeGenerate*
   * --genomeDir */path/to/index_dir* : directory where index files will be stored. The directory needs to be created beforehand
   * --genomeFastaFiles */path/to/fasta* : path to one or more genome fasta files
   * --sjbGTFfile */path/to/gtf* : path to GTF file with annotations
   * --sjdbOverhang *read length - 1* : set this maximum read length of the sequencing reads minus one base. Defaults to 99. Rebuild index if the read length in the project significantly deviates from the one used to build the index.
   * --sjdbGTFtagExonParent : use this option with GFF3 files!  

Create genome index:
    
    ```
    STAR --runThreadN num_threads --runMode genomeGenerate --genomeFastaFiles path/to/genome.fasta \ 
    --sjdbGTFfile path/to/annotation.gtf --sjdbOverhang max(readLength)-1 --genomeDir /path/to/index_directory
    
    ```

### 2. Read mapping (basic options)
    *   --runThreadN: number of threads.
    *   --genomeDir: /path/to/index_dir
    *   --readFilesIn: /path/to/read1 /path/to/read2
    *   --outSAMtype BAM SortedByCoordinate (many types of output are available, but this one is the most useful)
    *   --outFileNamePrefix: *prefix* (use the sample name as a prefix)
       
       Example of mapping command:
     
       ```
       STAR --runThreadN 16 --genomeDir /path/to/star_index_dir --readFilesIn sample_R1.fq sample_R2.fq \
            --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sample
       ```      

NOTE: STAR tends to open may temporary files and this problem is aggravated by number of threads used. If the program crashes due to high number of open files,
      you can either increase *ulimit* or decrease a number of threads.
NOTE: Another problem that occurs frequently is STAR crashing due to broken fastq files. For example paired-end files may contain singletons or any fastq files may 
      containe the reads that are too short.
      The problem with sigletons can be solved using *repair* tool from BBtools like so:
     
      ```
      repair.sh in1=s_1.fastq in2=s_2.fastq out1=s1_fixed_1.fq out2=s2_fixed_2.fq repair

      ```
     
      Short reads can be filtered out using *reformat* tool from BBtools:
     
      ```
      reformat.sh in=file in2=file2 out=outfile out2=outfile2 minlen=25
     
      ```

       
       
