# jnord
Java version of Nord (Target CNV detection)

## Citing jnord
- Efficient detection of copy-number variations using exome data: Batch- and sex-based analyses. Hum Mutat
. 2021 Jan;42(1):50-65. doi: 10.1002/humu.24129. Epub 2020 Nov 11.
  - https://onlinelibrary.wiley.com/doi/epdf/10.1002/humu.24129
## Requirements
- Maven:
  - For compilation
  -   https://maven.apache.org/
- JDK8 or later
  - For compilation and execution 
- Bait file
  - If you use SureSelect for DNA capture, use the padded bed file of the version of your capture kit.
- refGene.txt
  - For plotting graph

## Build
+ Compile
  - ```mvn compile```
+ Make a JAR binary file
  - ```mvn package```
  - This will produce a binary: target/jnord_project-0.0.1-SNAPSHOT.jar

## Installation
+ Copy target/jnord_project-0.0.1-SNAPSHOT.jar and dependent libraries(gral-core, commons-math3, htsjdk) to your environment.
+ Fix paths in the shortcut script "jnord" 
+ Put "jnord" script where you can invoke through PATH environment variable

# Usage
- Command usage:
   - ```jnord --threads [threads] --genes [genes' list] --bamdir [bam dir] --sureselect [sureselect bed] --refgene [refGene.txt] --samples [Sample ID1] [Sample_ID2]...```
+ --threads (optional)
   - default: 4
+ --genes (required)
   - A file contains list of target genes. 
   - File format: Gene names should be listed line by line.
      - These gene names should be included in files capture bed file and refGene.txt. 
+ --bamdir (required)
   - A directory contains BAM and BAI files
+ --sureselect (required)
   - BED formatted file of your capture kit
+ --refgene
   - refGene.txt file obtained from UCSC web site
   - This data is used for plotting graph
+ --samples
   - Put target sample names(written in SM fieled of the @RG record in bam) to output graph files.

## References
+ Accurate and exact CNV identification from targeted high-throughput sequence data
    https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-184
