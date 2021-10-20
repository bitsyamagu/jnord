# jnord
Java version of Nord's CNV analysis (Target CNV detection)

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
   - Path to BED formatted file of your capture kit
+ --refgene (required)
   - Path to refGene.txt file obtained from UCSC web site
   - This data is used for plotting graph
+ --samples (required)
   - Put target sample names(written in SM fieled of the @RG record in bam) to output graph files.

## Tips
If you feel jnord requires too many mandatory commmand line arguments, 
customize 16th line of "jnord" shortcut script like this:
   ```exec "$JAVA -cp $CLASSPATH jnord.Main2 --refGene /path/to/default/refGene.txt --sureselect /path/to/default/SureSelect.bed " . join(" ", @ARGV);```
   
When you put some command line options twice or more, only the last one is used for analysis.

## Recommendations
+ 30 samples or more 
   - At least 30 samples(i.e. 30 bam files) or more samples are required for this analysis
   - All of samples should be sequenced under the same or similar conditions(same flowcell, same sequencer, same date, same operator etc) as possible as
+ Check gene names
   - Are your target gene names included in both of refGene.txt and capture bed file? If these three files don't have consistent record(genes' name and regions), then jnord will fail.

## References
+ Accurate and exact CNV identification from targeted high-throughput sequence data
    https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-184
