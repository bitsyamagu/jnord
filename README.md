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

## Installation
+ Copy jar files from target directory and dependent libraries(gral-core, commons-math3, htsjdk) to your environment.
+ Fix paths in the shortcut script "jnord" 
+ Put "jnord" script where you can invoke through PATH environment variable
  
## References
+ Accurate and exact CNV identification from targeted high-throughput sequence data
    https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-12-184
