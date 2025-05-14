# serratia-code
The script enables systematic identification of Type VI Secretion System (T6SS) loci and their associated orthologous groups in *Serratia* genomes.

## Usage
The serratia_locus_finder.py can detect T6SS loci from *Proteinortho v6.3.1* output csv files. You can run it through python3 when giving it a csv file and it will directly output results on the screen.

```shell
python3 serratia_locus_finder.py <filename.tsv>
```
After collecting all the T6SS distribution in *Serratia*, you can canculate the P/A-Value of all the orthologous groups. First you use *CD-HIT v4.6* to cluster all the *Serratia* proteins, and then you can get the output as a txt file. Name this file `input.txt` as the cdhit-to-matrix.py input you will get a output file named `output_matrix.txt`, a TXT file as well.

## Usage

```shell
python3 cdhit-to-matrix.py <input.txt>
```
Rename the `output_matrix.txt` as `orthologue_file.txt`, it can be the Comparative_T6SS.py input file. This script can canculate the *P/A-Value*, means the ratio of the number of proteins exist in genomes have complete T6SS and have no T6SS.

```shell
python3 Comparative_T6SS.py <orthologue_file.txt> <presence_genomes.txt> <absence_genomes.txt>
```
The `presence_genomes.txt` and the `absence_genomes.txt` are the files contain which genomes contain complete T6SS and no T6SS. When you found all the T6SS-associated orthologous groups, you can find all the proteins connected in genomics through `findconnected.py`. You just need to put the proteins in a file as input then the output will show on the screen.
```shell
python3 findconnected.py <positive_ID>
```

You can find all the example files in Datasets documents.
