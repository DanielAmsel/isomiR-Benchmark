## isomiR tools - benchmark supportive scripts
The available scripts can be used to simulate sequencing reads with a mature microRNA dataset as reference, taking the individual lengths into account.
Furthermore, one can simulate microRNA isoforms, using a reference mature microRNA dataset. <br />

With the other scripts, one can evaluate the performance (TP,FP,FN) of isomiRID, isomiR-SEA and miraligner.<br />

If you find the scripts useful, please cite:<br />



## Simulating small RNA sequencing runs
The script `create_sequencing.pl` can be used to simulate a microRNA sequencing run.
#### Dependencies : 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/)
#### Installation : 
Since it is a pure perl file, no installation is needed. Tested on Linux. 
#### Usage : 
`./create_sequencing.pl ` <br/>
#### Details :
* The standard path of ART is specified as `/opt/art_bin_MountRainier/art_illumina` for illumina sequencing. If you whish to change the path or switch to another sequencing platform, you need to change the path at the beginnging of the source code
* The script is designed in a way that it only creates sequenced reads for microRNAs within a length between 17 and 30nt. If this does not fit your experimental design, you can change the array <br /> `@mir_len = (17,18,19,20,21,22,23,24,25,26,27,28,29,30);`<br /> by the lengths you may need or focus on.

## Simulating microRNA isoforms (isomiRs)
The script `create_isomiRs.pl` is designed to simulate equally distributed microRNA isoforms.
* 5' template additions (1-3)
* 5' deletions (1-3)
* 3' template additions (1-3)
* 3' deletions (1-3)
* 3' non-template additions (1-3)
* SNPs in seed region
* SNPs outside of seed region
#### Dependencies : 
none
#### Installation :
Since it is a pure perl file, no installation is needed. Tested on Linux. 
#### Usage :
`./create_isomiRs.pl <mature_mirs.fasta> <hairpin_mirs.fasta>` <br/>
`<mature_mirs.fasta>`  := The mature microRNA set of your species of interest  <br/>
`<hairpin_mirs.fasta>` := The hairpin microRNA set of your species of interst  <br/>
#### Output :
|output file                          |        description                                 |
| ---                                 | ---                                                |
|`mature_mirs.fasta-FIVE_ADD.fa`      | template 5' additions (1-3nt)  <br/>               |
|`mature_mirs.fasta-THREE_ADD.fa`     | template 3' additions (1-3nt)  <br/>               |
|`mature_mirs.fasta-THREE_DEL.fa`     | terminal 3' deletions (1-3nt)  <br/>               |
|`mature_mirs.fasta-FIVE_DEL.fa`      | terminal 5' deletions (1-3nt)  <br/>               |
|`mature_mirs.fasta-SNP_SEED.fa`      | SNPs in seed region (1-8)  <br/>                   |
|`mature_mirs.fasta-SNP_REST.fa`      | SNPs outside the seed region (9-end)  <br/>        |
|`mature_mirs.fasta-NON_TEMPLATE.fa`  | non-template 3' additions (1-3nt)  <br/>           |


## Evaluate isomiR mining tools

### Tool 1: eval_isomiRID_0.3.pl


### Tool 2: eval_miraligner_0.3.pl


### Tool 3: eval_isomiR-SEA_0.3.pl


