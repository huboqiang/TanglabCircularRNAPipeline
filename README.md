TanglabCircularRNAPipeline
==========================
This Pipeline process data from raw fastq data to circular RNA related results.

The following software were used:
* [TopHat](http://cufflinks.cbcb.umd.edu) 2.0.9
* [Cufflinks](http://cufflinks.cbcb.umd.edu) 2.1.1
* [TopHat-Fusion](http://ccb.jhu.edu/software/tophat/fusion_index.html) included in TopHat 2.0.9
* [bedtools](https://github.com/arq5x/bedtools2)
* [SAMtools](http://samtools.sourceforge.net)
* [CIRCExplorer](http:////github.com/Yanglab/CIRCexplorer/)
Running:
```bash
  python run\_cirRNA.py [options] sample\_lists 
```
