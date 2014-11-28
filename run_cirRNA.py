from __future__ import division
import re,sys,os
import module00_preprocessing       as m00
import module01_mapping_from_cln    as m01
import module02_RNA_Quantification  as m02
import module03_cirRNA_analysis     as m03
import module04_Stat_Sample         as m04
import cPickle as pickle
import numpy as np
from optparse   import OptionParser


class DirSystem(dict):
   def __init__(self):
      home_dir       =  os.path.abspath('./')

      self['dir'] = {                                                       \
         'raw_data'              : "%s/00.raw_data"         % (home_dir),   \
         'raw_split'             : "%s/00.raw_split"        % (home_dir),   \
         'clean_data'            : "%s/01.clean_data"       % (home_dir),   \
         'trim_data'             : "%s/02.trim_data"        % (home_dir),   \
         'tophat_dir'            : "%s/03.tophat"           % (home_dir),   \
         'tophat_fusion'         : "%s/04.tophat_fusion"    % (home_dir),   \
         'HTSeq_result_dir'      : "%s/05.HTSeq_result"     % (home_dir),   \
         'HTSeq_known_dir'       : "%s/05.1.HTSeq_known"    % (home_dir),   \
         'HTSeq_unknown_dir'     : "%s/05.2.HTSeq_unknown"  % (home_dir),   \
         'cufflinks_unknown_dir' : "%s/06.cufflinks_unknown"% (home_dir),   \
         'cuffquant_dir'         : "%s/07.cuffquant"        % (home_dir),   \
         'cuffnorm_dir'          : "%s/08.cuffnorm"         % (home_dir),   \
         'CIRCexplorer_dir'      : "%s/09.CIRCexplorer"     % (home_dir),   \
      }      

def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s --index_file /datc/huboqiang/cir_dyj_V2/sample_list.20141008.less.xls --index_list /datc/huboqiang/cir_dyj_V2/index_list.xls --bt_index_base /data/Analysis/huboqiang/database/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome_SpikeIn_ERCC_RGCPloyA --genome_gtf /datc/huboqiang/cir_dyj_V2/Database/refseqGene.ERCC_RGCPloyA.exon.sort.gtf --genome_gtf_with_lncRNA  /datc/huboqiang/cir_dyj_V2/Database/refGene_plus_NoncodeV4_nsmb/Early_Embroy_LncRNA_Pool.ERCC_RGCPolyA.sort.gtf --intragenic_bed /datc/huboqiang/cir_dyj_V2/Database/region.Intragenic.bed   --refGene_hg19 /datc/huboqiang/cir_dyj_V2/Database/ref.sort.txt   --genome_hg19 /datc/huboqiang/cir_dyj_V2/Database/hg19.fa --ercc_info /datc/huboqiang/cir_dyj_V2/Database/ercc.info.xls --cir_min_depth 2 sample_list.20141104.xls
   
   """ % (sys.argv[0],sys.argv[0])

   description = " Split fastq file by the barcode "

   optparser = OptionParser(version="%s v0.2 20140528" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-f", "--index_file",    default="/datc/huboqiang/cir_dyj_V2/sample_list.20141008.less.xls", help="\nPhreads uses, [default: %default]")
   optparser.add_option("-i", "--index_list",    default="/datc/huboqiang/cir_dyj_V2/index_list.xls",                help="\nRead length [default: %default]")
   optparser.add_option("-b", "--bt_index_base", default="/data/Analysis/huboqiang/database/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome_SpikeIn_ERCC_RGCPloyA",help="\nReference genome [default: %default]")
   optparser.add_option("-g", "--genome_gtf", default="/datc/huboqiang/cir_dyj_V2/Database/refseqGene.ERCC_RGCPloyA.exon.sort.gtf",help="\nGTF for transcriptome file [default: %default]")
   optparser.add_option("-G", "--genome_gtf_with_lncRNA", default="/datc/huboqiang/cir_dyj_V2/Database/refGene_plus_NoncodeV4_nsmb/Early_Embroy_LncRNA_Pool.ERCC_RGCPolyA.sort.gtf",help="\nGTF for transcriptome file with lncRNA annotation [default: %default]")
   optparser.add_option("--intragenic_bed", default="/datc/huboqiang/cir_dyj_V2/Database/region.Intragenic.bed",help="\nIntragenic region          [default: %default]")
   optparser.add_option("--refGene_hg19",   default="/datc/huboqiang/cir_dyj_V2/Database/ref.sort.txt"         ,help="\nhg19 refGene file          [default: %default]")
   optparser.add_option("--genome_hg19",    default="/datc/huboqiang/cir_dyj_V2/Database/hg19.fa"              ,help="\nhg19 genome file           [default: %default]")
   optparser.add_option("--ercc_info",      default="/datc/huboqiang/cir_dyj_V2/Database/ercc.info.xls"        ,help="\nERCC information           [default: %default]")
   optparser.add_option("--cir_min_depth",  default=2                                                          ,help="\nMin cirRNA-depth required  [default: %default]")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser

def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      samp_info      = args[0]
      index_file     = options.index_file
      index_list     = options.index_list
      bt_index_base  = options.bt_index_base
      genome_gtf     = options.genome_gtf
      genome_gtf_wl  = options.genome_gtf_with_lncRNA
      intragenic_bed = options.intragenic_bed
      refGene_hg19   = options.refGene_hg19
      genome_hg19    = options.genome_hg19
      ercc_info      = options.ercc_info
      cir_min_depth  = int(options.cir_min_depth)
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   dir_name = DirSystem()
   
   # pre_mapping process (done 20141104):
#   samp_prepare = m00.PreMapping( index_file,index_list,dir_name['dir'] )
#   samp_prepare.load_lists()
#   samp_prepare.load_index_seq()
#   samp_prepare.run_indexSplit()
#   samp_prepare.run_QC()
   
#   samp_process = m01.Map_From_cln( samp_info,bt_index_base,genome_gtf,dir_name['dir'] )
#   samp_process.load_samp()
#   samp_process.run_trim()
#   samp_process.run_tophat()
#   samp_process.run_tophat_fusion()
   
#   samp_mRNAQ = m02.RnaQuantification( samp_info,bt_index_base,genome_gtf_wl,intragenic_bed,dir_name['dir'] )
#   samp_mRNAQ.load_samp()
#   samp_mRNAQ.RNA_QuantPipe()

#   samp_cirRNA = m03.cirRNA(           samp_info,genome_hg19 ,refGene_hg19,dir_name['dir'] )
#   samp_cirRNA.run_CIRCexplorer()
   
   samp_Stat   = m04.SampStat( samp_info,ercc_info, genome_gtf,cir_min_depth ,dir_name['dir'] )
   samp_Stat.load_samp()
   samp_Stat.CIRC_Stat()
#   samp_Stat.Basic_Stat()
   samp_Stat.ERCC_Stat()
#   samp_Stat.plot_regression()
#   samp_Stat.CIRC_Stat()
   
if __name__ == '__main__':
   main()
