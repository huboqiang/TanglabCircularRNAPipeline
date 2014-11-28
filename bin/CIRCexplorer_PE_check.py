import re,sys,os
import cPickle as pickle
import subprocess,time
import numpy as np
from optparse   import OptionParser
import sys
import pysam

def prepare_optparser():
   usage ="""usage: %s [options] 

   Using -h or --help for more information

Example:
   python %s --raw_bam /datc/huboqiang/cir_dyj_V2/03.tophat/2cell_1_B/accepted_hits.bam  --out_prefix /datc/huboqiang/cir_dyj_V2/09.CIRCexplorer/2cell_1_B/CIRCexplorer_circ_PE /datc/huboqiang/cir_dyj_V2/09.CIRCexplorer/2cell_1_B/CIRCexplorer_circ.txt
   
   """ % (sys.argv[0],sys.argv[0])

   description = " Get the Pair-End circular-RNA information from CIRExplorer data. When detected a cirRNA, scanning the region between the two break-points, and using its qname-ID to see if this read's PE-reads could map to genome linear and completely."

   optparser = OptionParser(version="%s v0.1 20141117" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
   optparser.add_option("-f", "--raw_bam",    default="/datc/huboqiang/cir_dyj_V2/03.tophat/2cell_1_B/accepted_hits.bam",       help="\nPhreads uses, [default: %default]")
   optparser.add_option("-o", "--out_prefix", default="/datc/huboqiang/cir_dyj_V2/09.CIRCexplorer/2cell_1_B/CIRCexplorer_circ_PE",help="\nOutput file prefix[default: %default]")
   optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
   return optparser

def main():
   prepare_optparser()
   (options,args) = prepare_optparser().parse_args()
   try:
      cir_info    = args[0]
      raw_bam     = options.raw_bam     
      out_prefix  = options.out_prefix
   except IndexError:
      prepare_optparser().print_help()
      sys.exit(1)
   
   if not os.path.isfile( "%s.bai" % (raw_bam) ):
      shell_info = "samtools index %s" % (raw_bam)
      print >>sys.stderr, shell_info
      p = subprocess.Popen(shell_info,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(3)
         if run_cnt == 0:
            break
   
   f_sam = pysam.Samfile( raw_bam,"rb" )
   f_cir = open( cir_info,"r" )
   
   PE_circ        = "%s.txt"              % (out_prefix)
   PE_not_pass    = "%s.PE_not_pass.txt"  % (out_prefix)
   PE_log         = "%s.pass.log"         % (out_prefix)
   
   f_PE_circ      = open( PE_circ    ,"w" )
   f_PE_not_pass  = open( PE_not_pass,"w" )
   f_PE_log       = open( PE_log     ,"w" )
   
   total_circ   = 0
   pass_PE_circ = 0
   notP_PE_circ = 0

   total_read   = 0
   pass_PE_read = 0
   notP_PE_read = 0
   
   for line in f_cir:
      line  =  line.strip('\n')
      f     =  line.split()
      chrom =      f[0]
      beg   =  int(f[1])
      end   =  int(f[2])
      
      read_cnt    =  int(f[12])

      total_circ += 1
      total_read += read_cnt

      l_read_qname=  f[17].split(",")

      PE_read_cnt    =  0
      l_PE_read_qname=  []
      reads = f_sam.fetch( reference=chrom, start=beg, end=end )
      for r_line in reads:
         for qname in l_read_qname:
            if qname == r_line.qname:
               PE_read_cnt += 1
               l_PE_read_qname.append( qname )
               l_read_qname.remove( qname )
               continue

      if PE_read_cnt == read_cnt:
         pass_PE_circ += 1
         f[12] = str( PE_read_cnt )
         f[17] = ",".join( l_PE_read_qname )
         out_pass = "\t".join( f )         
         print >>f_PE_circ, out_pass
      else:
         if PE_read_cnt == 0:
            notP_PE_circ += 1
            f[12] = str( read_cnt - PE_read_cnt )
            f[17] = ",".join( set(l_read_qname) - set(l_PE_read_qname) )
            out_not_pass = "\t".join( f )
            print >>f_PE_not_pass, out_not_pass
         else:
            pass_PE_circ += 1
            
            f[12] = str( read_cnt - PE_read_cnt )
            f[17] = ",".join( set(l_read_qname) - set(l_PE_read_qname) )
            out_not_pass = "\t".join( f )
            print >>f_PE_not_pass, out_not_pass
            
            f[12] = str( PE_read_cnt )
            f[17] = ",".join( set(l_PE_read_qname) )
            out_pass = "\t".join( f )
            print >>f_PE_circ, out_pass
      
      pass_PE_read += PE_read_cnt
      notP_PE_read += read_cnt - PE_read_cnt
   
   print >>f_PE_log, "Total_circ\tTotal_pass_PE_check_circ\tTotal_notPass_PE_check_circ\tTotal_CircRead\tTotal_pass_PE_check_circ\tTotal_notPass_PE_check_circ\n%d\t%d\t%d\t%d\t%d\t%d" % \
   ( total_circ,pass_PE_circ,notP_PE_circ,  total_read,pass_PE_read,notP_PE_read )
   
   f_cir.close()
   f_sam.close()
               
   f_PE_circ.close()
   f_PE_not_pass.close()
   f_PE_log.close()

if __name__ == '__main__':
   main()
