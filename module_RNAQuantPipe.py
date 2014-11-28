#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle
import scipy   as sp
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import stats 
import module_running_jobs as m_jobs
import module_GTFFeature   as m_gtf
import module_CountInfo    as m_cnt

class QuantPipe(dict):
   def __init__( self,M_samp_info,l_samp,genome_file,anno_file,intragenic_bed,  dir_name ):
      self['samp_info'] = M_samp_info
      self['samp']      = l_samp
      self['infile'] = { 'genome_file':genome_file,'anno_file':anno_file,'intragenic_bed':intragenic_bed }
      self['dir_name']  = dir_name
      self.__load_name()
      
   def __load_name(self):
      self.tophat       = self['dir_name']['tophat_dir']
      self.HTS          = self['dir_name']['HTSeq_result_dir']
      self.HTS_k        = self['dir_name']['HTSeq_known_dir']
      self.HTS_u        = self['dir_name']['HTSeq_unknown_dir']
      self.cufflink_u   = self['dir_name']['cufflinks_unknown_dir']
      self.cuffquant    = self['dir_name']['cuffquant_dir']
      self.cuffnorm     = self['dir_name']['cuffnorm_dir']
      self.CIRCexplorer = self['dir_name']['CIRCexplorer_dir']
      home_dir          = os.path.abspath('./')
      self.script_dir   = "%s/scripts" % (home_dir)
      self.bin_dir      = "%s/bin"     % (home_dir)
      self.data_dir     = "%s/Database"% (home_dir)

   def run_HTSeq_known(self):
      
      sh_file        =  "%s/s04.HTSeq_known.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/s04.HTSeq_known_work.sh"  %  (self.script_dir)
      
      py_deseq       =  "/data/Analysis/huboqiang/bin/htseq-count"
      
      sh_info = """
tophat_dir=$1
samp_name=$2
py_deseq=$3
HTS_k_dir=$4
known_GTF=$5

samtools view  -H              $tophat_dir/$samp_name/accepted_hits.bam  > $tophat_dir/$samp_name/accepted_hits.header.sam
samtools sort  -n -m 200000000 $tophat_dir/$samp_name/accepted_hits.bam    $tophat_dir/$samp_name/accepted_hits.sort_name
samtools view  -o    $tophat_dir/$samp_name/accepted_hits.sort_name.sam    $tophat_dir/$samp_name/accepted_hits.sort_name.bam 

[ ! -d $HTS_k_dir/$samp_name ] && mkdir -p $HTS_k_dir/$samp_name

python $py_deseq                                                                                                                             \\
   -s no -f sam -a 10                                                                                                                        \\
   -o $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam                                                                                \\
   $tophat_dir/$samp_name/accepted_hits.sort_name.sam  $known_GTF            >$HTS_k_dir/$samp_name/$samp_name.dexseq.txt                 && \\
grep -v -P '^ERCC-|^RGC-|MIR|SNORD|Mir|Snord' $HTS_k_dir/$samp_name/$samp_name.dexseq.txt > $HTS_k_dir/$samp_name/$samp_name.dexseq_clean.txt			&& \\
grep    -P '^ERCC-|^RGC-'                     $HTS_k_dir/$samp_name/$samp_name.dexseq.txt > $HTS_k_dir/$samp_name/$samp_name.dexseq_ERCC_RGCPloyA.txt	&& \\
grep "__no_feature" $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam | grep -v chrM |                                                 \\
	cat           $tophat_dir/$samp_name/accepted_hits.header.sam /dev/stdin | 		                                                         \\
	samtools view -Sb /dev/stdin >$tophat_dir/$samp_name/accepted_hits.genome.bam                                                          && \\
samtools sort  -m 200000000 $tophat_dir/$samp_name/accepted_hits.genome.bam    $tophat_dir/$samp_name/accepted_hits.genome.sort          
rm             $tophat_dir/$samp_name/accepted_hits.sort_name.sam $tophat_dir/$samp_name/accepted_hits.sort_name.bam   $tophat_dir/$samp_name/accepted_hits.sort_name.gene.sam  $tophat_dir/$samp_name/accepted_hits.genome.bam
      """
      sh_work = ""
      for samp in self['samp']:
         tophat_dir  =  self.tophat
         samp_name   =  self['samp_info']['samp_brief'][samp]
         known_GTF   =  self['infile']['anno_file']

         sh_work += "sh %s  %s %s %s %s %s\n" % ( sh_file,  self.tophat, samp_name, py_deseq, self.HTS_k, known_GTF)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="1g",maxjob=100 )

   
   ### 
   def run_cufflinks_u(self):

      sh_file      = "%s/s05.cufflinks_GenomeMapped.sh"      % (self.script_dir)
      sh_work_file = "%s/s05.cufflinks_GenomeMapped_work.sh" % (self.script_dir)
      
      sh_info = """
in_bam=$1
gtf_file=$2
out_dir=$3

/data/Analysis/huboqiang/software/cufflinks-2.2.1.Linux_x86_64/cufflinks      \\
   -p 8  -u                                                                   \\
   -o $out_dir                                                                \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.genome.sort.bam"% ( self.tophat,    brief_name   ) 
         out_dir     = "%s/%s"                              % ( self.cufflink_u,brief_name   )
         sh_work += "sh %s  %s %s %s \n" % ( sh_file,  in_bam, self['infile']['anno_file'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="500m",maxjob=100 )

   def run_cuffcomp_novo_trans(self):

      sh_file      = "%s/s06.1.cuffcompare_novo.sh"      % (self.script_dir)
      sh_work_file = "%s/s06.1.cuffcompare_novo_work.sh" % (self.script_dir)
      
      sh_info = """
out_prefix=$1
shift 

/data/Analysis/huboqiang/software/cufflinks-2.2.1.Linux_x86_64/cuffcompare    \\
   -o  $out_prefix                                                            \\
   -T  $@                                                                     \\
      """
      sh_work = ""
      out_prefix  = "%s/novo_lnc_raw"           % ( self.data_dir )
      l_in_samp   = [  "%s/%s/transcripts.gtf"  % ( self.cufflink_u,self['samp_info']['samp_brief'][samp] ) for samp in self['samp'] ]
      sh_work = "sh %s   %s %s" % ( sh_file, out_prefix, " ".join(l_in_samp) )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="500m",maxjob=100 )


   def run_HTSeq_unknown(self):
      
      sh_file        =  "%s/s07.HTSeq_unknown.sh"       %  (self.script_dir)
      sh_work_file   =  "%s/s07.HTSeq_unknown_work.sh"  %  (self.script_dir)
      
      py_deseq       =  "/data/Analysis/huboqiang/bin/htseq-count"
      
      if not os.path.isdir( self.HTS_u ):
         os.mkdir( self.HTS_u )
      
      sh_info = """
tophat_dir=$1
samp_name=$2
py_deseq=$3
HTS_u_dir=$4
unknown_GTF=$5

samtools view  -H              $tophat_dir/$samp_name/accepted_hits.genome.sort.bam          > $tophat_dir/$samp_name/accepted_hits.header.sam
samtools sort  -n -m 200000000 $tophat_dir/$samp_name/accepted_hits.genome.sort.bam            $tophat_dir/$samp_name/accepted_hits.genome.sort_name
samtools view  -o              $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam  $tophat_dir/$samp_name/accepted_hits.genome.sort_name.bam 

[ ! -d $HTS_u_dir/$samp_name ] && mkdir -p $HTS_u_dir/$samp_name

python $py_deseq                                                                                                                             \\
   -s no -f sam -a 10                                                                                                                        \\
   $tophat_dir/$samp_name/accepted_hits.genome.sort_name.sam  $unknown_GTF            >$HTS_u_dir/$samp_name/$samp_name.dexseq_NeoRaw.txt
      """
      sh_work = ""
      for samp in self['samp']:
         tophat_dir  =  self.tophat
         samp_name   =  self['samp_info']['samp_brief'][samp]
         unknown_GTF =  "%s/novo_lnc_raw.combined.gtf" % ( self.data_dir )
         sh_work += "sh %s  %s %s %s %s %s\n" % ( sh_file,  self.tophat, samp_name, py_deseq, self.HTS_u, unknown_GTF)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="1g",maxjob=100 )
   
   def RPKM_novo_trans(self):
      l_brief_samp    = [ "%s" % ( self['samp_info']['samp_brief'][samp] ) for samp in self['samp'] ]
      unknown_GTF     = "%s/novo_lnc_raw.combined.gtf" % ( self.data_dir )
      
      Gtf_Info = m_gtf.GTFFeature( unknown_GTF  )
      Gtf_Info.get_intergenic( self['infile']['intragenic_bed'] )
      Cnt_Info = m_cnt.CountInfo(  self.HTS_u, l_brief_samp, "dexseq_NeoRaw", self.HTS  )

      Cnt_Info.generate_mat()
      Cnt_Info.load_mat()
      Cnt_Info.cal_RPKM( Gtf_Info.gene,self.tophat )
      
      rpkm_file = "%s/merge.%s.RPKM.xls" % ( self.HTS,"dexseq_NeoRaw" )
      
      Gtf_Info.load_gene_RPKM( rpkm_file )
      Gtf_Info.output_GTF()
      Gtf_Info.get_gene_info()
      
   def run_cuffquant(self):
      sh_file      = "%s/s08.cuffquant.sh"      % (self.script_dir)
      sh_work_file = "%s/s08.cuffquant_work.sh" % (self.script_dir)
      
      if not os.path.isdir( self.cuffquant ):
         os.mkdir( self.cuffquant )
      
      sh_info = """
in_bam=$1
gtf_file=$2
out_dir=$3

/data/Analysis/huboqiang/software/cufflinks-2.2.1.Linux_x86_64/cuffquant      \\
   -p 8  -u                                                                   \\
   -o $out_dir                                                                \\
   $gtf_file                                                                  \\
   $in_bam
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat,    brief_name   ) 
         out_dir     = "%s/%s"                     % ( self.cuffquant ,brief_name   )
         sh_work += "sh %s  %s %s %s \n" % ( sh_file,  in_bam, self['infile']['anno_file'], out_dir)
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="1000m",maxjob=100 )      


   def run_cuffnorm(self):
      sh_file      = "%s/s09.cuffnorm.sh"      % (self.script_dir)
      
      if not os.path.isdir( self.cuffnorm ):
         os.mkdir( self.cuffnorm )

      l_brief = []
      l_cxb   = []
      for samp in self['samp']:
         brief_name = self['samp_info']['samp_brief'][samp]
         l_brief.append( brief_name  )
         l_cxb.append(   "%s/%s/abundances.cxb" % (self.cuffquant,brief_name) )
      
      
      list_brief = ",".join( l_brief )
      list_cxb   = " ".join( l_cxb   )
      
      sh_info = """
/data/Analysis/huboqiang/software/cufflinks-2.2.1.Linux_x86_64/cuffnorm        \\
   -p 8  -o %s  -L %s     \\
   %s                     \\
   %s
      """ % ( self.cuffnorm, list_brief, self['infile']['anno_file'], list_cxb )
      
      f_sh_file = open( sh_file,"w" )
      print >>f_sh_file, sh_info
      f_sh_file.close()
      
      shell_info = "sh %s" % ( sh_file )
#      p = subprocess.Popen(shell_info,shell='True')
#      while 1:
#         run_cnt = 0
#         if p.poll() is None:
#            run_cnt += 1
#            time.sleep(1)
#         if run_cnt == 0:
#            break
   
   
      
      