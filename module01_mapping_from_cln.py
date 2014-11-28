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

class Map_From_cln(dict):
   
   def __init__(self,samp_info, genome_file,anno_file, dir_name ):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info, 'genome_file':genome_file,'anno_file':anno_file }
      self['stage'] = { 'name':[] }
      self['dir_name']  = dir_name
      
   def load_samp(self):
      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['type']     = {}
      self['sam_info']['stage']    = {}
      self['sam_info']['dilute']   = {}
      self['sam_info']['stage_sam']= {}
      self['sam_info']['RFP_mols'] = {}
      self['sam_info']['GFP_mols'] = {}
      self['sam_info']['CRE_mols'] = {}
      self['sam_info']['stage_sam']= {}
      
      info_file = self['infile']['info_file']
      file = open(info_file,"r")
      in_h = file.readline()
      for line in file:
         line = line.strip('\n')
         f = line.split()
         samp       = f[0]
         brief_name = f[1]
         ltype      = f[2]
         stage      = f[3]
         ERCC_dilute= float( f[4] )
         RFP_mols   = float( f[5] )
         GFP_mols   = float( f[6] )
         CRE_mols   = float( f[7] )

         self['sample'].append( samp )
         self['sam_info']['samp_brief'][ samp ] = brief_name
         self['sam_info']['type'][  samp ]      = ltype
         self['sam_info']['stage'][ samp ]      = stage
         self['sam_info']['dilute'][samp ]      = ERCC_dilute
         self['sam_info']['RFP_mols'][ samp ]   = RFP_mols
         self['sam_info']['GFP_mols'][ samp ]   = GFP_mols
         self['sam_info']['CRE_mols'][ samp ]   = CRE_mols


         
         if stage not in self['stage']['name']:
            self['stage']['name'].append( stage )
         if stage not in self['sam_info']['stage_sam']:
            self['sam_info']['stage_sam'][stage] = []
            
         self['sam_info']['stage_sam'][stage].append( samp )
      file.close()
   
   def run_trim(self):
      home_dir     = os.path.abspath('./')
      
      cln_dir      = self['dir_name']['clean_data']
      trim_dir     = self['dir_name']['trim_data']
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s01.trim.sh"      % (script_dir)
      sh_work_file = "%s/s01.trim_work.sh" % (script_dir)
      py_trim      = "%s/step1.trim.py"    % (bin_dir)
      
      sh_info = """
py_trim=$1
in_fq1=$2
in_fq2=$3
out_dir=$4
out_prefix=$5

python $py_trim $in_fq1 $in_fq2 $out_dir $out_prefix
      """
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         in_fq1      = "%s/%s/1.cln.fq.gz" % ( cln_dir ,samp )
         in_fq2      = "%s/%s/2.cln.fq.gz" % ( cln_dir ,samp )
         out_dir     = "%s"                % ( trim_dir      )
         out_prefix  = brief_name
         sh_work += "sh %s  %s %s %s %s %s\n" % ( sh_file,  py_trim, in_fq1, in_fq2, out_dir, out_prefix )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="400m",maxjob=100 )
   
   def run_tophat(self):
      home_dir     = os.path.abspath('./')
      
      cln_dir      = self['dir_name']['clean_data']
      trim_dir     = self['dir_name']['trim_data']
      tophat_dir   = self['dir_name']['tophat_dir']
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s02.tophat.sh"      % (script_dir)
      sh_work_file = "%s/s02.tophat_work.sh" % (script_dir)
      
      sh_info = """
trim_dir=$1
brief_name=$2
tophat_dir=$3
genome=$4
gtf_file=$5

/data/Analysis/huboqiang/software/tophat-2.0.12.Linux_x86_64/tophat  \\
   -a 6 --microexon-search -m 2                                      \\
   -p 8 -G $gtf_file                                                 \\
   --library-type fr-unstranded                                      \\
   --transcriptome-index /datc/huboqiang/cir_dyj_V2/Database/refseqGene.ERCC_RGCPloyA.exon.sort \\
   -o $tophat_dir/$brief_name                                        \\
   $genome                                                           \\
   $trim_dir/TRIMED_${brief_name}.1.clean.fq.gz                      \\
   $trim_dir/TRIMED_${brief_name}.2.clean.fq.gz
      """
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         sh_work += "sh %s  %s %s %s %s %s\n" % ( sh_file,   trim_dir, brief_name, tophat_dir, self['infile']['genome_file'],self['infile']['anno_file']  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="7g",maxjob=100 )
      
   def run_tophat_fusion(self):
      home_dir     = os.path.abspath('./')
      
      cln_dir      = self['dir_name']['clean_data']
      trim_dir     = self['dir_name']['trim_data']
      tophat_dir   = self['dir_name']['tophat_dir']
      fusion_dir   = self['dir_name']['tophat_fusion']
      
      script_dir   = "%s/scripts"         % (home_dir)
      bin_dir      = "%s/bin"             % (home_dir)
      
      sh_file      = "%s/s03.tophat_fusion.sh"      % (script_dir)
      sh_work_file = "%s/s03.tophat_fusion_work.sh" % (script_dir)
      
      sh_info = """
tophat_dir=$1
brief_name=$2
fusion_dir=$3
genome=$4

/data/Analysis/huboqiang/software/bedtools-2.17.0/bin/bedtools bamtofastq -i $tophat_dir/$brief_name/unmapped.bam -fq /dev/stdout | gzip - >$tophat_dir/$brief_name/unmapped.fq.gz

/data/Analysis/huboqiang/software/tophat-2.0.12.Linux_x86_64/tophat  \\
   --fusion-search --keep-fasta-order --bowtie1                      \\
   --no-coverage-search -p 8                                         \\
   -o $fusion_dir/$brief_name                                        \\
   $genome                                                           \\
   $tophat_dir/$brief_name/unmapped.fq.gz
      """
      sh_work = ""
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         sh_work += "sh %s  %s %s %s %s \n" % ( sh_file,   tophat_dir,brief_name,fusion_dir,self['infile']['genome_file']  )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
      my_job.running_SGE( vf="7g",maxjob=5 )

      
