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

class cirRNA(dict):
   
   def __init__( self, samp_info,genome_file,ref_file,  dir_name ):

      self['infile'] = { 'info_file':samp_info,'genome_file':genome_file,'ref_file':ref_file }
      self['dir_name']  = dir_name
      self['sam_info']  = {}
      self['samp']      = []
      self.__load_name()
      self.__load_samp()

   def __load_samp(self):      
      self['stage'] = { 'name':[] }
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

         self['samp'].append( samp )
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
      
   def __load_name(self):
      self.tophat       = self['dir_name']['tophat_dir']
      self.tophat_fusion= self['dir_name']['tophat_fusion']
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
   
   def run_CIRCexplorer(self):
      sh_file      = "%s/s10.CIRCexplorer.sh"      % (self.script_dir)
      sh_work_file = "%s/s10.CIRCexplorer_work.sh" % (self.script_dir)
      
      if not os.path.isdir( self.CIRCexplorer ):
         os.mkdir( self.CIRCexplorer )
      
      py_CIRCexplorer         = "/data/Analysis/huboqiang/software/CIRCexplorer/CIRCexplorer_PE.py"
      py_CIRCexplorer_PE_Check= "/datc/huboqiang/cir_dyj_V2/bin/CIRCexplorer_PE_check.py"
      sh_info = """
py_CIRCexplorer=$1
in_bam=$2
genome=$3
ref_file=$4
out_file=$5
samp=$6

py_CIRCexplorer_PE_check=$7
in_raw_bam=$8

[ ! -d $out_file/$samp ] && mkdir -p $out_file/$samp

#python $py_CIRCexplorer                                                 \\
#   -f $in_bam                                                           \\
#   -g $genome                                                           \\
#   -r $ref_file                                                         \\
#   --tmp                                                                \\
#   -o $out_file/$samp/CIRCexplorer
   
python $py_CIRCexplorer_PE_check                                        \\
   --raw_bam      $in_raw_bam                                           \\
   --out_prefix   $out_file/$samp/CIRCexplorer_circ_PE                  \\
   $out_file/$samp/CIRCexplorer_circ.txt
      """
      sh_work = ""
      for samp in self['samp']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         in_bam      = "%s/%s/accepted_hits.bam"   % ( self.tophat_fusion, brief_name )
         genome      =  self['infile']['genome_file']
         ref_file    =  self['infile']['ref_file']
         out_file    = self.CIRCexplorer
         in_raw_bam  = "%s/%s/accepted_hits.bam"   % ( self.tophat, brief_name )
         
         sh_work += "sh %s  %s %s %s %s %s %s   %s %s \n" % ( sh_file,  py_CIRCexplorer, in_bam, genome,  ref_file, out_file, brief_name,  py_CIRCexplorer_PE_Check,in_raw_bam )
      
      my_job = m_jobs.running_jobs(sh_file,sh_work_file)
      my_job.load_sh_file(      sh_info )
      my_job.load_sh_work_file( sh_work )
#      my_job.running_SGE( vf="1000m",maxjob=100 )
