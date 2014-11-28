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
import module_RNAQuantPipe as m_cuff
import module_GTFFeature   as m_gtfF



class RnaQuantification(dict):
   def __init__(self,samp_info,genome_file,anno_file,intragenic_bed,dir_name ):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info, 'genome_file':genome_file,'intragenic_bed':intragenic_bed,'anno_file':anno_file }
      self['stage'] = { 'name':[] }
      self['dir_name'] = dir_name

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
   
   def RNA_QuantPipe(self):
      new_trans = m_cuff.QuantPipe( self['sam_info'], self['sample'], self['infile']['genome_file'], self['infile']['anno_file'], self['infile']['intragenic_bed'], self['dir_name'] )
      new_trans.run_HTSeq_known()
      new_trans.run_cufflinks_u()
      new_trans.run_cuffcomp_novo_trans()
      new_trans.run_HTSeq_unknown()
#      new_trans.RPKM_novo_trans()
      new_trans.run_cuffquant()
      new_trans.run_cuffnorm()

      
      
      
      
      