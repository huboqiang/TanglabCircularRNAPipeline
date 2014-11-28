#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os
import subprocess
import numpy   as np
import cPickle as pickle

class LoadSamp(object):
   def __init__(self,infile):
      self.infile = infile
      self.sample = []
      self.sam_info = {}
      self.stage    = { 'name':[] }
   def load_samp(self):
      
      self.sam_info['samp_brief'] = {}
      self.sam_info['type']     = {}
      self.sam_info['stage']    = {}
      self.sam_info['dilute']   = {}
      self.sam_info['stage_sam']= {}
      self.sam_info['RFP_mols'] = {}
      self.sam_info['GFP_mols'] = {}
      self.sam_info['CRE_mols'] = {}
      self.sam_info['stage_sam']= {}
      
      info_file = self.infile
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

         self.sample.append( samp )
         self.sam_info['samp_brief'][ samp ] = brief_name
         self.sam_info['type'][  samp ]      = ltype
         self.sam_info['stage'][ samp ]      = stage
         self.sam_info['dilute'][samp ]      = ERCC_dilute
         self.sam_info['RFP_mols'][ samp ]   = RFP_mols
         self.sam_info['GFP_mols'][ samp ]   = GFP_mols
         self.sam_info['CRE_mols'][ samp ]   = CRE_mols


         
         if stage not in self.stage['name']:
            self.stage['name'].append( stage )
         if stage not in self.sam_info['stage_sam']:
            self.sam_info['stage_sam'][stage] = []
            
         self.sam_info['stage_sam'][stage].append( samp )
      file.close()