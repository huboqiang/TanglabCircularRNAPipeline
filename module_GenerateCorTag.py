from __future__ import division
import re,sys,os
import module_StatInfo as Stat
import scipy
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import itertools

def delta_sort(val):
   out_value = 0
   if len(val) > 1:
      for i in range( 0,len(val)-1 ):
         out_value += val[i+1] - val[i]
   return out_value

class CorTag(object):
   '''
   Generate tissule_level tags and its 0/1 correlation matrix
   '''
   
   def __init__( self,samp_info ):
      self.l_samp       = samp_info.sample
      self.sam_info     = samp_info.sam_info
      self.l_stage_name = samp_info.stage
      self.tag_cor      = {}
      self.l_tags       = []
      
   def generate_tags(self):
      for stage in self.l_stage_name['name']:
         for i,samp in enumerate(self.sam_info['stage_sam'][stage]):
            self.__load_tag( stage,i+1 )
            
   def __load_tag( self,stage,i ):
      tag = "%s_%d" % ( stage,i )
      self.tag_cor[ tag ] = []
      l_name = self.sam_info['stage_sam'][ stage ]
      stage_beg_idx = self.l_samp.index( l_name[  0 ] )
      stage_end_idx = self.l_samp.index( l_name[ -1 ] )
      
      l_idx_comb  = list( itertools.combinations( range( stage_beg_idx,stage_end_idx+1 ), i ) )
      l_idx_delta = [ delta_sort(val) for val in l_idx_comb ]
      
      np_idx_comb = np.array( l_idx_comb  )
      np_idx_delta= np.array( l_idx_delta )
      
#      print stage_beg_idx,stage_end_idx,l_idx_comb
      np_idx_comb_sort = np_idx_comb[ np_idx_delta.argsort() ]
      
      for idx_comb in np_idx_comb_sort:
         l_index = np.zeros( len(self.l_samp),dtype=int )
#         print idx_comb
         l_index[ idx_comb ] = 1
         self.tag_cor[ tag ].append( list(l_index) )
#         print tag, l_index
      self.l_tags.append( tag )
         