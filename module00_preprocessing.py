from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import time,subprocess


class PreMapping(dict):
   def __init__( self,index_file,index_list,dir_name ):
      self['infile'] = { "index":index_file, "lists":index_list }
      self['index_seq'] = {}
      self['samp_info'] = { 'in_fq':[],'fq_sam':{},'fq_idx':{} }
      self['dir_name']  = dir_name
      
   def load_lists(self):
      f_lists = open( self['infile']["lists"],"r" )
      for line in f_lists:
         line = line.strip('\n')
         f    = line.split()
         self['index_seq'][ f[0] ] = f[1]
      f_lists.close()
      
   def load_index_seq(self):   

      home_dir = os.path.abspath('./')
      f_index = open(  self['infile']["index"],"r"  )

      f_index.readline()
      for line in f_index:
         line = line.strip('\n')
         f    = line.split()
         samp = f[0]
         idx  = f[1]
         in_fq= f[2]
         
         l_in_fq = [ "%s/00.raw_data/%s" % (home_dir,sam) for sam in in_fq.split(',') ]
         in_fq = ",".join( l_in_fq )
         
         if in_fq not in self['samp_info']['in_fq']:
            self['samp_info']['in_fq'].append( in_fq )
         
         if in_fq not in self['samp_info']['fq_sam']:
            self['samp_info']['fq_sam'][in_fq] = []
         self['samp_info']['fq_sam'][in_fq].append( samp )
         
         if samp not in self['samp_info']['fq_idx']:
            self['samp_info']['fq_idx'][ samp ] = { 'sam_fq':[],'use_idx':[] }
         
         self['samp_info']['fq_idx'][ samp ]['use_idx'] = idx.split(',')
         self['samp_info']['fq_idx'][ samp ]['sam_fq'] = [ "%s/00.raw_data/%s" % (home_dir,sam) for sam in in_fq.split(',') ]
      f_index.close()
   
   def run_indexSplit(self):
      home_dir = os.path.abspath('./')
      
      if not os.path.isdir( "%s/scripts" % (home_dir) ):
         os.mkdir(           "%s/scripts" % (home_dir) )
      
      py_splitIndex       = "%s/bin/split_index/index_split.py" % ( home_dir )
      
      indexSplit_sh       = "%s/scripts/indexSplit.sh"      % (home_dir)
      f_indexSplit_sh     = open( indexSplit_sh    ,"w" )
      run_indexSplit_sh   = "%s/scripts/indexSplit_work.sh" % (home_dir)
      f_run_indexSplit_sh = open( run_indexSplit_sh,"w" )
      
      shell_info = """
py_splitIndex=$1
in_fq1=$2
in_fq2=$3
idx_cod=$4

python $py_splitIndex -t 4 -l 125 -b 0 $in_fq1 $in_fq2 $idx_cod 
      """
      
      print >>f_indexSplit_sh, shell_info
      f_indexSplit_sh.close()
         
      for in_fq in self['samp_info']['in_fq']:
         in_fq1,in_fq2 = in_fq.split(',')
         
         l_idx = []
         for samp in self['samp_info']['fq_sam'][in_fq]:
            for idx in self['samp_info']['fq_idx'][ samp ]['use_idx']:
               l_idx.append( self['index_seq'][ idx ] )
         
         idx_cod = ",".join( l_idx )
         
         print >>f_run_indexSplit_sh, " sh %s  %s %s %s %s " % ( indexSplit_sh,   py_splitIndex, in_fq1, in_fq2, idx_cod )
      
      f_run_indexSplit_sh.close()
   
      shell_work = "perl /data/Analysis/huboqiang/bin/qsub-sge.pl --resource vf=200m --maxjob 100 %s" % ( run_indexSplit_sh )
#      p = subprocess.Popen(shell_work,shell='True')
#      while 1:
#         run_cnt = 0
#         if p.poll() is None:
#            run_cnt += 1
#            time.sleep(3)
#         if run_cnt == 0:
#            break
      
      self.__rename()



   def run_QC(self):
      home_dir = os.path.abspath('./')
      
      if not os.path.isdir( "%s/scripts" % (home_dir) ):
         os.mkdir(           "%s/scripts" % (home_dir) )
      
      clean_dir = self['dir_name']['clean_data']
      
      if not os.path.isdir( clean_dir ):
         os.mkdir(          clean_dir )
      
      pl_QC       = "%s/bin/QC.pl" % ( home_dir )
      
      QC_sh       = "%s/scripts/QC.sh"      % (home_dir)
      f_QC_sh     = open( QC_sh    ,"w" )
      run_QC_sh   = "%s/scripts/QC_work.sh" % (home_dir)
      f_run_QC_sh = open( run_QC_sh,"w" )
      
      shell_info = """
pl_QC=$1
in_dir=$2
out_dir=$3

perl $pl_QC $in_dir $out_dir 2 
      """
      print >>f_QC_sh, shell_info
      f_QC_sh.close()

      for in_fq in self['samp_info']['in_fq']:
         for samp in self['samp_info']['fq_sam'][in_fq]:
            
            if not os.path.isdir( "%s/01.clean_data/%s" % (home_dir,samp) ):
               os.mkdir(          "%s/01.clean_data/%s" % (home_dir,samp) )
            
            in_dir  = "%s/00.raw_split/%s"  % ( home_dir,samp )
            out_dir = "%s/01.clean_data/%s" % ( home_dir,samp )
            
            print >>f_run_QC_sh, " sh %s  %s %s %s  " % ( QC_sh,   pl_QC, in_dir,out_dir )
            
      f_run_QC_sh.close()
   
      shell_work = "perl /data/Analysis/huboqiang/bin/qsub-sge.pl --resource vf=500m --maxjob 100 %s" % ( run_QC_sh )
      p = subprocess.Popen(shell_work,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(3)
         if run_cnt == 0:
            break


   
   def __rename(self):                          # run_indexSplit
      
      home_dir = os.path.abspath('./')
      
      if not os.path.isdir( "%s/00.raw_split" % (home_dir) ):
         os.mkdir(          "%s/00.raw_split" % (home_dir) )
      
      rename_sh       = "%s/scripts/rename.sh"      % (home_dir)
      f_rename_sh     = open( rename_sh    ,"w" )
      
      for in_fq in self['samp_info']['in_fq']:
         
         in_fq1,in_fq2 = in_fq.split(',')
         fq1_prefix = ".".join( in_fq1.split(".")[:-2] )
         fq2_prefix = ".".join( in_fq2.split(".")[:-2] )
         print fq1_prefix, fq2_prefix
         
         for samp in self['samp_info']['fq_sam'][in_fq]:
            
            if not os.path.isdir( "%s/00.raw_split/%s" % (home_dir,samp) ):
               os.mkdir(          "%s/00.raw_split/%s" % (home_dir,samp) )
            
            mrg_samp_file_1 = "%s/00.raw_split/%s/%s.1.fq.gz" % ( home_dir,samp,samp )
            mrg_samp_file_2 = "%s/00.raw_split/%s/%s.2.fq.gz" % ( home_dir,samp,samp )
            
            l_in_file1 = []
            l_in_file2 = []
            for idx in self['samp_info']['fq_idx'][ samp ]['use_idx']:
               # "%s/00.raw_data/%s" % (home_dir,sam) for sam in in_fq.split(',')
               # R0140100433_lane1/lane1_L001_R1_001.fastq.gz,  R0140100433_lane1/lane1_L001_R2_001.fastq.gz
               in_file1 = "%s.%s.fq.gz" % ( fq1_prefix, self['index_seq'][ idx ] )
               in_file2 = "%s.%s.fq.gz" % ( fq2_prefix, self['index_seq'][ idx ] )
               l_in_file1.append( in_file1 )
               l_in_file2.append( in_file2 )
            
            if len( l_in_file1 ) == 1:
#               print >>f_rename_sh, " mv %s %s && mv %s %s" % ( " ".join(l_in_file1),mrg_samp_file_1,  " ".join(l_in_file2),mrg_samp_file_2 )
               print >>f_rename_sh, " ln -s %s %s && ln -s %s %s" % ( " ".join(l_in_file1),mrg_samp_file_1,  " ".join(l_in_file2),mrg_samp_file_2 )
            else:
#               print >>f_rename_sh, " cat %s >%s && cat %s >%s && rm %s %s" % ( " ".join(l_in_file1),mrg_samp_file_1,  " ".join(l_in_file2),mrg_samp_file_2,   " ".join(l_in_file1), " ".join(l_in_file2))
               print >>f_rename_sh, " cat %s >%s && cat %s >%s" % ( " ".join(l_in_file1),mrg_samp_file_1,  " ".join(l_in_file2),mrg_samp_file_2  )

      f_rename_sh.close()
      
         
         
         