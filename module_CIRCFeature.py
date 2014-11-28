from __future__ import division
import re,sys,os

class CirFeature(object):

   def __init__(self,root_dir,l_samp,name_prefix = "CIRCexplorer_circ_PE.txt"):
      self.dir    = root_dir
      self.l_samp = l_samp
      self.name_prefix = name_prefix
      self.count_info = {}
      self.cir_info   = {}
      
      self.samp_pos   = {}
      self.gen2cir_pos= {}
      
      self.gene_count_info = {}
      
   def load_CIRC(self):
      for samp in self.l_samp:
         infile = "%s/%s/%s" % ( self.dir,samp,self.name_prefix )
         f_infile= open( infile,"r" )
         for line in f_infile:
            line = line.strip('\n')
            f    = line.split()

            '''
            ignore ciRNA (Intron only circular RNA)
            '''
            if f[13] == "Yes":
               continue
            
            self.__load_feature( line )
            self.__load_count( samp,line )
         f_infile.close()
         
   def __load_feature(self,line):
      f = line.split()
      chrom = f[0]
      beg   = int(f[1])
      end   = int(f[2])
      strand= f[5]
      pos_info = "%s:%d-%d,%s" % ( chrom,beg,end,strand )
      
      idx_realign = f[4]
      exon_cnt    = int( f[9] )
      l_exon_size = f[10].split(',')
      l_exon_ofst = f[11].split(',')
      gene_name   = f[14]
      intron_pos  = f[16]
      
      l_exon_size = [ int(exon_size) for exon_size in l_exon_size ]
      l_exon_ofst = [ int(exon_ofst) for exon_ofst in l_exon_ofst ]
      
      if pos_info not in self.cir_info:
         self.cir_info[ pos_info ] = {       \
            'realign'      :  idx_realign,   \
            'exon_cnt'     :  exon_cnt,      \
            'l_exon_size'  :  l_exon_size,   \
            'l_exon_ofst'  :  l_exon_ofst,   \
            'gene_name'    :  gene_name,     \
            'intron_pos'   :  intron_pos     }
      if gene_name not in self.gen2cir_pos:      
         self.gen2cir_pos[ gene_name ] = []
      if pos_info not in self.gen2cir_pos[ gene_name ]:
         self.gen2cir_pos[ gene_name ].append( pos_info )
   
   def __load_count(self,samp,line):
      if samp not in self.samp_pos:
         self.samp_pos[samp] = []
      f = line.split()
      chrom = f[0]
      beg   = int(f[1])
      end   = int(f[2])
      strand= f[5]
      cir_cnt=int(f[12])
      pos_info = "%s:%d-%d,%s" % ( chrom,beg,end,strand )
      if pos_info not in self.count_info:
         self.count_info[ pos_info ] = {}

      self.count_info[ pos_info ][ samp ] = cir_cnt
      if pos_info not in self.samp_pos[ samp ]:
         self.samp_pos[ samp ].append( pos_info )
      
   def merge_gene_count(self):
      for gene in self.gen2cir_pos:
         
         if gene not in self.gene_count_info:
            self.gene_count_info[gene] = {}
         
         for pos_info in self.gen2cir_pos[gene]:
            for samp in self.l_samp:
               
               if samp not in self.gene_count_info[gene]:
                  self.gene_count_info[gene][samp] = []
               
               if samp in self.count_info[pos_info]:
                  cir_cnt = self.count_info[pos_info][samp]
                  self.gene_count_info[gene][samp].append( cir_cnt )
                  
   def filter_min_depth(self,pos_info,min_depth):
      """
      Return a bool value for cirRNA positions for whether there is a sample have more than (min_depth) reads. 
      """
      is_pass = 0
      for samp in self.count_info[ pos_info ]:
         if self.count_info[pos_info][samp] >= min_depth:
            is_pass = 1
      return is_pass
      
   def filter_min_depth_geneLevel(self,gene,min_depth):
      """
      Return a bool value for cirRNA gene for whether there is a sample have more than (min_depth) reads in each position for this given gene. 
      """
      is_pass = 0
      M_gene_filter = {}
      for samp in self.gene_count_info[gene]:
         M_gene_filter[samp] = []
         for cir_cnt in self.gene_count_info[gene][samp]:
            if cir_cnt >= min_depth:
               M_gene_filter[samp].append( cir_cnt )
               is_pass = 1
      return is_pass,M_gene_filter
      
         
      
      
      
      