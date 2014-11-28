from __future__ import division
import re,sys,os
import module_StatInfo as Stat
import scipy
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import module_CountInfo  as m_cnt
import module_GTFFeature   as m_gtf
import module_CIRCFeature  as m_cir
import module_Matrix       as m_mat

class SampStat(dict):
   def __init__(self,samp_info,ercc_info,genome_gtf,cir_min_read,dir_name):
      self['sam_info'] = {}
      self['sample'] = []
      self['infile'] = {'info_file':samp_info,'ercc_info':ercc_info,'anno_file':genome_gtf}
      self['stage'] = { 'name':[] }
      self['dir_name'] = dir_name
      self.ERCC_info   = { 'len':{}, 'mol':{} }
      self.min_depth   = cir_min_read
      self.__load_name()
      self.__load_ERCC_info()
      
   def __load_ERCC_info(self):
      f_file = open( self['infile']['ercc_info'],"r")
      f_file.readline()
      for line in f_file:
         line = line.strip('\n')
         f    = line.split()
         ERCC_id = f[0]
         self.ERCC_info['len'][ERCC_id] = int(f[1])
         self.ERCC_info['mol'][ERCC_id] = float(f[2])
      f_file.close()
      
   def __load_name(self):
      self.cln          = self['dir_name']['clean_data']
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
      self.statInfo     = "%s/StatInfo"% (home_dir)
      if not os.path.isdir( self.statInfo ):
         os.mkdir( self.statInfo )

   def load_samp(self):      
      self['sam_info']['samp_brief'] = {}
      self['sam_info']['brief_samp'] = {}
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
         self['sam_info']['brief_samp'][ brief_name ] = samp
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
   
   def __get_refseq_reads(self, Refseq_info, samp):
      idx = self['sample'].index( samp )
      return Refseq_info.sam_tot_reads[idx]
      
      
   def Basic_Stat(self):
      """
      Stat for QC, Tophat mapping, ERCC RGC count
      """
      out_file   = "%s/01.BasicInfo_QC_map_SpikeIn.xls" % (self.statInfo)
      f_out_file = open( out_file,"w" )
      out_info   = "Sample\tBrief_samp\tRaw_Reads\tClean_Reads\t"    +\
                   "Pre_Map_Reads\tAligned_Reads\tRefseq_Reads\t"    +\
                   "Circular_genes\tCircular_exons\tCircular_reads\t"+\
                   "RFP_Reads\tGFP_Reads\tCRE_Reads\tERCC_Reads\t"   +\
                   "RFP_Mols\tGFP_Mols\tCRE_Mols\tERCC_Mols"
                   
      print >>f_out_file, out_info

      l_breif_samp = [ self['sam_info']['samp_brief'][samp] for samp in self['sample'] ]

      '''
         Load refseq reads
      '''
      Refseq_info = m_cnt.CountInfo( self.HTS_k,l_breif_samp,"dexseq_clean",self.HTS )
      Refseq_info.load_mat()
      Refseq_info.sam_tot_reads()
      
      '''
         Load circular reads. Run this step after CIRC_Stat.
      '''
      exon_level_MinDepth = "%s/04.CIRC_info/CIRC_PE_result.merge_exon_level.MinDepth_%d.xls" % ( self.statInfo,self.min_depth )
      gene_level_MinDepth = "%s/04.CIRC_info/CIRC_PE_result.merge_gene_level.MinDepth_%d.xls" % ( self.statInfo,self.min_depth )

      exon_level_mat = m_mat.Matrix_info( exon_level_MinDepth,2 )
      gene_level_mat = m_mat.Matrix_info( gene_level_MinDepth,1 )

      exon_level_mat.load_mat()
      gene_level_mat.load_mat()
      
      # How many genes have cirRNA with depth>2 in one junction of this gene for a given sample?
      # How many exons have cirRNA with depth>2 in one junction for a given sample?
      # Sum( reads ) for a given sample?
      np_gene_cirs_samp = np.sum( gene_level_mat.matrix >= self.min_depth ,axis=0)
      np_exon_cirs_samp = np.sum( exon_level_mat.matrix >= self.min_depth ,axis=0)
      np_exon_read_samp = np.sum( exon_level_mat.matrix                   ,axis=0)
      
      '''
         Load other information
      '''
      for idx,samp in enumerate(self['sample']):
         brief_name  = self['sam_info']['samp_brief'][samp]
         QC_log      =  "%s/%s/log"                         % ( self.cln   ,samp )
         Tophat_log  =  "%s/%s/align_summary.txt"           % ( self.tophat,brief_name )
         HTSeq_SpikeIn= "%s/%s/%s.dexseq_ERCC_RGCPloyA.txt" % ( self.HTS_k ,brief_name, brief_name )
         
         
         QcStat_info = Stat.QcStat( QC_log )
         MapStat_info= Stat.TophatStat( Tophat_log )
         SpikeIn_info= Stat.SpikeIn( HTSeq_SpikeIn, self['infile']['ercc_info'] )
         
      
         QcStat_info.read_infile()
         MapStat_info.read_infile()
         SpikeIn_info.load_HTS_file()
         
         pre_map_read = MapStat_info['statInfo']['totalRead']
         aligned_read = MapStat_info['statInfo']['mappedRead']
         refseq_read  = self.__get_refseq_reads( Refseq_info,samp )
         
         cir_genes   = np_gene_cirs_samp[idx]
         cir_exons   = np_exon_cirs_samp[idx]
         cir_reads   = np_exon_read_samp[idx]
         
         read_RFP   =  SpikeIn_info.RGC_count['RGC-mRFP']
         read_GFP   =  SpikeIn_info.RGC_count['RGC-GFP' ]
         read_CRE   =  SpikeIn_info.RGC_count['RGC-CRE' ]
         read_ERCC  =  SpikeIn_info.ERCC_total
         
         mol_RFP  = self['sam_info']['RFP_mols'][ samp ]
         mol_GFP  = self['sam_info']['GFP_mols'][ samp ]
         mol_CRE  = self['sam_info']['CRE_mols'][ samp ]
         mol_ERCC = self['sam_info']['dilute'][   samp ] * 6.023*10**10
         
         out_info =  "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%1.2e\t%1.2e\t%1.2e\t%1.2e"   \
            % ( samp, brief_name, QcStat_info.raw_reads, QcStat_info.cln_reads,  \
                pre_map_read, aligned_read, refseq_read,                         \
                cir_genes,    cir_exons   , cir_reads,                           \
                read_RFP, read_GFP,read_CRE, read_ERCC,                          \
                mol_RFP , mol_GFP ,mol_CRE , mol_ERCC )
         print  >>f_out_file, out_info
      f_out_file.close()
   
   def ERCC_Stat(self):
      """
      Using ERCC mols to estimate the total amount of refGene RNA-mols in a cell. 
      (ERCC_MOLs) / (ERCC_FPKM) = (mRNA_MOLs) / (mRNA_FPKM)
      """
      self.l_ERCC_name = []
      self.l_RGCs_name = []
      self.l_mRNA_name = []
      self.l_ERCC_FPKM = {}
      self.l_RGCs_FPKM = {}
      self.l_mRNA_FPKM = {}
      self.l_cirRNA_FPKM={}

      self.l_ERCC_HTSname = []
      self.l_RGCs_HTSname = []
      self.l_mRNA_HTSname = []
      self.l_ERCC_RPKM = {}
      self.l_RGCs_RPKM = {}
      self.l_mRNA_RPKM = {}


      self.l_ERCC_MOLs = {}
      self.l_RGCs_MOLs = {}
      self.l_mRNA_MOLs = {}
      self.l_cirRNA_MOLs={}
      self.l_mRNA_MOLs_HTSname = {}

      
      self.regression = {}
      
      self.__load_FPKM()
      self.__load_MOLs()      # ERCC RGC mols
      self.__get_mRNA_MOLs()  # get mRNA mols using ERCC_FPKM, ERCC_MOLs and mRNA_FPKM
      
      exon_level_MinDepth = "%s/04.CIRC_info/CIRC_PE_result.merge_exon_level.MinDepth_%d.xls" % ( self.statInfo,self.min_depth )
      self.__cir_mols(  exon_level_MinDepth, self.min_depth )
      
      out_file   = "%s/02.ERCC_Mols.xls" % (self.statInfo)
      f_out_file = open( out_file,"w" )
      
      out_info = "Sample\tBrief_samp\tERCC_MOLs\tRGC_MOLs\tmRNA_MOLs\tcirRNA_MOLs\tRefSeq_mRNA_MOLs\tRefSeq_mRNA_FPKM>0.1"
      print >>f_out_file, out_info
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         ERCC_MOLs  = sum(    self.l_ERCC_MOLs[  brief_name ] )
         RGC_MOLs   = sum(    self.l_RGCs_MOLs[  brief_name ] )
         mRNA_MOLs  = np.sum( self.l_mRNA_MOLs[  brief_name ] )
         cirRNA_MOLs= np.sum( self.l_cirRNA_MOLs[brief_name ] )
         RefSeq_mRNA_MOLs           = np.sum(  self.l_mRNA_MOLs[ brief_name ][ self.mRNA_refSeq_index ] )
         RefSeq_mRNA_lFPKM          = np.array(self.l_mRNA_FPKM[ brief_name ],dtype=float )
         RefSeq_mRNA_lFPKM          = RefSeq_mRNA_lFPKM[ self.mRNA_refSeq_index ] 
         RefSeq_mRNA_Exps           = np.shape(  RefSeq_mRNA_lFPKM[RefSeq_mRNA_lFPKM >= 0.1] )[0]
         
         out_info = "%s\t%s\t%6.3e\t%6.3e\t%6.3e\t%6.3e\t%6.3e\t%d" % ( samp,brief_name, ERCC_MOLs,RGC_MOLs,mRNA_MOLs,cirRNA_MOLs, RefSeq_mRNA_MOLs,RefSeq_mRNA_Exps )
         print >>f_out_file, out_info
      
      f_out_file.close()
   
   def CIRC_Stat(self):
      
      if not os.path.isdir( "%s/04.CIRC_info" % ( self.statInfo) ):
         os.mkdir( "%s/04.CIRC_info" % ( self.statInfo) )
      
      exon_level = "%s/04.CIRC_info/CIRC_PE_result.merge_exon_level.xls"   % ( self.statInfo )
      circ_infor = "%s/04.CIRC_info/merge_exon.statInfo"                   % ( self.statInfo )
      gene_level = "%s/04.CIRC_info/CIRC_PE_result.merge_gene_level.xls"   % ( self.statInfo )
      gene_cirs  = "%s/04.CIRC_info/CIRC_PE_result.merge_gene_CirCount.xls"% ( self.statInfo )
      
      l_brief_name = [  self['sam_info']['samp_brief'][samp] for samp in self['sample']  ]
      
      circ_info  = m_cir.CirFeature(self.CIRCexplorer, l_brief_name, "CIRCexplorer_circ_PE.txt")
      circ_info.load_CIRC()
      circ_info.merge_gene_count()
      
      self.__mrg_exon( circ_info,exon_level           )
      self.__mrg_gene( circ_info,gene_level,gene_cirs )
      
      exon_level_MinDepth = "%s/04.CIRC_info/CIRC_PE_result.merge_exon_level.MinDepth_%d.xls" % ( self.statInfo,self.min_depth )
      gene_level_MinDepth = "%s/04.CIRC_info/CIRC_PE_result.merge_gene_level.MinDepth_%d.xls" % ( self.statInfo,self.min_depth )
      gene_cirs_MinDepth  = "%s/04.CIRC_info/CIRC_PE_result.merge_gene_CirCount.MinDepth%d.xls"% (self.statInfo,self.min_depth )
      
      self.__mrg_exon( circ_info,exon_level_MinDepth, self.min_depth                    )
      self.__mrg_gene( circ_info,gene_level_MinDepth,gene_cirs_MinDepth, self.min_depth )
      self.__report_cirFeature( circ_info, circ_infor, self.min_depth )
      
      
      
      
   def __cir_mols(self, exon_level_MinDepth, min_depth=0):
      l_brief_name = [  self['sam_info']['samp_brief'][samp] for samp in self['sample']  ]
      ltype  = "circ_mols"
      inFile = "%s/04.CIRC_info"%(self.statInfo)
      outFile= "%s/04.CIRC_info"%(self.statInfo)
      circ_Mols = m_cnt.CountInfo( inFile, l_brief_name, ltype, outFile )
      circ_Mols.load_mat( exon_level_MinDepth,gen_col=2 )
      circ_Mols.sam_tot_reads()
      
      M_gen_len = {}
      for gene in circ_Mols.gene:
         M_gen_len[gene] = { 'max_len':200 }    #### Please ONLY USE samples with less index, that sequenced by illumina X10 !!!!!!
      
      circ_Mols.cal_RPKM( M_gen_len,self.tophat )
      out_file_RPKM    = "%s/merge.%s.RPKM.xls" % ( outFile, ltype )
      
      cirRPKM_mat = m_mat.Matrix_info( out_file_RPKM,inf_column=1,in_dtype="float" )
      cirRPKM_mat.load_mat()
      l_cir_genes = cirRPKM_mat.colname
      
      for brief_name in l_brief_name:
         idx_rowname = cirRPKM_mat.rowname.index( brief_name )
         self.l_cirRNA_FPKM[brief_name] = cirRPKM_mat.matrix[:,idx_rowname]
      
      self.__get_cirRNA_MOLs()  # get cirRNA mols using ERCC_FPKM, ERCC_MOLs and cirRNA_FPKM
      
      
      
      
      
   
   def __mrg_exon(self, circ_info, exon_level, min_depth=0):
      f_exon_level = open( exon_level,"w" )
      out = "chrom:beg-end\tgene"
      for samp in circ_info.l_samp:
         out += "\t%s" % (samp)
      print >>f_exon_level, out
   
      for gene in sorted(circ_info.gen2cir_pos):
         for pos in sorted(circ_info.gen2cir_pos[gene]):
            is_pass_pos = circ_info.filter_min_depth( pos,min_depth )
            if not is_pass_pos:
               continue
            out = "%s\t%s" % (pos,gene)
            for samp in circ_info.l_samp:
               val = 0
               if samp in circ_info.count_info[ pos ]:
                  val = circ_info.count_info[ pos ][ samp ]
                  if val < min_depth:
                     val = 0
               out += "\t%d" % (val)
            
            print >>f_exon_level, out
      f_exon_level.close()
      
   def __mrg_gene( self,circ_info, gene_level,gene_cirs, min_depth=0 ):
      f_gene_level = open( gene_level,"w" )
      f_gene_cirs  = open( gene_cirs ,"w" )
      out = "gene"
      for samp in circ_info.l_samp:
         out += "\t%s" % (samp)
      print >>f_gene_level, out
      print >>f_gene_cirs , out
      
      for gene in sorted(circ_info.gen2cir_pos):
         is_pass_pos,M_gen_sam_cnt = circ_info.filter_min_depth_geneLevel( gene,min_depth )
         if not is_pass_pos:
            continue
         out = "%s" % (gene)
         out1= "%s" % (gene)
         for samp in circ_info.l_samp:
            val = 0
            val1= 0
            if samp in circ_info.gene_count_info[gene]:
               val = sum(M_gen_sam_cnt[samp])
               val1= len(M_gen_sam_cnt[samp])
            out += "\t%d" % ( val )
            out1+= "\t%d" % ( val1)
         print >>f_gene_level, out
         print >>f_gene_cirs , out1
      
      f_gene_level.close()
      f_gene_cirs.close()
            
   def __report_cirFeature( self, circ_info, circ_infor, min_depth ):
      f_cir_infor_exon_cnt_len   = open( "%s.exon_cnt_len.xls" % (circ_infor),"w" )
      f_cir_infor_intron         = open( "%s.intron.bed"       % (circ_infor),"w" )
      f_cir_infor_intron_len     = open( "%s.intron_len.xls"   % (circ_infor),"w" )
      f_cir_infor_each_samp      = open( "%s.sample_sum.xls"   % (circ_infor),"w" )
      
      
      head_exon_cnt_len = "chr_beg_end_strand\texon_count\texon_length\tavg_exon_len"
      head_intron_len   = "chr_beg_end_strand\tintron_gene\tintron_info\tup_stream_intron_len\tdown_stream_intron_len"
      
      print >>f_cir_infor_exon_cnt_len, head_exon_cnt_len
      print >>f_cir_infor_intron_len  , head_intron_len
      
      pat      = "(\w+)\:(\w+)\-(\w+),(\S+)"
      pattern  = re.compile( pat )
      
      for pos_info in circ_info.cir_info:
         is_pass_pos = circ_info.filter_min_depth( pos_info,min_depth )
         if not is_pass_pos:
            continue
         exon_cnts = circ_info.cir_info[pos_info]['exon_cnt']
         exon_lens = ",".join( [ str(size) for size in circ_info.cir_info[pos_info]['l_exon_size'] ] )
         avg_exon_len = sum(circ_info.cir_info[pos_info]['l_exon_size'])/len(circ_info.cir_info[pos_info]['l_exon_size'])
         print >>f_cir_infor_exon_cnt_len, "%s\t%d\t%s\t%6.2f" % ( pos_info,exon_cnts, exon_lens, avg_exon_len )
         
         match = pattern.search( pos_info )
         chrom,beg,end,strand = match.group(1),match.group(2),match.group(3),match.group(4)
         intron_1,intron_2 = circ_info.cir_info[ pos_info ][ 'intron_pos' ].split('|')

         chr_pos  = "%s:%s-%s" % (chrom,beg,end) 
         beg_i1,end_i1,beg_i2,end_i2   = ["0","0","0","0"]
         
         intron_gene = circ_info.cir_info[pos_info]['gene_name']
         ltype = "left"
         if strand == "-":
            ltype = "right"
            
         l_out_info = []
         if intron_1 != "None":
            chrom_i1,beg_i1,end_i1 = re.split('[\:\-\,]',intron_1)
            l_out_info.append( "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( chrom_i1,beg_i1,end_i1,chr_pos,strand,intron_gene,ltype ) )

         ltype = "right"
         if strand == "-":
            ltype = "left"            
         out_info_2 = ""
         if intron_2 != "None":
            chrom_i2,beg_i2,end_i2 = re.split('[\:\-\,]',intron_2)
            l_out_info.append( "%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( chrom_i2,beg_i2,end_i2,chr_pos,strand,intron_gene,ltype ) )

         out_info = "\n".join( l_out_info )
         print >>f_cir_infor_intron, out_info

         len_l = int( end_i1 ) - int( beg_i1 )
         len_r = int( end_i2 ) - int( beg_i2 )
         if strand == "-":
            len_r = int( end_i1 ) - int( beg_i1 )
            len_l = int( end_i2 ) - int( beg_i2 )
         
         print >>f_cir_infor_intron_len, "%s\t%s\t%s\t%d\t%d" % ( pos_info, intron_gene, circ_info.cir_info[pos_info]['intron_pos'], len_l, len_r )
      
      f_cir_infor_exon_cnt_len.close()
      f_cir_infor_intron.close()
      f_cir_infor_intron_len.close()
      
            
   def __load_FPKM(self):
      Cuffquant_file = "%s/genes.fpkm_table" % ( self.cuffnorm )
      f_Cuffquant    = open( Cuffquant_file,"r" )
      in_head        = f_Cuffquant.readline()
      f_head         = in_head.split()
      l_samp_brief   = [ sam.split("_0")[0] for sam in f_head[1:] ]
      
      for brief_name in l_samp_brief:
         self.l_ERCC_FPKM[brief_name] = []
         self.l_RGCs_FPKM[brief_name] = []
         self.l_mRNA_FPKM[brief_name] = []
      
      for line in f_Cuffquant:
         line = line.strip('\n')
         f    = line.split()
         gene = f[0]
         
         if gene[0:3] == "MIR":
            continue
         if gene[0:5] == "SNORD":
            continue
         
         if gene[0:5] == "ERCC-":
            self.l_ERCC_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_ERCC_FPKM[brief_name].append( f[i+1] )
               
         elif gene[0:4] == "RGC-":
            self.l_RGCs_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_RGCs_FPKM[brief_name].append( f[i+1] )
         
         else:
            self.l_mRNA_name.append( gene )
            for i,brief_name in enumerate( l_samp_brief ):
               self.l_mRNA_FPKM[brief_name].append( f[i+1] )
      
      f_Cuffquant.close()
      
      self.__mRNA_refSeq_index()   
   
   def __mRNA_refSeq_index(self):
      self.mRNA_refSeq_index = []
      for i,gene in enumerate(self.l_mRNA_name):
         if gene[0:7] == "NONHSAG" or gene[0:5] == "XLOC_":
            continue
         self.mRNA_refSeq_index.append( i )
      self.mRNA_refSeq_index = np.array( self.mRNA_refSeq_index,dtype=int )
   
   def __load_MOLs(self):
      for samp in self['sample']:

         brief_name = self['sam_info']['samp_brief'][samp]
         if brief_name not in self.l_ERCC_MOLs:
            self.l_ERCC_MOLs[brief_name] = []
         if brief_name not in self.l_RGCs_MOLs:
            self.l_RGCs_MOLs[brief_name] = []
            
         for ERCC_id in self.l_ERCC_name:
            self.l_ERCC_MOLs[brief_name].append( self.ERCC_info['mol'][ERCC_id]*(6.02*10**23/10**18)*self['sam_info']['dilute'][samp] )
         
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['CRE_mols'][ samp ] )
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['GFP_mols'][ samp ] )
         self.l_RGCs_MOLs[ brief_name ].append( self['sam_info']['RFP_mols'][ samp ] )
         
   def __get_mRNA_MOLs(self):
      self.regression = {}
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         np_ERCC_FPKM = np.array( self.l_ERCC_FPKM[brief_name],dtype=float )
         np_ERCC_MOLs = np.array( self.l_ERCC_MOLs[brief_name],dtype=float )
         
         np_ERCC_FPKM = np.log10( np_ERCC_FPKM+0.1 )
         np_ERCC_MOLs = np.log10( np_ERCC_MOLs )

         mol_idx = ( np_ERCC_MOLs - np.log10(6.02*10**23/10**18) > -3 )

         np_ERCC_FPKM = np_ERCC_FPKM[ mol_idx ]
         np_ERCC_MOLs = np_ERCC_MOLs[ mol_idx ]
         
         slope,intercept,r_value,p_value,slope_std_error =  scipy.stats.linregress( np_ERCC_MOLs,np_ERCC_FPKM )
         
         self.regression[brief_name] = {'slope'    :slope,           \
                                        'inter'    :intercept,       \
                                        'r_value'  :r_value,         \
                                        'p_value'  :p_value,         \
                                        'std_err'  :slope_std_error  \
         }
         
         np_mRNA_FPKM = np.array( self.l_mRNA_FPKM[brief_name],dtype=float )
         np_mRNA_FPKM = np.log10( np_mRNA_FPKM+0.1 )
         np_mRNA_MOLs = ( np_mRNA_FPKM - intercept ) / slope
         np_mRNA_MOLs = np.power( np_mRNA_MOLs, 10 )
         np_mRNA_MOLs[ np_mRNA_MOLs<0.0            ] = 0.0
         np_mRNA_MOLs[ np_mRNA_FPKM<=np.log10(0.1) ] = 0.0
         
         self.l_mRNA_MOLs[ brief_name ] = np_mRNA_MOLs
   
   def __get_cirRNA_MOLs(self):
      self.regression = {}
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         np_ERCC_FPKM = np.array( self.l_ERCC_FPKM[brief_name],dtype=float )
         np_ERCC_MOLs = np.array( self.l_ERCC_MOLs[brief_name],dtype=float )
         
         np_ERCC_FPKM = np.log10( np_ERCC_FPKM+0.1 )
         np_ERCC_MOLs = np.log10( np_ERCC_MOLs )

         mol_idx = ( np_ERCC_MOLs - np.log10(6.02*10**23/10**18) > -3 )
         
         np_ERCC_FPKM = np_ERCC_FPKM[ mol_idx ]
         np_ERCC_MOLs = np_ERCC_MOLs[ mol_idx ]
         
         slope,intercept,r_value,p_value,slope_std_error =  scipy.stats.linregress( np_ERCC_MOLs,np_ERCC_FPKM )
         
         self.regression[brief_name] = {'slope'    :slope,           \
                                        'inter'    :intercept,       \
                                        'r_value'  :r_value,         \
                                        'p_value'  :p_value,         \
                                        'std_err'  :slope_std_error  \
         }
         
         np_cirRNA_FPKM = np.array( self.l_cirRNA_FPKM[brief_name],dtype=float )
         np_cirRNA_FPKM = np.log10( np_cirRNA_FPKM+0.1 )
         np_cirRNA_MOLs = ( np_cirRNA_FPKM - intercept ) / slope
         np_cirRNA_MOLs = np.power( np_cirRNA_MOLs, 10 )
         np_cirRNA_MOLs[ np_cirRNA_MOLs<0.0            ] = 0.0
         np_cirRNA_MOLs[ np_cirRNA_FPKM<=np.log10(0.1) ] = 0.0
         
         self.l_cirRNA_MOLs[ brief_name ] = np_cirRNA_MOLs
   
   
   
   
   def plot_regression(self):
      pdfname = "All_samples.rpkm_cnt.pdf" 
      fig = plt.figure(figsize=(72,72))
      cnt = 0
      for samp in self['sample']:
         brief_name = self['sam_info']['samp_brief'][samp]
         
         slope          =  self.regression[brief_name]['slope']
         intercept      =  self.regression[brief_name]['inter']
         r_value        =  self.regression[brief_name]['r_value']
         p_value        =  self.regression[brief_name]['p_value']
         slope_std_error=  self.regression[brief_name]['std_err']
         
         np_ERCC_FPKM = np.array( self.l_ERCC_FPKM[brief_name],dtype=float )
         np_ERCC_MOLs = np.array( self.l_ERCC_MOLs[brief_name],dtype=float )/(6.02*10**23/10**18)
      
         np_ERCC_FPKM = np.log10( np_ERCC_FPKM+0.1 )
         np_ERCC_MOLs = np.log10( np_ERCC_MOLs )

         mol_idx = ( np_ERCC_MOLs > -3 )

         np_ERCC_FPKM_use = np_ERCC_FPKM[ mol_idx ]
         np_ERCC_MOLs_use = np_ERCC_MOLs[ mol_idx ]
         
         MOLs_predict = np.linspace( -4,2,100 )
         FPKM_predict = intercept + slope * (MOLs_predict + np.log10(6.02*10**23/10**18) )
      
         cnt += 1
         ax = plt.subplot( 8,8,cnt)
         ax.plot( np_ERCC_MOLs,     np_ERCC_FPKM      ,".r")
         ax.plot( np_ERCC_MOLs_use, np_ERCC_FPKM_use  ,"ro")
         ax.plot( MOLs_predict,     FPKM_predict           )
         title = brief_name + ' ERCC-mols vs RPKM in samples'
         ax.set_title(title,fontsize=12)
         ax.set_xlabel('Spike-in mols per reaction(attlmole,log10)')
         ax.set_ylabel('Mean read-coverage of spike-in molecules(RPKM,log10)')
         ax.text(-5,5,r'$ y = %f + %fx $' % (intercept,slope))
         ax.text(-5,4.6,r'$ Pearson R = %f $' % (r_value))
         ax.text(-5,4.2,r'$ p = %6.2e $' % (p_value))
         ax.text(-5,3.8,r'$ mRNA genes = %6.2e $' % (  np.sum( self.l_mRNA_MOLs[ brief_name ][ self.mRNA_refSeq_index ] )  ) )
         for tick in ax.xaxis.get_major_ticks():
            tick.tick1On = True
            tick.tick2On = False
         for tick in ax.yaxis.get_major_ticks():
            tick.tick1On = True
            tick.tick2On = False
         ax.set_xlim(-6,2)
         ax.set_ylim(-2,7)
         
      plt.savefig(pdfname,format='pdf')
   
