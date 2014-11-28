from __future__ import division
import re,sys,os
import subprocess,time

class running_jobs(object):
   def __init__( self, sh_file,sh_work_file ):
      self.sh_file      = sh_file
      self.sh_work_file = sh_work_file
   
   def load_sh_file(self,sh_scripts):
      f_out = open( self.sh_file,"w" )
      print >>f_out, sh_scripts
      f_out.close()
   
   def load_sh_work_file(self,sh_scripts):
      f_out = open( self.sh_work_file,"w" )
      print >>f_out, sh_scripts
      f_out.close()
   
   def running_SGE(self,vf,maxjob=100):
      shell_work = 'perl /data/Analysis/huboqiang/bin/qsub-sge.pl --resource vf=%s  --maxjob %d %s ' % (vf, maxjob,  self.sh_work_file)
      p = subprocess.Popen(shell_work,shell='True')
      while 1:
         run_cnt = 0
         if p.poll() is None:
            run_cnt += 1
            time.sleep(3)
         if run_cnt == 0:
            break
      
