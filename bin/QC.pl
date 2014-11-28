#!/usr/bin/perl

use strict;
use Getopt::Long;
use PerlIO::gzip;

my $usage = <<__EOF__;
USAGE
   perl $0 <indir> <outdir> <end(1|2)> [options]
OPTION
   adapt   Adapter sequence list (default: /gpfsdata/Analysis/zhuping/bin/adapter.list).
   ascii   ASCII number that stands for quality 0 (default: 33).
   Nrate   N base rate to discard read (default: 0.1).
   Qmin    Minimal quality be regarded as too low (default: <ascii> + 5).
   Qrate   Low quality base rate (default: 0.5).
   main    Identifier in main title of graphs (default:<outdir>).
__EOF__

die $usage if @ARGV<3;
my $indir = shift @ARGV;
my $outdir = shift @ARGV;
chop($outdir) if $outdir =~ /\/$/;
my $end = shift @ARGV;
my ($adaptlist,$asc,$Nrate,$Qmin,$Qrate,$mainID);
GetOptions(
   "adapt=s"=>\$adaptlist,
   "ascii=i"=>\$asc,
   "Nrate=i"=>\$Nrate,
   "Qmin=i"=>\$Qmin,
   "Qrate=i"=>\$Qrate,
   "main=s"=>\$mainID,
);

$adaptlist ||= "/data/Analysis/huboqiang/database/adapter.list";
$asc ||= 33;
$Nrate ||= 0.1;
$Qmin ||= ($asc + 5);
$Qrate ||= 0.5;
$end ||= 2;
$mainID ||= $outdir;

system("mkdir -p $outdir") unless -d $outdir;
open LOG,">$outdir/log" or die $!;
open A,"<$adaptlist" or die $!; 
chomp(my @adapter = <A>);
close A;

my ($total_reads,$total_bases,$remanent_reads,$remanent_bases,$without_adapter_reads,$read_length,$adapter_num,$remove_N_num,$low_quality_num) = (0) x 9;
my (%hash_base,%hash_quality);
my ($gc_1,$Q20_1,$Q30_1,$error_1) = (0) x 4;
my ($gc_2,$Q20_2,$Q30_2,$error_2,$remove_duplication_num) = (0) x 5;
my %duplication;
my ($Q20,$Q30) = ($asc+20,$asc+30);

if($end == 2){
   my ($file_1,$file_2) = glob("$indir/*.gz");
   open IN_1,"gzip -dc $file_1|" or die $!;
   open IN_2,"gzip -dc $file_2|" or die $!;
   open OUT_1,">:gzip", "$outdir/1.cln.fq.gz" or die $!;
   open OUT_2,">:gzip", "$outdir/2.cln.fq.gz" or die $!;
	# Read the two raw-data files (<:gzip)
	# QC
	# Print the clean-data which passed QC(<:gzip)
	# QC will :
	# 
	#	1:	remove_adapter										sub remove_adapter
	#	2:	too much N base in Read, discard.			sub count_bases
	#	3:	too much low-quality in Read, discard.    sub count_quality
	
   while(1){
      my $line1_1 = <IN_1>;
      my $line1_2 = <IN_1>;
      my $line1_3 = <IN_1>;
      my $line1_4 = <IN_1>;
      my $line2_1 = <IN_2>;
      my $line2_2 = <IN_2>;
      my $line2_3 = <IN_2>;
      my $line2_4 = <IN_2>;
		print "$total_reads\t$line1_1\t$line2_1\n";
      last unless (defined($line1_1) && defined($line2_1));
      chomp($line1_1,$line1_2,$line1_3,$line1_4,$line2_1,$line2_2,$line2_3,$line2_4);
      ++$total_reads;
      $read_length = length($line1_2) if ($read_length == 0);
      my $remove_a1 = &remove_adapter($line1_2);
      my $remove_a2 = &remove_adapter($line2_2);
      if($remove_a1 or $remove_a2){
         $adapter_num++;
         next;
      }
      my $remove_n1 = &count_bases($line1_2,0);
      my $remove_n2 = &count_bases($line2_2,$read_length);
      ++$remove_N_num if($remove_n1 or $remove_n2);
      my $low1 = &count_quality($line1_4,\$Q20_1,\$Q30_1,0);
      my $low2 = &count_quality($line2_4,\$Q20_2,\$Q30_2,$read_length);
      ++$low_quality_num if($low1 or $low2);
      next if ($remove_n1 or $remove_n2 or $low1 or $low2);
      ++$remanent_reads;
      print OUT_1 "$line1_1\n$line1_2\n$line1_3\n$line1_4\n";
      print OUT_2 "$line2_1\n$line2_2\n$line2_3\n$line2_4\n";
   }
   $without_adapter_reads = $total_reads - $adapter_num;
   close IN_1;close IN_2;close OUT_1;close OUT_2;
   &error_rate(2);
   &gc_content(2);
   &caculate_rates(2);
}else{
   my $file = glob("$indir/*.gz");
   open IN_1,"<:gzip", $file or die $!;
   open OUT_1,">:gzip", "$outdir/1.cln.fq.gz" or die $!;
   while (1) {
      my $line1 = <IN_1>;
      my $line2 = <IN_1>;
      my $line3 = <IN_1>;
      my $line4 = <IN_1>;
      last unless (defined($line1));
      chomp ($line1,$line2,$line3,$line4);
      ++$total_reads;
      ($read_length = length($line2)) if ($read_length == 0);
      if(&remove_adapter($line2)){
         ++$adapter_num;
         next;
      }
      my $remove_n1 = &count_bases($line2,0);
      ++$remove_N_num if($remove_n1);
      my $low1 = &count_quality($line4,\$Q20_1,\$Q30_1,0);
      ++$low_quality_num if($low1);
      next if ($remove_n1 or $low1);
      ++$remanent_reads;
      print OUT_1 "$line1\n$line2\n$line3\n$line4\n";
   }
   $without_adapter_reads = $total_reads - $adapter_num;
   close IN_1;
   close OUT_1;
   &error_rate(1);
   &gc_content(1);
   &caculate_rates(1);
}

my $X_axis;
($end == 2) ? ($X_axis = $read_length * 2) : ($X_axis = $read_length);
my $vertical_bar;
($end == 2) ? ($vertical_bar = "abline(v=$read_length,col='darkblue',lty=2)") : ($vertical_bar = "");

my $GC_figure = <<__EOF__;
gc<-read.table("$outdir/NC")
site<-gc[,1]
base_a<-gc[,4]
base_t<-gc[,7]
base_g<-gc[,10]
base_c<-gc[,13]
base_n<-gc[,16]
total_sites<-$X_axis
half_sites<-$read_length/2
pdf("$outdir/NC.pdf",width=8,height=6)
plot(site,base_a,xlim=c(0,total_sites),ylim=c(0,50),axes=FALSE,col="red",type="l",xlab="Position along reads",ylab="percent",main="Base percentage composition of\\n$mainID",lty=1,lwd=1.5)
lines(site,base_t,col="magenta",type="l",lty=2,lwd=1.5)
lines(site,base_g,col="darkblue",type="l",lty=4,lwd=1.5)
lines(site,base_c,col="green",type="l",lty=5,lwd=1.5)
lines(site,base_n,col="cyan3",type="l",lty=6,lwd=1.5)
legend("topright",legend=c("A","T","G","C","N"),col=c("red","magenta","darkblue","green","cyan3"),lty=c(1,2,4,5,6))
$vertical_bar
axis(side=1,at=seq(from=0,to=total_sites,by=half_sites))
axis(side=2,at=seq(from=0,to=50,by=10))
dev.off()
__EOF__

my $meanQ_errorR = <<__EOF__;
table<-read.table("$outdir/BQ")
site<-table[,1]
quality<-table[,2]
error<-table[,3]
total_sites<-$X_axis
pdf("$outdir/BQ.pdf",width=8,height=6)
plot(site,quality,xlim=c(0,total_sites),ylim=c(0,40),axes=FALSE,col="red",type="p",pch=".",cex=1.5,xlab="Position along reads",ylab="Quality",main="Base quality distribution of\\n$mainID")
axis(side=1,at=seq(from=0,to=total_sites,by=20))
axis(side=2,at=seq(from=0,to=40,by=10))
abline(h=20,col="darkblue",lty=2)
abline(v=seq(0,total_sites, by=10),col="darkblue",lty=3 )
pdf("$outdir/ER.pdf",width=8,height=6)
plot(site,error,xlim=c(0,total_sites),col="red",type="h",xlab="Position along reads",ylab="% Error-Rate", main="Error rate distribution of\\n$mainID")
axis(side=1,at=seq(from=0,to=total_sites,by=20))
abline(v=seq(0,total_sites, by=10),col="darkblue",lty=3 )
dev.off()
__EOF__

open R,"|R --vanilla --slave >/dev/null" or die $!;
print R $GC_figure;
print R $meanQ_errorR;
close R;

sub remove_adapter{ # 0-sequence
   my $seq = $_[0];
   for(my $i=0;$i<@adapter;++$i){
      return 1 if ($seq =~ /$adapter[$i]/);
   }   
   return 0;
}
sub count_bases{
   my ($seq,$start_site) = @_;
   my $length = length($seq);
   for my $i(0 .. $length-1) {
      my $base = substr($seq,$i,1);
      ++$hash_base{$i+$start_site}{$base};
   }
   my $N_num = ($seq =~ tr/N/N/);
   ($N_num >= $length * $Nrate)? return 1: return 0;
}
sub count_quality{ # 0-quality line, 1-q20, 2-q30, 3-start
   my ($seq,$Q_20,$Q_30,$start_site) = @_;
   my ($low_q_site,$base_quality);
   my $length = length($seq);
   for my $i(0 .. $length-1){
      my $base_asc = substr($seq,$i,1);
      $base_quality = ord($base_asc);
      $hash_quality{$i+$start_site} += $base_quality;
      ++$low_q_site if ($base_quality <= $Qmin);
      ++$$Q_20 if ($base_quality >= $Q20);
      ++$$Q_30 if ($base_quality >= $Q30);
   }
   ($low_q_site >= $length*$Qrate) ? (return 1) : (return 0);
}
sub error_rate{ # 0-end
   my $end = $_[0];
   open OUT_3,">$outdir/BQ" or die $!;
   my @keys = sort {$a<=>$b} keys %hash_quality;
   my $minus = $asc;
   for(my $i=0;$i<@keys;++$i) {
      my $mean_quality = ($hash_quality{$keys[$i]}/$without_adapter_reads) - $minus;
      my $index = 0 -($mean_quality/10);
      my $error_rate = (10**$index)*100;
      ($i < $read_length)? $error_1 += $error_rate: $error_2 += $error_rate;
      printf OUT_3 "%d\t%.5f\t%f\n",$keys[$i],$mean_quality,$error_rate;
   }
   close OUT_3;
}
sub gc_content{ # 0-end
   my $end = $_[0];
   open OUT_4,">$outdir/NC" or die $!;
   my @keys = sort {$a<=>$b} keys %hash_base;
   my @bases = qw/A T G C N/;
   for(my $i=0;$i<@keys;++$i) {
      print OUT_4 "$keys[$i]\t";
      for(my $j=0;$j<@bases;++$j) {
         if(exists $hash_base{$keys[$i]}{$bases[$j]}){
            my $frequency = ($hash_base{$keys[$i]}{$bases[$j]}/$without_adapter_reads)*100;
            printf OUT_4 "%s\t%d\t%.3f\t",$bases[$j],$hash_base{$keys[$i]}{$bases[$j]},$frequency;
         }else{
            print OUT_4 "$bases[$j]\t0\t0\t";
         }
      }
      print OUT_4 "\n";
      my ($g,$c) = (0) x 2;
      $g = $hash_base{$keys[$i]}{"G"} if exists($hash_base{$keys[$i]}{"G"});
      $c = $hash_base{$keys[$i]}{"C"} if exists($hash_base{$keys[$i]}{"C"}); 
      ($i < $read_length) ? ($gc_1 += $g + $c) : ($gc_2 += $g + $c);
   }
   close OUT_4;
}
sub caculate_rates{ # 0-end
   my $end = $_[0];
   my ($gc_rate_2,$Q20_rate_2,$Q30_rate_2,$error_rate_2,$duplication_rate);
   $remanent_bases = $remanent_reads * $read_length;
   $total_bases = $total_reads * $read_length;
   my $without_adapter_bases = $without_adapter_reads * $read_length;
   my $gc_rate_1 = ($gc_1/$without_adapter_bases)*100;
   my $Q20_rate_1 = ($Q20_1/$without_adapter_bases)*100;
   my $Q30_rate_1 = ($Q30_1/$without_adapter_bases)*100;
   my $error_rate_1 = $error_1/100;
   if ($end == 2) {
      $total_reads = $total_reads * 2;
      $total_bases = $total_reads * $read_length;
      $remanent_reads = $remanent_reads * 2;
      $remanent_bases = $remanent_reads * $read_length;
      $gc_rate_2 = ($gc_2/$without_adapter_bases)*100;      
      $Q20_rate_2 = ($Q20_2/$without_adapter_bases)*100;   
      $Q30_rate_2 = ($Q30_2/$without_adapter_bases)*100;      
      $error_rate_2 = $error_2/100;
   }
   
   my $title = "Raw reads\tRaw bases\tClean reads\tClean bases\tErrorRate\tQ20\tQ30\tGC content\n";
   printf LOG $title;
   my $output1 = "$total_reads\t$total_bases\t$remanent_reads\t$remanent_bases\t";
   printf LOG $output1;
   if ($end == 2){
      printf LOG "%.2f;%.2f\t%.2f;%.2f\t%.2f;%.2f\t%.2f;%.2f\n",$error_rate_1,$error_rate_2,$Q20_rate_1,$Q20_rate_2,$Q30_rate_1,$Q30_rate_2,$gc_rate_1,$gc_rate_2;
   }else{
      printf LOG "%.2f\t%.2f\t%.2f\t%.2f\n",$error_rate_1,$Q20_rate_1,$Q30_rate_1,$gc_rate_1;
   }
   my $output2 = "N remove $remove_N_num\nQuality remove $low_quality_num\nAdapter remove $adapter_num\n";
   printf LOG $output2;
   close LOG;
}
