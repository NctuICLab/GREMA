#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';
my $ans = $ARGV[0];
my $outdir = $ARGV[1];
my $program = abs_path($0);
my ($line,@ele);
my (%ans);
if(@ARGV<2){
    print STDERR "perl ".$program." gold_standard output_dir\n";die;
}
my $src = dirname($program);
my $grema = $src."/GREMA_main.pl";
if(!-e $grema){
    print STDERR "GREMA <".$grema."> does not exist\n";die;
}
$outdir =~ s/\/$//;
if(!-d $outdir){
    my $command = "mkdir ".$outdir;
    print STDERR $command."\n";
    `$command`;
}
print STDERR "run GREMA\n";
my $command = "perl ".$grema." -i input/Dream4_10_1_timeseries_expression.txt -o ".$outdir." ";
$command .= "-k input/insilico_size10_1_know_knowledge.txt -m HFODE -t 2";
print STDERR $command."\n";
`$command`;
my $results = $outdir."/final_results.txt";
if(!-e $results){
    print STDERR "results<".$results."> does not exist\n";die;
}
open ANS,"<",$ans;
while($line=<ANS>){
    chomp $line;
    @ele = split(/\t/,$line);
    my $role = $ele[0]."_".$ele[1];
    $ans{$role}++;
}
close ANS;
my $tp = 0;
my $tn = 0;
my $fp = 0;
my $fn = 0;
open RES,"<",$results;
$line=<RES>;
while($line=<RES>){
    chomp $line;
    @ele = split(/\t/,$line);
    my $role = $ele[0]."_".$ele[1];
    if($ele[2] eq '0'){
        if(!$ans{$role}){
            $tn++;
        }else{
            $fn++;
        }
    }else{
        if($ans{$role}){
            $tp++;
        }else{
            $fp++;
        }
    }	
}
close RES;
print STDERR "tp:".$tp."\ntn:".$tn."\nfp:".$fp."\nfn:".$fn."\n";
my $acc = ($tp+$tn)/($tp+$tn+$fp+$fn);
my $sen = $tp/($tp+$fn);
my $spe = $tn/($tn+$fp);
print STDERR "acc:".$acc."\n";
print STDERR "spe:".$spe."\n";
print STDERR "sen:".$sen."\n";
