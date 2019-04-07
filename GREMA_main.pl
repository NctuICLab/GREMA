#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
my $expression;
my $knowledge;
my $threads;
my $model;
my $output_dir;
my $help;
my $program = abs_path($0);
#======================
my %gene_index;

#my $total_point = 0;
#my %confidence_level;

sub Usage {
	print STDERR "Usage: perl $program [Options]
Options:
	-i	[FILE]	time-series profile
	-o	[PATH]	Output Directory
	-k	[FILE]	knowledge of regulatory
	-m	[model] {s-system,HFODE}
	-t	[No]	Number of threads
	-h	Show the usage
	";
}
sub read_expression {
	my ($profile) = @_;
	my ($line,@ele);
	my $profile_tag = 0;
	my $total_rep = 0;
	my $total_point = 0;
	my %profile_expression;
	open PROFILE,"<",$profile;
	while($line=<PROFILE>){
		chomp $line;
		if($line =~ /repeat_number=/){
			@ele = split(/=/,$line);
			$total_rep = $ele[1];
			print STDERR "total replication:".$total_rep."\n";
		}elsif($line =~ /^rep/){
			@ele = split(/\t/,$line);
			$total_point = scalar @ele - 2;
			#print STDERR "total point:".$total_point."\n";
			my $rep_no = $ele[0];
			my $gene_name = $ele[1];
			my $profile_info;
			$rep_no =~ s/rep//;
			#print STDERR "rep".$rep_no."\n";
			if($rep_no > $total_rep){
				print STDERR "Error, the replication number is larger total replication\n";die;
			}
			for(my $i=2;$i< scalar @ele; $i++){#profile from column2
				$profile_info .= $ele[$i]."\t";
			}
			$profile_info =~ s/\t$//;#remove the last tab
			$profile_expression{$rep_no}{$gene_name} = $profile_info;
		}
	}
	close PROFILE;
	my @all_gene = (sort (keys $profile_expression{1}));#get the gene name
	my $gene_no = 0;
	foreach my $i (@all_gene){
		print STDERR $i."\n";
		$gene_index{$gene_no} = $i;
		$gene_no++;
	}
	return ($gene_no,$total_rep,$total_point,%profile_expression);
}

sub read_knowledge {
	my ($know,$total_gene_no,$know_init_file) = @_;
	my ($line,@ele);
	my $format_tag = 0;
	my %knowledge;
	open KNOW,"<",$know;
	while($line=<KNOW>){
		chomp $line;
		#print STDERR $line."\n";
		if($line eq "TF\tGENE\tREGULATORY"){
			print STDERR "format is correct\n";
			$format_tag = 1;
			next;
		}
		if($format_tag){
			@ele = split(/\t/,$line);
			my ($tf,$gene,$regulation) = ($ele[0],$ele[1],$ele[2]);
			$knowledge{$tf}{$gene} = $regulation;
			#print STDERR "TF:".$tf."\tGene:".$gene."\treg:".$regulation."\n";
		}
	}
	close KNOW;
	open INIT,">",$know_init_file;
	for(my $i=0;$i<$total_gene_no;$i++){#each target gene
		my $gene_name = $gene_index{$i};
		#print STDERR $gene_name.": ";
		for(my $j=0;$j<$total_gene_no;$j++){#each TF
			my $tf_name = $gene_index{$j};
			#print STDERR "TF:".$tf_name."\n";
			if($knowledge{$tf_name}{$gene_name} eq '0'){
				print INIT "0 ";
			}elsif(!$knowledge{$tf_name}{$gene_name}){
				print INIT "? "
			}else{
				print INIT $knowledge{$tf_name}{$gene_name}." ";
			}
		}
		print INIT "\n";
	}
	close INIT;
}
sub generate_config {
	my ($configure,$total_gene_no,$total_rep,$total_point,%profile_expression) = @_;
	open CONF,">",$configure;
	print CONF "num_var ".$total_gene_no."\n";
	print CONF "num_tf ".$total_gene_no."\n";
	print CONF "num_run ".$total_rep."\n";
	print CONF "num_trial ".$total_point."\n";
	print CONF "nochange 1000\n";
	print CONF "delta_t 1\n";
	print CONF "b_up 1\n";
	print CONF "b_down 0\n";
	print CONF "transMax_up 8\n";
	print CONF "transMax_down 0\n";
	print CONF "n_up 4\n";
	print CONF "n_down -4\n";
	print CONF "k_up 4\n";
	print CONF "k_down 0\n";
	print CONF "degrade_up 1\n";
	print CONF "degrade_down 0\n";
	print CONF "list_data\n";
	for(my $i=1;$i<=$total_rep;$i++){
		print CONF $i."\n";
		for(my $j=0;$j<$total_gene_no;$j++){
			my $gene_name = $gene_index{$j};
			print CONF $profile_expression{$i}{$gene_name}."\n";
		}
		#print CONF "\n";
	}

	close CONF;
}
sub run_quantification {
	my ($gen,$use_know,$use_config,$total_gene_no,%confidence_level) = @_;
	my $n_start;
	my $n_end;
	my @threads;
	my $thr_count = 0;
	my $gene_name;
	if($model eq "HFODE"){
		$n_start = 2;
		$n_end = 2+$total_gene_no;
	}
	my $final_results = $output_dir."/final_results.txt";
	open FINAL,">",$final_results;
	print FINAL "TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\n";
	for(my $i=0; $i<$total_gene_no; $i++){
		$threads[$thr_count] = threads->new(\&run_iga, $i,$gen,$use_know,$use_config);
		sleep(1) while(scalar threads->list(threads::running) >= $threads);
		$thr_count++;
	}
	foreach my $i (@threads){
		$i->join();
	}
	for(my $i=0; $i<$total_gene_no; $i++){
		$gene_name = $gene_index{$i};
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results(1,$iga_results);#select top 1
		print STDERR "Gene:".$gene_name." ".$select_iga_results[0]."\n";
		my @iga = split(/\t/,$select_iga_results[0]);
		for(my $r=$n_start; $r<$n_end; $r++){
			my $tf_no = $r-2;
			my $tf_name = $gene_index{$tf_no};
			my $role;
			if($iga[$r]>0){
				$role = "+";
			}elsif($iga[$r]<0){
				$role = "-";
			}else{
				$role = 0;
			}
			my $relationship = $tf_name."_".$gene_name;
			print FINAL $tf_name."\t".$gene_name."\t".$role."\t".$confidence_level{$relationship}."\n";
		}
	}
	close FINAL;
}
sub run_EMA {
	my ($gen,$type,$use_know,$use_config,$total_gene_no,$fix_ref,%confidence_level) = @_;
	my @fix_value = @{$fix_ref};
	my @threads;
	my $thr_count = 0;
	my $gene_name;
	my $select_top_no = 5;
	my $cutoff = 0.8;
	my $next_generation = $gen + 1;
	my $new_knowledge = $knowledge."_knowledge_ForStep".$next_generation;
	open KNOW,">",$new_knowledge;
	print STDERR "Step2:GRN decomposition\n";
	print STDERR "Step3:Parallel solving\n";
	for(my $i=0; $i<$total_gene_no; $i++){
		if($fix_value[$i]){
			$gene_name = $gene_index{$i};
			print STDERR $gene_name ." is finished\n";
		}else{
			$threads[$thr_count] = threads->new(\&run_iga, $i,$gen,$use_know,$use_config);
			sleep(1) while(scalar threads->list(threads::running) >= $threads);
			$thr_count++;
		}
	}
	foreach my $i (@threads){
		$i->join();
	}
	for(my $i=0; $i<$total_gene_no; $i++){
		$gene_name = $gene_index{$i};
		my (@regulatory_p,@regulatory_n,@regulatory_z);
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results($select_top_no,$iga_results);#select top 5
		for(my $j=0; $j<$total_gene_no; $j++){
			$regulatory_p[$j] = 0;
			$regulatory_n[$j] = 0;
			$regulatory_z[$j] = 0;
		}
		if($model eq "HFODE"){
			#b(i),TransMax(i),n(i,j),k(i,j),deg(i)
			my $n_start = 2;
			my $n_end = 2+$total_gene_no;
			for(my $j=0; $j<$select_top_no; $j++){
				my @iga = split(/\t/,$select_iga_results[$j]);
				for(my $r=$n_start; $r<$n_end; $r++){
					my $gene_no = $r-2;
					if($iga[$r] > 0){
						$regulatory_p[$gene_no]++;
					}elsif($iga[$r] < 0){
						$regulatory_n[$gene_no]++;
					}else{
						$regulatory_z[$gene_no]++;
					}
				}
			}
		}else{
			#S-system
			my $n_start = 1;
			my $n_end = 1+$total_gene_no;
			for(my $j=0; $j<$select_top_no; $j++){
				my @iga = split(/\t/,$select_iga_results[$j]);
				for(my $r=$n_start; $r<$n_end; $r++){
					my $gene_no = $r-1;
					if($iga[$r] > 0){
						$regulatory_p[$gene_no]++;
					}elsif($iga[$r] < 0){
						$regulatory_n[$gene_no]++;
					}else{
						$regulatory_z[$gene_no]++;
					}
				}
			}
		}#S-system model
		my $know_info;
		print STDERR "Step4:Regulation determination\n";
		for(my $j=0; $j<$total_gene_no; $j++){
			my $tf_name = $gene_index{$j};
			my $each_gene_knowledge;
			my $relationship = $tf_name."_".$gene_name;
			if($regulatory_p[$j]){
				if(($regulatory_p[$j]/$select_top_no) >= $cutoff){
					$each_gene_knowledge = "+ ";
					if(!$confidence_level{$relationship}){
						$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
					}
				}
			}
			if($regulatory_n[$j]){
				if(($regulatory_n[$j]/$select_top_no) >= $cutoff){
					$each_gene_knowledge = "- ";
					if(!$confidence_level{$relationship}){
						$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
					}
				}
			}
			if($regulatory_z[$j]){
				if(($regulatory_z[$j]/$select_top_no) >= $cutoff){
					$each_gene_knowledge = "0 ";
					if(!$confidence_level{$relationship}){
						$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
					}
				}
			}
			if(!$each_gene_knowledge){
				$each_gene_knowledge = "? ";
			}
			print STDERR "TF-gene:".$tf_name."-".$gene_name."\n";
			print STDERR "regulation:".$each_gene_knowledge."\n";
			#print STDERR "P:".$regulatory_p[$j]."\n";
			#print STDERR "N:".$regulatory_n[$j]."\n";
			#print STDERR "Z:".$regulatory_z[$j]."\n";
			print STDERR "=============================\n";
			$know_info .= $each_gene_knowledge;
		}
		$know_info =~ s/ $/\n/;
		print KNOW $know_info;	
	}
	return %confidence_level;
}
sub run_iga {
	my ($gene_no,$gen,$know,$conf) = @_;
	my $src_dir = dirname($program);
	if($model eq "HFODE"){
		my $ema_HFODE = $src_dir."/EMA_HFODE/EMA_HFODE";
		if(!-e $ema_HFODE){
			print STDERR $ema_HFODE." does not exist\n";die;
		}
		if(!-e $conf){
			print STDERR $conf." does not exist\n";die;
		}
		my $gene_name = $gene_index{$gene_no};
		my $command = $ema_HFODE." -i ".$gene_no." -n 30 -m 0 -G 10 -I ".$gen." -F 2 ".$conf." ".$know." > ".$output_dir."/".$gene_name.".txt";
		print STDERR $command."\n";
		`$command`;
	}

}
sub select_results {
	my ($top_no,$results) = @_;
	my (@unsort_fitness,@sort_fitness,%fitness_value);
	my @select_results;
	my ($line,@ele);
	my $no = 0;
	open IGA,"<",$results;
	while($line=<IGA>){
		chomp $line;
		@ele = split(/ /,$line);
		my $all_parameter_no = scalar @ele;
		my $iga_info;
		for(my $j=0; $j< ($all_parameter_no - 2); $j++){
			$iga_info .= $ele[$j]."\t";
		}
		$iga_info =~ s/\t$//;
		my $fitness_value = $ele[$all_parameter_no - 2];
		#print STDERR "fitness:".$fitness_value."\n";
		$unsort_fitness[$no] = $fitness_value;
		$fitness_value{$fitness_value} = $iga_info;
		$no++;
	}
	close IGA;
	@sort_fitness = sort {$a <=> $b} @unsort_fitness;
	for(my $i=0; $i<$top_no; $i++){
		my $fitness_value = $sort_fitness[$i];
		print STDERR $i." fitness:".$fitness_value."\t".$fitness_value{$fitness_value}."\n";
		$select_results[$i] = $fitness_value{$fitness_value};
		#print STDERR "select result:".$select_results[$i]."\n";
	}
	return @select_results;
}
sub check_knowledge {
	my ($know,$gene_no) = @_;
	my ($line,@ele);
	my $no = 0;
	open KNOW,"<",$know;
	while($line=<KNOW>){
		chomp $line;
		if($no == $gene_no){
			my $unfix_no = 0;
			@ele = split(/ /,$line);
			foreach my $i(@ele){
				if($i eq "?"){
					$unfix_no++;
				}
			}
			if($unfix_no == 0){
				return 1;
			}else{
				return 0;
			}
		}
		$no++;
	}
}
sub main {
	my ($total_gene,$total_repeat,$total_data_points,%profile) = read_expression($expression);
	print STDERR "Number of gene:".$total_gene."\n";
	my $know_init = $knowledge."_knowledge_init";
	print STDERR "know_init:".$know_init;
	read_knowledge($knowledge,$total_gene,$know_init);
	if(!-e $know_init){
		print STDERR "generate initial knowledge failed <".$know_init.">\n";die;
	}
	my $config = $expression."_config";
	generate_config($config,$total_gene,$total_repeat,$total_data_points,%profile);
	print STDERR "Step1:Initialisation\n";
	my @fix;
	my $all_fix = 0;
	my $total_fix_no = 0;
	my %fix_status;
	my %confidence;
	my $generation = 1;#first generation;
	for(my $i=0;$i<$total_gene;$i++){
		$fix[$i] = check_knowledge($know_init,$i);
		print STDERR "fix[".$i."]:".$fix[$i]."\n";
		$total_fix_no += $fix[$i];
	}
	print STDERR "Generation".$generation." fix:".$total_fix_no."\n";
	if($total_fix_no == $total_gene){
		$all_fix = 1;
	}
	
	#=======finished initial=================================================
	if(!$all_fix){
		do{
			print STDERR "Generation ".$generation."\n";
			my $use_knowledge;
			if($generation == 1){
				$use_knowledge = $know_init;
			}else{
				$use_knowledge = $knowledge."_knowledge_ForStep".$generation;
			}
			if(!-e $use_knowledge){
				print STDERR "use knowledge does not exist <".$use_knowledge.">\n";die;
			}
			%confidence = run_EMA($generation, "evolutionary",$use_knowledge,$config,$total_gene,\@fix,%confidence);
			$generation++;
			my $new_knowledge_file = $knowledge."_knowledge_Forstep".$generation;
			$total_fix_no = 0;
			print STDERR "Step5:GRN combination\n";
			for(my $i=0;$i<$total_gene;$i++){
				$fix[$i] = check_knowledge($new_knowledge_file,$i);
				print STDERR "fix[".$i."]:".$fix[$i]."\n";
				$total_fix_no += $fix[$i];
			}
			print STDERR "Generation".$generation." fix:".$total_fix_no."\n";
			print STDERR "Step6:Qualitative termination test\n";
			if($total_fix_no == $total_gene){
				$all_fix = 1;
			}else{
				print STDERR "Step7:Inheritance\n";
			}
		}while(!$all_fix)
	}
	my $final_knowledge = $knowledge."_knowledge_Forstep".$generation;
	if(!-e $final_knowledge){
		print STDERR "use final knowledge does not exist <".$final_knowledge.">\n";die;
	}
	#run_EMA($generation, "final",$final_knowledge,$config,$total_gene,\@fix);
	print STDERR "GRN quantification\n";
	run_quantification($generation,$final_knowledge,$config,$total_gene,%confidence);
	
	my $tmp_dir = $output_dir."/tmp";
	if(!-d $tmp_dir){
		my $command = "mkdir ".$tmp_dir;
		print STDERR "create tmp folder\n";
		`$command`;
	}
	my $command = "mv calprofile* ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
	$command = "mv progress* ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
	$command = "mv printdata.txt ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
}
#=========subroutine end================================
GetOptions(
	'i=s'	=>\$expression,
	'o=s'	=>\$output_dir,
	'k=s'	=>\$knowledge,
	't=i'	=>\$threads,
	'm=s'	=>\$model,
	'h'	=>\$help,
); 
if(!$expression ||!$knowledge || !$threads || $help ||$threads < 1 || !$model || !$output_dir){
	Usage();die;
}
if(!-e $expression){
	print STDERR $expression." does not exist\n";
	Usage();die;
}
if(!-e $knowledge){
	print STDERR $knowledge." does not exist\n";
	Usage();die;
}
if(($model ne "s-system")and($model ne "HFODE")){
	print STDERR "type of model error\n";
	Usage();die;
}
if(-d $output_dir){
	$output_dir =~ s/\/$//;
}else{
	print STDERR "Directory of output does not exist\n";
	Usage();die;
}

main();

