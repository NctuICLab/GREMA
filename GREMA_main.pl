#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Scalar::Util qw(looks_like_number);
use File::Copy qw(move copy);
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
my $expression;
my $knowledge;
my $threads;
my $model = "HFODE";
my $fitness_type = 1;
my $generation_no = 1000;
my $output_dir;
my $tmp_dir;
my $cc = 0;
my $help;
my $program = abs_path($0);
#======================
my %gene_index;
my $min = 0;
my $perturb_no = 0;	# number of extra series made from the data by adding noise (HFODE only)
my $time_list = "";	# real sampling times of the time points (HFODE only), comma separated; default: equal steps
my $time_unit = 1;	# the times are divided by this value
my $n_max = 4;	# bound of the Hill coefficients n_ij of HFODE: -n_max <= n_ij <= n_max
my $noise_degree = 0;	# standard deviation of that noise, relative to the value
my $select_top_no = 5;#number of top runs that vote for a regulation
my $h_know_file;#S-system only
my $input_list = "";	# names of the rows of the profile that are external inputs (HFODE only), comma separated
my %is_input;	# input name => 1
my $n_fit = 0;	# number of fitted genes (the targets); the variables are the fitted genes G.. followed by the inputs
my %vote_score;#TF_GENE => fitness-weighted fraction of the runs of the first iteration that use the regulation



sub Usage {
	print STDERR "Usage: perl $program [Options]
Options:
	-i	[FILE]	time-series profile
	-o	[PATH]	Output Directory
	-k	[FILE]	knowledge of regulatory
	-m	[model] {s-system,HFODE} Default is HFODE
	-f	[Fitness] {0,1,2,3,4} Default is 1
	-c	[cc]	{0,1} Default is 0, consider correlation coefficient
	-g	[Generation] Number of generation, default is 1000
	-t	[No]	Number of threads
	-min	[value] ex: 0.01 [min value] default is 0
	-perturb	[No] number of extra series made from the data by adding noise (HFODE only), default 0
	-noise	[value] ex: 0.05 standard deviation of that noise relative to the value, needed with -perturb
	-nmax	[value] HFODE only: bound of the Hill coefficients, -value <= n <= value, default is 4
	-input	[list] HFODE only: names of the rows of the profile that are external inputs (a stimulus given for every series, not fitted, ex: U1,U2); an input can regulate a gene
		(write it in the TF column of the knowledge file) but is never a target; the values must be on the scale of the genes (0.01 to 1); not with -perturb
	-time	[list] HFODE only: sampling times of the time points, comma separated and increasing, ex: 0,5,10,15,20,30,45,90; default: equal steps of 1
	-tunit	[value] the times of -time are divided by this value (the parameter ranges are made for steps of about 1), default is 1
	-h	Show the usage
	";
}
#Column range of the regulation parameters in the EMA output (one line per run)
#HFODE: b,TransMax,n(1..N),k(1..N),degrade ; S-system: alpha,g(1..N),beta,h(1..N),delta_t
sub regulation_columns {
	my ($total_gene_no) = @_;
	my $n_start = ($model eq "HFODE") ? 2 : 1;
	return ($n_start, $n_start + $total_gene_no);
}
#Natural order of gene names: G1,G2,...,G10 (not G1,G10,G2)
sub natural_cmp {
	my ($x,$y) = @_;
	my @a = split(/(\d+)/,$x);
	my @b = split(/(\d+)/,$y);
	while(@a && @b){
		my ($p,$q) = (shift @a, shift @b);
		next if($p eq $q);
		if($p =~ /^\d+$/ && $q =~ /^\d+$/){
			return $p <=> $q if($p != $q);
			next;
		}
		return $p cmp $q;
	}
	return scalar(@a) <=> scalar(@b);
}
#true if the value is a real, finite number (rejects NaN, inf, text)
sub is_finite {
	my ($x) = @_;
	return 0 if(!defined($x) || !looks_like_number($x));
	return 0 if($x != $x);#NaN
	return 0 if($x == 9**9**9 || $x == -9**9**9);#inf
	return 1;
}
#quote a string for the shell
sub shq {
	my ($s) = @_;
	$s =~ s/'/'\\''/g;
	return "'".$s."'";
}
sub read_expression {
	my ($profile,$log_file) = @_;
	my ($line,@ele);
	my $profile_tag = 0;
	my $total_rep = 0;
	my $total_point = 0;
	my %profile_expression;
	open PROFILE,"<",$profile or die "cannot open the profile ".$profile.": ".$!."\n";
	while($line=<PROFILE>){
		chomp $line;
		$line =~ s/\r$//;
		if($line =~ /repeat_number=/){
			@ele = split(/=/,$line);
			$total_rep = $ele[1];
			print STDERR "total replication:".$total_rep."\n";
		}elsif($line =~ /timepoint_number=/){
			@ele = split(/=/,$line);
			$total_point = $ele[1];
			print STDERR "total timepoint:".$total_point."\n";
		}
		elsif($line =~ /^rep/){
			@ele = split(/\t/,$line);
			my $tps = scalar @ele - 2;
			if($tps != $total_point){
				print STDERR $line."\n";die;
			}
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
	#sorted gene order, so the gene numbers (and the run.log mapping) are reproducible
	my @all_name = sort { natural_cmp($a,$b) } keys %{$profile_expression{1}};#get the gene name
	foreach my $in (keys %is_input){
		die "-input: ".$in." is not a row of the profile ".$profile."\n" if(!defined $profile_expression{1}{$in});
	}
	#the fitted genes first, then the inputs (the numbers of the genes do not change when an input is added)
	my @all_gene = ((grep { !$is_input{$_} } @all_name), (grep { $is_input{$_} } @all_name));
	$n_fit = scalar(@all_name) - scalar(keys %is_input);
	die "-input: no gene is left to fit\n" if($n_fit < 1);
	my $gene_no = 0;
	#the log is created here (mapping first); main() appends the run information
	open LOG,">",$log_file or die "cannot write ".$log_file.": ".$!."\n";
	print LOG "RUN_ID\tGENE_NAME\n";
	foreach my $i (@all_gene){
		#print STDERR $i."\n";
		$gene_index{$gene_no} = $i;
		print LOG $gene_no."\t".$i.($is_input{$i} ? "\tinput" : "")."\n";
		$gene_no++;
	}
	close LOG;
	if($cc == 1){
		#-c 1 uses the correlation coefficient, which is undefined (NaN) for a flat profile
		foreach my $rep (sort { $a <=> $b } keys %profile_expression){
			foreach my $g (@all_gene){
				next if(!defined $profile_expression{$rep}{$g});
				my %uniq = map { $_ => 1 } split(/\t/,$profile_expression{$rep}{$g});
				if(scalar(keys %uniq) <= 1){
					print STDERR "Warning: the profile of ".$g." in replication ".$rep." is flat, the correlation coefficient (-c 1) is undefined for it; use -c 0 if the runs fail\n";
				}
			}
		}
	}
	return ($gene_no,$total_rep,$total_point,%profile_expression);
}

sub read_knowledge {
	my ($total_gene_no,$know_init_file) = @_;
	my ($line,@ele);
	my $format_tag = 0;
	my %knowledge;
	my %fix_status_generation;
	my %known_gene = map { $gene_index{$_} => 1 } keys %gene_index;
	open KNOW,"<",$knowledge or die "cannot open the knowledge file ".$knowledge.": ".$!."\n";
	my $line_no = 0;
	while($line=<KNOW>){
		$line_no++;
		chomp $line;
		$line =~ s/\r$//;#CRLF files
		#print STDERR $line."\n";
		if($line =~ /^TF\tGENE\tREGULATORY\s*$/){
			print STDERR "format is correct\n";
			$format_tag = 1;
			next;
		}
		if($format_tag){
			next if($line =~ /^\s*$/);
			@ele = split(/\t/,$line);
			if(scalar @ele < 3){
				die "knowledge file ".$knowledge." line ".$line_no.": expected TF<tab>GENE<tab>REGULATORY\n";
			}
			my ($tf,$gene,$regulation) = ($ele[0],$ele[1],$ele[2]);
			$regulation =~ s/^\s+|\s+$//g;
			if($regulation !~ /^[?+\-0]$/){
				die "knowledge file ".$knowledge." line ".$line_no.": regulation must be one of ? + - 0, got '".$regulation."'\n";
			}
			if(!$known_gene{$tf} || !$known_gene{$gene}){
				print STDERR "Warning: knowledge file line ".$line_no.": ".$tf." or ".$gene." is not a gene of the profile, ignored\n";
				next;
			}
			if($is_input{$gene}){
				die "knowledge file ".$knowledge." line ".$line_no.": ".$gene." is an input (-input), it cannot be a target (GENE column)\n";
			}
			$knowledge{$tf}{$gene} = $regulation;
			#print STDERR "TF:".$tf."\tGene:".$gene."\treg:".$regulation."\n";
		}
	}
	close KNOW;
	if(!$format_tag){
		die "knowledge file ".$knowledge." has no header line 'TF<tab>GENE<tab>REGULATORY'\n";
	}
	open INIT,">",$know_init_file or die "cannot write ".$know_init_file.": ".$!."\n";
	for(my $i=0;$i<$n_fit;$i++){#each target gene (an input is not a target)
		my $gene_name = $gene_index{$i};
		my $gene_name_generation = $gene_name."_0";#0 is the initial generation
		my $regulation_knowledge;
		#print STDERR $gene_name.": ";
		for(my $j=0;$j<$total_gene_no;$j++){#each TF (the genes and the inputs)
			my $tf_name = $gene_index{$j};
			#print STDERR "TF:".$tf_name."\n";
			my $v = $knowledge{$tf_name}{$gene_name};
			$v = ($is_input{$tf_name} ? '0' : '?') if(!defined($v));#an input that is not in the knowledge file does not regulate
			print INIT $v." ";
			$regulation_knowledge .= $v." ";
		}
		$regulation_knowledge =~ s/ $//;
		$fix_status_generation{$gene_name_generation} = $regulation_knowledge;
		#print STDERR "fix_status_generation{".$gene_name_generation."}:".$fix_status_generation{$gene_name_generation}."\n";
		print INIT "\n";
	}
	print_input_rows(\*INIT,$total_gene_no);
	close INIT;
	return %fix_status_generation;
}
#Rows of the knowledge matrix of the inputs (they are not targets): no connection. The solver reads one row per variable.
sub print_input_rows {
	my ($fh,$total_gene_no) = @_;
	for(my $i=$n_fit;$i<$total_gene_no;$i++){
		print $fh join(" ",("0") x $total_gene_no)." \n";
	}
}
#S-system only: Hij are fixed to 0 (only Gij is used to determine the regulations)
sub generate_h_knowledge {
	my ($total_gene_no,$h_file) = @_;
	open HKNOW,">",$h_file or die "cannot write ".$h_file."\n";
	for(my $i=0;$i<$total_gene_no;$i++){
		print HKNOW join(" ",(0) x $total_gene_no)."\n";
	}
	close HKNOW;
}
sub generate_config {
	my ($configure,$total_gene_no,$total_rep,$total_point,%profile_expression) = @_;
	open CONF,">",$configure or die "cannot write the config file ".$configure.": ".$!."\n";
	print CONF "num_var ".$total_gene_no."\n";
	print CONF "num_tf ".$total_gene_no."\n";
	print CONF "num_run ".$total_rep."\n";
	print CONF "num_trial ".$total_point."\n";
	if($time_list ne ""){
		my @tp = split(/[,\s]+/,$time_list);
		die "-time has ".scalar(@tp)." values but the time series have ".$total_point." time points\n" if(scalar @tp != $total_point);
		for(my $q=0; $q<scalar @tp; $q++){
			die "-time: '".$tp[$q]."' is not a number\n" if(!looks_like_number($tp[$q]));
			die "-time must be strictly increasing\n" if($q > 0 && $tp[$q] <= $tp[$q-1]);
		}
		print CONF "time_points ".join(" ",map { sprintf("%.10g",$_/$time_unit) } @tp)."\n";
	}
	print CONF "add_perturbation yes ".$perturb_no." ".$noise_degree."\n" if($perturb_no > 0);
	print CONF "nochange 1000\n";
	if($model eq "HFODE"){
		print CONF "delta_t 1\n";
		print CONF "b_up 1\n";
		print CONF "b_down 0\n";
		print CONF "transMax_up 8\n";
		print CONF "transMax_down 0\n";
		print CONF "n_up ".$n_max."\n";
		print CONF "n_down -".$n_max."\n";
		print CONF "k_up 4\n";
		print CONF "k_down 0\n";
		print CONF "degrade_up 1\n";
		print CONF "degrade_down 0\n";
	}else{
		#S-system: dX/dt = alpha*prod(Xj^gij) - beta*prod(Xj^hij), time step is fixed to 1
		print CONF "alpha_up 4\n";
		print CONF "alpha_down 0\n";
		print CONF "beta_up 4\n";
		print CONF "beta_down 0\n";
		print CONF "g_up 3\n";
		print CONF "g_down -3\n";
		print CONF "h_up 3\n";
		print CONF "h_down -3\n";
		print CONF "deltat_up 1\n";
		print CONF "deltat_down 1\n";
	}
	print CONF "list_data\n";
	for(my $i=1;$i<=$total_rep;$i++){
		print CONF $i."\n";
		for(my $j=0;$j<$total_gene_no;$j++){
			my $gene_name = $gene_index{$j};
			if(!$profile_expression{$i}{$gene_name}){
			    print STDERR "i:".$i."\tgene_name:".$gene_name."\n";die;
			}
			print CONF $profile_expression{$i}{$gene_name}."\n";
		}
		#print CONF "\n";
	}

	close CONF or die "cannot write the config file ".$configure.": ".$!."\n";
}
#Read the result file of one gene (one line per run: parameters, fitness, one trailing field).
#Only valid runs are returned: enough columns, finite regulation parameters and a finite fitness
#(NaN/inf, e.g. from -c 1 with a flat profile, are skipped). Each run is [fitness, @columns].
sub read_valid_runs {
	my ($results) = @_;
	my ($n_start,$n_end) = regulation_columns(scalar keys %gene_index);
	my @runs;
	open my $fh,"<",$results or die "cannot open the result file ".$results.": ".$!."\n";
	while(my $line=<$fh>){
		chomp $line;
		$line =~ s/\r$//;
		my @ele = split(/ /,$line);
		next if(scalar @ele < ($n_end + 2));
		next if(!is_finite($ele[scalar(@ele) - 2]));
		my $ok = 1;
		for(my $c=$n_start; $c<$n_end; $c++){
			if(!is_finite($ele[$c])){ $ok = 0; last; }
		}
		next if(!$ok);
		push @runs, [$ele[scalar(@ele) - 2], @ele];
	}
	close $fh;
	return @runs;
}
#Fitness-weighted vote of all runs (not only the top ones) of the current iteration.
#weight of a run = best fitness / fitness of the run (0..1]; a regulation is used when |value| > min.
#VOTE_SCORE is sign-agnostic: it is the weighted fraction of runs in which the regulation is used,
#whatever its sign (+ or -), so a regulation that flips sign between runs still scores high.
#Returns a hash reference TF_GENE => score; the caller stores it in %vote_score when the regulation
#gets fixed (the score is kept because the result file <gene>.txt is overwritten by the quantification).
sub compute_vote_scores {
	my ($gene_name,$iga_results,$total_gene_no) = @_;
	my ($n_start,$n_end) = regulation_columns($total_gene_no);
	my @runs = read_valid_runs($iga_results);
	return {} if(!@runs);
	my $best;
	foreach my $r (@runs){ $best = $r->[0] if(!defined($best) || $r->[0] < $best); }
	my $eps = 1e-12;
	my ($wsum,@use) = (0);
	foreach my $r (@runs){
		my $w = ($best+$eps)/($r->[0]+$eps);
		$wsum += $w;
		for(my $c=$n_start; $c<$n_end; $c++){
			$use[$c-$n_start] += $w if(abs($r->[$c+1]) > $min);#$r->[0] is the fitness
		}
	}
	my %score;
	for(my $c=$n_start; $c<$n_end; $c++){
		my $tf_name = $gene_index{$c-$n_start};
		$score{$tf_name."_".$gene_name} = sprintf("%.3f",($use[$c-$n_start]||0)/$wsum);
	}
	return \%score;
}
#Run the EMA of the given genes in parallel; die (naming the gene) if any of the runs failed.
sub run_genes {
	my ($genes_ref,$gen,$use_know,$use_config) = @_;
	my @threads;
	foreach my $i (@{$genes_ref}){
		push @threads, [$i, threads->new(\&run_iga, $i,$gen,$use_know,$use_config)];
		sleep(1) while(scalar threads->list(threads::running) >= $threads);
	}
	my @errors;
	foreach my $t (@threads){
		my $err = $t->[1]->join();
		if(!defined($err)){
			$err = $gene_index{$t->[0]}.": the EMA thread terminated abnormally";
		}
		push @errors, $err if($err ne "");
	}
	die join("\n",@errors)."\n" if(@errors);
}
sub run_quantification {
	my ($gen,$use_know,$use_config,$total_gene_no,%confidence_level) = @_;
	my $n_start;
	my $n_end;
	my $gene_name;
	($n_start,$n_end) = regulation_columns($total_gene_no);
	my $final_results = $output_dir."/final_results.txt";
	#NOTE: the quantification rewrites <gene>.txt; the vote scores were therefore computed (and stored in
	#%vote_score) in the iteration where the regulation was fixed, they are not recomputed here.
	run_genes([0..($n_fit-1)],$gen,$use_know,$use_config);
	open FINAL,">",$final_results or die "cannot write ".$final_results.": ".$!."\n";
	print FINAL "TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\tVOTE_SCORE\n";
	for(my $i=0; $i<$n_fit; $i++){
		$gene_name = $gene_index{$i};
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results(1,$iga_results);#select top 1
		print STDERR "Gene:".$gene_name." ".$select_iga_results[0]."\n";
		my @iga = split(/\t/,$select_iga_results[0]);
		for(my $r=$n_start; $r<$n_end; $r++){
			my $tf_no = $r-$n_start;
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
			#regulations given by the knowledge file: confidence 1.000 (same as a user-given +/- in a partially unknown gene)
			my $conf = defined($confidence_level{$relationship}) ? $confidence_level{$relationship} : "1.000";
			my $vote = defined($vote_score{$relationship}) ? $vote_score{$relationship} : ($role eq "0" ? "0.000" : "1.000");#regulations given by the knowledge file have no vote
			print FINAL $tf_name."\t".$gene_name."\t".$role."\t".$conf."\t".$vote."\n";
		}
	}
	close FINAL or die "cannot write ".$final_results.": ".$!."\n";
}
sub run_EMA {
	my ($gen,$know_filename,$use_config,$total_gene_no,$fix_ref,$confidence_ref,$fix_status_ref,$fix_know_ref) = @_;
	my @fix_value = @{$fix_ref};
	my %confidence_level = %{$confidence_ref};
	my %fix_status_generation = %{$fix_status_ref};
	my %fix_know = %{$fix_know_ref};
	my $gene_name;
	my $cutoff = 0.8;
	my $use_know = $tmp_dir."/".$know_filename."_ForStep".$gen.".txt";
	my $next_generation = $gen + 1;
	my $new_know = $tmp_dir."/".$know_filename."_ForStep".$next_generation.".txt";

	print STDERR "Step2: GRN decomposition\n";
	print STDERR "Step3: Parallel solving\n";
	my @to_run;
	for(my $i=0; $i<$n_fit; $i++){
		if($fix_value[$i]){
			$gene_name = $gene_index{$i};
			print STDERR $gene_name ." is finished\n";
		}else{
			push @to_run, $i;
		}
	}
	run_genes(\@to_run,$gen,$use_know,$use_config);
	print STDERR "Step4: Regulation determination\n";
	open KNOW,">",$new_know or die "cannot write ".$new_know.": ".$!."\n";
	my ($n_start,$n_end) = regulation_columns($total_gene_no);
	for(my $i=0; $i<$n_fit; $i++){
		$gene_name = $gene_index{$i};
		if($fix_value[$i]){
			#fully fixed gene: it was not run (no <gene>.txt), its knowledge is carried forward
			#print STDERR "fix_knowledge{".$i."}:".$fix_know{$i}."\n";
			print KNOW $fix_know{$i}."\n";
			next;
		}
		my (@regulatory_p,@regulatory_n,@regulatory_z);
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results($select_top_no,$iga_results);#select top 5
		#vote scores of this iteration (kept for the regulations fixed in this iteration)
		my $vs = compute_vote_scores($gene_name,$iga_results,$total_gene_no);
		for(my $j=0; $j<$total_gene_no; $j++){
			$regulatory_p[$j] = 0;
			$regulatory_n[$j] = 0;
			$regulatory_z[$j] = 0;
		}
		#HFODE: n(i,j) ; S-system: g(i,j) (Hij are fixed to 0)
		for(my $j=0; $j<$select_top_no; $j++){
			my @iga = split(/\t/,$select_iga_results[$j]);
			for(my $r=$n_start; $r<$n_end; $r++){
				my $gene_no = $r-$n_start;
				if($iga[$r] > 0){
					if(abs($iga[$r]) > $min){
						$regulatory_p[$gene_no]++;
					}else{
						$regulatory_z[$gene_no]++;
					}
				}elsif($iga[$r] < 0){
					if(abs($iga[$r]) > $min){
						$regulatory_n[$gene_no]++;
					}else{
						$regulatory_z[$gene_no]++;
					}
				}else{
					$regulatory_z[$gene_no]++;
				}
			}
		}
		my $gene_generation_key = $gene_name."_".$gen;
		my $last_gene_generation_key = $gene_name."_".($gen - 1);
		#knowledge of the previous iteration: fixed regulations (+,-,0) are inherited, only '?' can change
		my @prev = split(/ /,$fix_status_generation{$last_gene_generation_key});
		my @row;
		for(my $j=0; $j<$total_gene_no; $j++){
			my $tf_name = $gene_index{$j};
			my $relationship = $tf_name."_".$gene_name;
			if(defined($prev[$j]) && $prev[$j] ne "?"){
				$row[$j] = $prev[$j];#Step7: inheritance, a fixed regulation cannot regress
				next;
			}
			my $each_gene_knowledge;
			if($regulatory_p[$j] && ($regulatory_p[$j]/$select_top_no) >= $cutoff){
				$each_gene_knowledge = "+";
			}
			if($regulatory_n[$j] && ($regulatory_n[$j]/$select_top_no) >= $cutoff){
				$each_gene_knowledge = "-";
			}
			if($regulatory_z[$j] && ($regulatory_z[$j]/$select_top_no) >= $cutoff){
				$each_gene_knowledge = "0";
			}
			if(defined($each_gene_knowledge)){
				if(!$confidence_level{$relationship}){
					$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
				}
			}else{
				$each_gene_knowledge = "?";
			}
			$row[$j] = $each_gene_knowledge;
			#print STDERR "TF-gene:".$tf_name."-".$gene_name." regulation:".$each_gene_knowledge."\n";
		}
		my $know_info = join(" ",@row);
		$fix_status_generation{$gene_generation_key} = $know_info;
		print STDERR "Check the status of the determined regulations\n";
		print STDERR "fix_status_generation{".$gene_generation_key."}:".$fix_status_generation{$gene_generation_key}."\n";
		print STDERR "fix_status_generation{".$last_gene_generation_key."}:".$fix_status_generation{$last_gene_generation_key}."\n";
		if($know_info eq $fix_status_generation{$last_gene_generation_key}){
			print STDERR "Get the same fix regulations, so need to fix more regulations\n";
			#fix the remaining '?' from the top-1 run (same abs(value) > min rule as the voting),
			#the regulations that are already fixed are kept
			my @select_top1 = select_results(1,$iga_results);#select top 1
			my @iga = split(/\t/,$select_top1[0]);
			for(my $j=0; $j<$total_gene_no; $j++){
				next if($row[$j] ne "?");
				my $v = $iga[$n_start+$j];
				my $role = "0";
				if(abs($v) > $min){
					$role = ($v > 0) ? "+" : "-";
				}
				$row[$j] = $role;
				my $relationship = $gene_index{$j}."_".$gene_name;
				if(!$confidence_level{$relationship}){
					$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
				}
			}
			$know_info = join(" ",@row);
			$fix_status_generation{$gene_generation_key} = $know_info;
			print STDERR "Update fix_status_generation{".$gene_generation_key."}:".$fix_status_generation{$gene_generation_key}."\n";
		}
		#vote score of every regulation that is fixed (first time seen as fixed) in this iteration
		for(my $j=0; $j<$total_gene_no; $j++){
			next if($row[$j] eq "?");
			my $relationship = $gene_index{$j}."_".$gene_name;
			if(!defined($vote_score{$relationship}) && defined($vs->{$relationship})){
				$vote_score{$relationship} = $vs->{$relationship};
			}
		}
		if(!$know_info){
			print STDERR "no value, error\n";die;
		}
		print KNOW $know_info."\n";
	}
	print_input_rows(\*KNOW,$total_gene_no);
	close KNOW or die "cannot write ".$new_know.": ".$!."\n";
	return (\%confidence_level, \%fix_status_generation);
}
#Returns "" on success, otherwise an error message naming the gene (die in a thread would be lost).
sub run_iga {
	my ($gene_no,$gen,$know,$conf) = @_;
	my $gene_name = $gene_index{$gene_no};
	my $out = $output_dir."/".$gene_name.".txt";
	my $ok = eval {
		my $src_dir = dirname($program);
		if($model eq "HFODE"){
			my $ema_HFODE = $src_dir."/EMA_HFODE/EMA_HFODE";
			die $ema_HFODE." does not exist\n" if(!-e $ema_HFODE);
			die $conf." does not exist\n" if(!-e $conf);
			my $command = shq($ema_HFODE)." -i ".$gene_no." -n 30 -m 0 -G ".$generation_no." -I ".$gen." -F ".$fitness_type." -c ".$cc." ".shq($conf)." ".shq($know)." > ".shq($out);
			#print STDERR $command."\n";
			my $status = system($command);
			die $gene_name.": ".$ema_HFODE." failed (exit status ".($status >> 8).", signal ".($status & 127).")\n" if($status != 0);
		}else{
			my $ema_ssystem = $src_dir."/EMA_SSYSTEM/EMA_ssystem";
			die $ema_ssystem." does not exist, please run make in EMA_SSYSTEM\n" if(!-e $ema_ssystem);
			die $conf." does not exist\n" if(!-e $conf);
			#EMA_ssystem takes the Gij knowledge file and the Hij knowledge file (all 0)
			#and writes the results to IGA_result<gene_no>_step<gen>.txt in the current directory (= tmp_dir, see main)
			my $result = "IGA_result".$gene_no."_step".$gen.".txt";
			unlink $result;
			my $command = shq($ema_ssystem)." -i ".$gene_no." -n 30 -m 0 -G ".$generation_no." -I ".$gen." -F ".$fitness_type." -c ".$cc." ".shq($conf)." ".shq($know)." ".shq($h_know_file)." > /dev/null";
			my $status = system($command);
			die $gene_name.": ".$ema_ssystem." failed (exit status ".($status >> 8).", signal ".($status & 127).")\n" if($status != 0);
			if(!-s $result){
				die $gene_name.": ".$ema_ssystem." produced no result (".$result.")\n";
			}
			#rename fails across filesystems, fall back to copy+unlink
			move($result,$out) or die $gene_name.": cannot move ".$result." to ".$out.": ".$!."\n";
		}
		die $gene_name.": result file ".$out." is missing or empty\n" if(!-s $out);
		my $valid = scalar(read_valid_runs($out));
		if($valid < $select_top_no){
			die $gene_name.": only ".$valid." valid result lines (finite fitness, all parameters) in ".$out.", at least ".$select_top_no." are needed".(($cc == 1) ? " (-c 1 gives NaN for a flat profile, try -c 0)" : "")."\n";
		}
		#keep the GA results of every iteration (<output>/iga_results/<gene>_step<k>.txt) because <gene>.txt is overwritten by the next iteration;
		#GREMA_KEEP_ITER=0 switches this off
		if(!defined($ENV{GREMA_KEEP_ITER}) || $ENV{GREMA_KEEP_ITER} ne "0"){
			my $keep_dir = $output_dir."/iga_results";
			mkdir $keep_dir if(!-d $keep_dir);
			copy($out,$keep_dir."/".$gene_name."_step".$gen.".txt") or die $gene_name.": cannot keep the result of iteration ".$gen.": ".$!."\n";
		}
		1;
	};
	return "" if($ok);
	my $err = $@;
	$err = "unknown error" if(!defined($err) || $err eq "");
	chomp $err;
	return $err;
}
#Top runs of a result file: the lines (tab separated parameters, without fitness) sorted by fitness.
#Runs with equal fitness are all kept (ties keep the run order); invalid (NaN/inf/short) lines are skipped.
sub select_results {
	my ($top_no,$results) = @_;
	my @select_results;
	my @runs = read_valid_runs($results);
	if(scalar @runs < $top_no){
		die "result file ".$results." has only ".scalar(@runs)." valid runs, ".$top_no." needed\n";
	}
	my @pairs;#[fitness, line, run order]
	my $no = 0;
	foreach my $r (@runs){
		my @ele = @{$r};
		shift @ele;#fitness (already known)
		my $all_parameter_no = scalar @ele;
		push @pairs, [$r->[0], join("\t",@ele[0..($all_parameter_no - 3)]), $no++];
	}
	my @sorted = sort { $a->[0] <=> $b->[0] || $a->[2] <=> $b->[2] } @pairs;
	for(my $i=0; $i<$top_no; $i++){
		$select_results[$i] = $sorted[$i]->[1];
	}
	return @select_results;
}
sub check_knowledge {
	my ($know,$gene_no, $fix_know_ref) = @_;
	my %fix_know = %{$fix_know_ref};
	my ($line,@ele);
	my $no = 0;
	open (KNOW,"<",$know) or die "open file ".$know." error\n";
	while($line=<KNOW>){
		chomp $line;
		if($no == $gene_no){
			#print STDERR $line."\n";
			my $unfix_no = 0;
			@ele = split(/ /,$line);
			foreach my $i(@ele){
				if($i eq "?"){
					$unfix_no++;
				}
			}
			close KNOW;
			if($unfix_no == 0){
				$fix_know{$no} = $line;
				return (1, \%fix_know);
			}else{
				$fix_know{$no} = "unfix";
				return (0, \%fix_know);
			}
		}
		$no++;
	}
	close KNOW;
	die "knowledge file ".$know." has no row for gene ".$gene_no."\n";
}
sub main {
	my $command;
	my $generation = 1;#first generation;
	my $log = $output_dir."/run.log";
	my $knowledge_filename_tmp = basename($knowledge);
	my $knowledge_filename;
	my %fix_status;
	my $fix_ref;
	my %fix_knowledge;
	my $fix_knowledge_ref;
	#generated side files go to the tmp folder of the output (no race between runs, input folder may be read-only)
	my $config = $tmp_dir."/".basename($expression)."_config";
	my @fix;
	my $all_fix = 0;
	my $total_fix_no = 0;
	my %confidence;
	my $confidence_ref;
	if($knowledge_filename_tmp =~ /\./){
		my @filename = split(/\./,$knowledge_filename_tmp);
		$knowledge_filename = $filename[0];
	}else{
		$knowledge_filename = $knowledge_filename_tmp;
	}
	my $use_knowledge = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
	print STDERR "know_init:".$use_knowledge."\n";
	my ($total_gene,$total_repeat,$total_data_points,%profile) = read_expression($expression,$log);
	print STDERR "Number of gene:".$n_fit."\n";
	print STDERR "Number of inputs:".($total_gene-$n_fit)."\n" if($total_gene > $n_fit);
	%fix_status = read_knowledge($total_gene,$use_knowledge);#generate the initial knowledge
	if(!-e $use_knowledge){
		print STDERR "generate initial knowledge failed <".$use_knowledge.">\n";die;
	}
	generate_config($config,$total_gene,$total_repeat,$total_data_points,%profile);
	if($model eq "s-system"){
		$h_know_file = $tmp_dir."/h_knowledge_zero.txt";
		generate_h_knowledge($total_gene,$h_know_file);
	}

	#the mapping RUN_ID->GENE_NAME was written by read_expression, append the run information
	open LOG,">>",$log or die "cannot write ".$log.": ".$!."\n";
	print LOG "Total gene:".$n_fit."\n";
	print LOG "Inputs:".join(",",map { $gene_index{$_} } ($n_fit..($total_gene-1)))."\n" if($total_gene > $n_fit);
	print LOG "Total repeat:".$total_repeat."\nTotal data points:".$total_data_points."\n";
	print LOG "Use know init:".$use_knowledge."\nUse config:".$config."\nUse input data:".$expression."\n";
	print LOG "Use Fitness".$fitness_type."\nUse threads:".$threads."\nUse generation:".$generation_no."\nUse CC:".$cc."\nUse min:".$min."\n";
	close LOG;

	print STDERR "Step1: Initialisation\n";
	for(my $i=0;$i<$n_fit;$i++){
		($fix[$i], $fix_knowledge_ref) = check_knowledge($use_knowledge,$i,\%fix_knowledge);
		%fix_knowledge = %{$fix_knowledge_ref};
		$total_fix_no += $fix[$i];
	}
	print STDERR "Initial Generation ".$generation." fix: ".$total_fix_no."\n";
	if($total_fix_no == $n_fit){
		$all_fix = 1;
	}

	#=======finished initial=================================================
	if(!$all_fix){
		do{
			print STDERR "Generation ".$generation." start running\n";
			if(!-e $use_knowledge){
				print STDERR "use knowledge does not exist <".$use_knowledge.">\n";die;
			}
			($confidence_ref,$fix_ref) = run_EMA($generation,$knowledge_filename,$config,$total_gene,\@fix,\%confidence,\%fix_status,\%fix_knowledge);
			%confidence = %{$confidence_ref};
			%fix_status = %{$fix_ref};
			$generation++;
			my $new_knowledge_file = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
			$total_fix_no = 0;
			print STDERR "Step5: GRN combination\n";
			for(my $i=0;$i<$n_fit;$i++){
				($fix[$i],$fix_knowledge_ref) = check_knowledge($new_knowledge_file,$i,\%fix_knowledge);
				%fix_knowledge = %{$fix_knowledge_ref};
				#print STDERR "fix[".$i."]:".$fix[$i]."\n";
				#print STDERR "fix_knowledge{".$i."}:".$fix_knowledge{$i}."\n";
				$total_fix_no += $fix[$i];
			}
			my $last_gen = $generation - 1;
			print STDERR "Generation ".$last_gen." finish, fix: ".$total_fix_no."\n";
			print STDERR "Step6: Qualitative termination test\n";
			if($total_fix_no == $n_fit){
				$all_fix = 1;
				print STDERR "All regulations are determined, go to Step8\n";
			}else{
				print STDERR "Step7: Inheritance\n";
			}
		}while(!$all_fix)
	}
	my $final_knowledge = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
	if(!-e $final_knowledge){
		print STDERR "use final knowledge does not exist <".$final_knowledge.">\n";die;
	}
	print STDERR "Step8:GRN quantification\n";
	run_quantification($generation,$final_knowledge,$config,$total_gene,%confidence);
	#calprofile*, progress* and printdata.txt written by the EMA are already in the tmp folder (main runs inside it)
}
#=========subroutine end================================
GetOptions(
	'i=s'	=>\$expression,
	'o=s'	=>\$output_dir,
	'k=s'	=>\$knowledge,
	't=i'	=>\$threads,
	'm=s'	=>\$model,
	'g=i'	=>\$generation_no,
	'f=i'	=>\$fitness_type,
	'c=i'	=>\$cc,
	'min=s'	=>\$min,
	'perturb=i'	=>\$perturb_no,
	'noise=f'	=>\$noise_degree,
	'nmax=f'	=>\$n_max,
	'time=s'	=>\$time_list,
	'input=s'	=>\$input_list,
	'tunit=f'	=>\$time_unit,
	'h'	=>\$help,
); 
if(!$expression ||!$knowledge || !$threads || $help ||$threads < 1 || !$output_dir){
	Usage();die;
}
if($perturb_no < 0 || $noise_degree < 0 || ($perturb_no > 0 && $noise_degree <= 0)){
	print STDERR "-perturb must be >= 0 and needs -noise > 0\n";
	Usage();die;
}
if($perturb_no > 0 && $model ne "HFODE"){
	print STDERR "-perturb is only supported with the HFODE model\n";
	Usage();die;
}
if($input_list ne ""){
	if($model ne "HFODE"){
		print STDERR "-input is only supported with the HFODE model\n";
		Usage();die;
	}
	if($perturb_no > 0){
		print STDERR "-input cannot be used with -perturb\n";
		Usage();die;
	}
	foreach my $in (split(/[,\s]+/,$input_list)){
		next if($in eq "");
		die "-input: ".$in." is given twice\n" if($is_input{$in});
		$is_input{$in} = 1;
	}
	if(!%is_input){
		print STDERR "-input needs at least one name\n";
		Usage();die;
	}
}
if($time_unit <= 0){
	print STDERR "-tunit must be > 0\n";
	Usage();die;
}
if($time_list ne "" && $model ne "HFODE"){
	print STDERR "-time is only supported with the HFODE model\n";
	Usage();die;
}
if($n_max <= 0){
	print STDERR "-nmax must be > 0\n";
	Usage();die;
}
if($min < 0){
	print STDERR "min value must be >= 0\n";
	Usage();die;
}
if(($cc != 0)and($cc != 1)){
	print STDERR "type of cc is error ".$cc."\n";
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
	print STDERR "type of model is error\n";
	Usage();die;
}
if($generation_no < 100){
	print STDERR "generation must be > 100\n";
	Usage();die;
}
if(($fitness_type != 0)and($fitness_type != 1)and($fitness_type != 2)and($fitness_type != 3)and($fitness_type != 4)){
	print STDERR "type of fitness function is error\n";
	Usage();die;
}
if(-d $output_dir){
	$output_dir =~ s/\/$//;
	my $final = $output_dir."/final_results.txt";
	if(-e $final){
		print STDERR "This output folder already has the final_results.txt\n";die;
	}
	$tmp_dir = $output_dir."/tmp";
	if(!-d $tmp_dir){
		print STDERR "create tmp folder:".$tmp_dir."\n";
		mkdir $tmp_dir or die "cannot create ".$tmp_dir.": ".$!."\n";
	}
	#absolute paths, because main() runs inside the tmp folder (the EMA programs write their side files to the cwd)
	$expression = abs_path($expression);
	$knowledge = abs_path($knowledge);
	$output_dir = abs_path($output_dir);
	$tmp_dir = abs_path($tmp_dir);
	chdir($tmp_dir) or die "cannot enter ".$tmp_dir.": ".$!."\n";
}else{
	print STDERR "Directory of output does not exist\n";
	die;
}
main();

