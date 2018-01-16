#!/soft/bin/perl -w                                                                                                                                                    
use strict;
use warnings;

my ($file) = @ARGV;

our $verbose   = 0 unless $verbose;

unless($file){
    print STDERR "Usage: $0 file\n";
    exit(0);
}

open(IN,"<$file") or die("cannot open $file");
my %list;
while(<IN>){
    chomp;
    #chr17   59763856        59763863        BRIP1   AGGCAATT        -       chr17   59758627        59940882        7       675     7.92288573222   1.48020882342e-15	0.5
    #chr17   59763857        59763864        BRIP1   AAGGCAAT        -       chr17   59758627        59940882        7       675     7.9453201786     1.48020882342e-15  0.5
    my ($chr, $start, $end, $gene, $kmer, $strand, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut, $score, $pval, $pval_adj)  = split /\t/, $_;
    push( @{$list{$chr}{$strand}}, 
	 [$chr, $start, $end, $gene, $kmer, $strand, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut, $score,$pval, $pval_adj] );
}	

foreach my $chr (keys %list ){
    foreach my $strand (keys %{$list{$chr}}){
	my $ranges = $list{$chr}{$strand};
	my ($clusters, $cluster_starts, $cluster_ends) = cluster($ranges,1,2);
		
	for(my $i=0; $i<scalar(@$clusters); $i++){
	    my $cluster       = $clusters->[$i];
	    my $cluster_start = $cluster_starts->[$i];
	    my $cluster_end   = $cluster_ends->[$i];
	    
	    if ($verbose){
		print "Cluster $cluster_start $cluster_end\n";
	    }

	    my @pvals;
	    my @scores;
	    my %labels;
	    my %gene;
	    # $chr, $start, $end, $gene, $kmer, $strand, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut, $score,$pval, $pval_adj
	    foreach my $range (@$cluster){
		if (my $verbose2){
		    my $s = join "\t", ("range", @$range);
		    print $s."\n";
		}
		my ($chr, $start, $end, $gene, $kmer, $strand, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut, $score, $pval, $pval_adj) = @$range;
		
		#my $gene_id = $gene.",".$gene_chr.":".$gene_start."-".$gene_end.":".$gene_mut;
		my $gene_id = $gene.",".$gene_chr.":".$gene_start."-".$gene_end.":".$mut.":".$gene_mut;
		$gene{$gene_id}++;
		push(@pvals, $pval_adj);
		push(@scores, $score);
	    }
	    my $mean_score = mean(@scores);
	    my $simes_p    = simes_p(@pvals);
	    my @labels     = sort {$a cmp $b} keys %labels;
	    
	    #my ($gene, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut);
	    my ($gene, $gene_chr, $gene_start, $gene_end, $mut, $gene_mut);
	    my @genes;
	    foreach my $id (keys %gene){
		$id =~/(\S+),(chr\S+):(\d+)-(\d+):(\d+):(\d+)/;
		#$id =~/(\S+),(chr\S+):(\d+)-(\d+):(\d+)/;
		push( @genes, $1);
		$gene_chr   = $2;
		$gene_start = $3 unless $gene_start;
		$gene_start = $3 if $3 < $gene_start;
		$gene_end   = $4 unless $gene_end;
		$gene_end   = $4 if $4 > $gene_end;
		$mut = $5 unless $mut ;  #added 
		$mut   = $5 if $5 > $mut; #added
		$gene_mut   = $6 unless $gene_mut; #added
		$gene_mut   = $6 if $6 > $gene_mut; #added
		#$gene_mut   = $5 unless $gene_mut;
		#$gene_mut   = $5 if $5 > $gene_mut;
	    }
	    if (scalar(@genes) > 1){
		$gene = join ":",@genes;
	    }
	    else{
		$gene = $genes[0];
	    }
	    
	    #my $s = join "\t", ($chr, $cluster_start, $cluster_end, $gene, ".", $strand, $gene_chr, $gene_start, $gene_end, $gene_mut ,$mean_score, $simes_p);
	    my $s = join "\t", ($chr, $cluster_start, $cluster_end, $gene, ".", $strand, $gene_chr, $gene_start, $gene_end, $gene_mut, $mut ,$mean_score, $simes_p);
	    print $s."\n";
	}
    }
}

sub simes_p{
    my @pvals = @_;
    #Rank the p-values of the overlapping 7mer windows in increasing order
    #	Pr, r=1,2,…n. Calculate Ps = min{ n P1/ 1,  n P2/2,  n P3/3, ... ,   nPn/n}
    #(where P1 is the lowest and Pn is the highest p-value in the cluster)
    my @sorted_pvals = sort {$a <=> $b} @pvals;
    my @rescaled_pvals;
    for(my $i=0; $i<scalar(@sorted_pvals); $i++){
	$sorted_pvals[$i] = scalar(@sorted_pvals) * $sorted_pvals[$i] / ($i+1);
    }
    return min(@sorted_pvals);
}



sub min{
    my @v = @_;
    my $min = $v[0];
    foreach my $v (@v){
        $min = $v if $v < $min;
    }
    return $min;
}

sub mean{
    my (@v) = @_;
    my $sum = 0;
    foreach my $v (@v){
        $sum += $v;
    }
    return $sum/scalar(@v);
}

sub cluster{
    my ($loci, $s, $e) = @_;

    # we generate cluster of loci                                                                                                                            
    my $locus_clusters;

    # sort loci by start coordinate in ascending order                                                                                                     
    # and when equal, in descending end coordinate
    my @sorted_loci = sort { $a->[$s] <=> $b->[$s] } @$loci;
    
    # start and end coordinates of the clusters                                                                                                            
    my @cluster_start;
    my @cluster_end;

    # Create the first locus_cluster with the first locus                                                                                                    
    my $locus_cluster = [];
    push( @$locus_clusters, $locus_cluster);

    # we go over all of them in sorted order (from left to right)                                                                                            
    my $cluster_count = 0;
    my $count = 0;

  LOCUS:
    foreach my $locus ( @sorted_loci ){
        if ($count == 0){
            # Create the first locus_cluster with the first locus                                                                                            
            push( @$locus_cluster, $locus);
            $cluster_start[0]  =  $locus->[$s];
            $cluster_end[0]    =  $locus->[$e];
            $count++;
            if ($verbose){
                print  "cluster_count: $cluster_count\n";
                print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
            }
            next LOCUS;
        }
        if ($verbose){
            print  "cluster_count: $cluster_count\n";
            print  "cluster_start: $cluster_start[$cluster_count] cluster_end: $cluster_end[$cluster_count]\n";
        }
        
        # test overlap                                                                                                                                       
        if ( !( $locus->[$e] < $cluster_start[$cluster_count] || $locus->[$s] > $cluster_end[$cluster_count]) ){
            # add locus to cluster                                                                                                                           
            print "added\n" if $verbose;
            push( @$locus_cluster, $locus );
            
            # update start and end of cluster if necessary                                                                                                   
            if ($locus->[$s] < $cluster_start[$cluster_count]) {
                $cluster_start[$cluster_count] = $locus->[$s];
            }
            if ($locus->[$e] > $cluster_end[$cluster_count]) {
                $cluster_end[$cluster_count]   = $locus->[$e];
            }
        }
        else{
            # then we proceed with the next one:                                                                                                             
            # create new cluster (same variable name, new memory address!!!!)                                                                                
            $locus_cluster = [];
            push( @$locus_clusters, $locus_cluster);
            
            # add locus in new cluster                                                                                                                       
            push( @$locus_cluster, $locus );
            $cluster_count++;
            $cluster_start[$cluster_count] = $locus->[$s];
            $cluster_end[$cluster_count]   = $locus->[$e];
        }
    }
    return( $locus_clusters, \@cluster_start, \@cluster_end);
}
