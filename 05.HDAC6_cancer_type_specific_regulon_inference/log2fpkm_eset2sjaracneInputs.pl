#!/usr/bin/perl

use strict;
use warnings;

my $dir = "/research/rgs01/project_space/yu3grp/Network_JY/yu3grp/generatedNetworks/TCGA/sjaracne.OL201905";
open (OUT1, "> $dir/log2fpkm_eset2sjaracneInputs.stat") or die;
print OUT1 "DataSet\tnSamples\texp_nIsoform\texp_nGene\ttf_nIsoform\ttf_nGene\tsig_nIsoform\tsig_nGene\n";
open (OUT2, "> $dir/log2fpkm_runSJARACNe.sh") or die;

opendir (DIR, "$dir/log2fpkm") or die;
while (my $file = readdir(DIR)) {
    next unless ($file =~ /\w/);
    print "$file is in progress...\n";

    my @F = split(/\./, $file); my $line1 = $F[0];
    my @G = split(/_/, $F[1]); $line1 .= "\t$G[2]\t$G[0]\t$G[1]";

    opendir (DIR1, "$dir/log2fpkm/$file/tf") or die;
    my $tf = "";
    while (my $dir1 = readdir(DIR1)) {
        next unless $dir1 =~ /\.tf\.txt$/;
        $tf = $dir1;
        my @F1 = split(/\./, $dir1);
        my @G1 = split(/\_/, $F1[1]); $line1 .= "\t$G1[0]\t$G1[1]";
    }
    closedir DIR1;

    opendir (DIR2, "$dir/log2fpkm/$file/sig") or die;
    my $sig = "";
    while (my $dir2 = readdir(DIR2)) {
        next unless $dir2 =~ /\.sig\.txt$/;
        $sig = $dir2;
        my @F2 = split(/\./, $dir2);
        my @G2 = split(/\_/, $F2[1]); $line1 .= "\t$G2[0]\t$G2[1]";
    }
    closedir DIR2;

    print OUT1 "$line1\n";
    print OUT2 "sjaracne $F[0] $dir/log2fpkm/$file/$file\.exp $dir/log2fpkm/$file/tf/$tf $dir/log2fpkm/$file/tf/ --host LSF --bootstrap 300 --c_threshold 1e-7 --queue standard --resource 2000 16000 100000 8000\n";
    print OUT2 "sjaracne $F[0] $dir/log2fpkm/$file/$file\.exp $dir/log2fpkm/$file/sig/$sig $dir/log2fpkm/$file/sig/ --host LSF --bootstrap 300 --c_threshold 1e-7 --queue standard --resource 2000 16000 100000 8000\n";
}
closedir DIR;

