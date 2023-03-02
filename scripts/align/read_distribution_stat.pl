#/usr/bin/perl
# @sample=(S5_5);
# foreach $cell (@sample){
# print("in");
open(IN, "$ARGV[0]");
open(OUT, ">$ARGV[1]");
<IN>;
<IN>;
$a = <IN>;
@a = split(/\s+/, $a);
$total = $a[3];
print($a[3]);
<IN>;
<IN>;

chomp($a = <IN>);
@a = split(/\s+/, $a);
$exon_cds = $a[2];

chomp($a = <IN>);
@a = split(/\s+/, $a);
$exon5_count = $a[2];

chomp($a = <IN>);
@a = split(/\s+/, $a);
$exon3_count = $a[2];

chomp($a = <IN>);
@a = split(/\s+/, $a);
$introns_count = $a[2];
<IN>;
<IN>;

chomp($a = <IN>);
@a = split(/\s+/, $a);
$tss_up = $a[2];
<IN>;
<IN>;

chomp($a = <IN>);
@a = split(/\s+/, $a);
$tss_down = $a[2];
print OUT "total_count\t".$total."\n";
print OUT "exon_cds\t".$exon_cds."\n";
print OUT "exon5_count\t".$exon5_count."\n";
print OUT "exon3_count\t".$exon3_count."\n";
print OUT "introns_count\t".$introns_count."\n";
print OUT "tss_up_up\t".$tss_up."\n";
print OUT "tss_down\t".$tss_down."\n";
close(IN);
close(OUT);
