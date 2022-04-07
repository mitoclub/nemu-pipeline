#!/usr/bin/perl -w

use LWP::Simple;
use LWP::UserAgent;
use File::Temp qw/ tempdir /;
use Cwd qw(getcwd);

#$pwd = getcwd();
#print "$pwd\n";

$query = "$ARGV[0]\[orgn\]+AND+$ARGV[1]";
$db1 = 'nucleotide';
$dir = "$ARGV[2]";

#assemble the esearch URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=$db1&term=$query&usehistory=y&retmax=100"; #5000 max

#post the esearch URL
$output = get($url);

#parse IDs retrieved
while ($output =~ /<Id>(\d+?)<\/Id>/sg) {
   push(@ids, $1);
   }

#assemble  the elink URL as an HTTP POST call
$url = $base . "elink.fcgi";

$url_params = "efetch.fcgi?db=$db1&rettype=gb";
foreach $id (@ids) {      
  $url_params .= "&id=$id";
  }
#print "$url_params\n";

$url = $base . $url_params;
$seqs = get($url);

$tmdir = tempdir( CLEANUP => 0, DIR => $dir );
chdir "$tmdir";
$pwd = getcwd();
print "temporary dir:: $pwd ::\n";

open (OUT, ">gbk") || die "Can't open file!\n";
print OUT "$seqs\n";
close OUT;
