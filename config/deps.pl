#! /usr/bin/perl

# Syntax: 
#       cd ppcode/bin
#       ./deps.pl
# or (better)
#       make deps

$verbose=0;
if ($ARGV[0] eq "-v") {
        shift;
        $verbose=1;
}

# Read and process a list of files names from stdin
$n = 0;
while (<>) {
        chomp;
        $f = $_;
        open FILE, "<", $f;
        s:.*/::;
        s/f90$/o/;
        $o=$_;
        while (<FILE>) {
                chomp;
                s/,.*//;
                @w=split;
                $h = lc($w[0]);
                if ($h eq "module") {
                        $m = lc($w[1]);
                        $defd{$o} = $m;
                        $defs{$m} = $o;
                        print "  $o defines module $m\n" if $verbose > 0;
                }
                if ($h eq "use") {
                        $m = lc($w[1]);
                        $uses{$o} = $m;
                        $obj[$n] = $o;
                        $mod[$n] = $m;
                        $n++;
                        print "  $o uses $m\n" if $verbose > 0;
                }
        }
}

print "uses:\n" if $verbose > 0;
print "# This dependency file was created automatically by \$(TOP)/config/deps.pl.\n";
print "# Do not edit manually -- instead run 'make deps > Makefile.dep' !\n\n";

foreach $f (sort(keys %uses)) {                       # $f is the object file we will build a list for
        print "building list for $f\n" if $verbose > 0;
        %list = ();
        for ($i=0; $i<$n; $i++) {
                $o = $obj[$i];
                $m = $mod[$i];
                if ($o eq $f) {
                        $d = $defs{$m};
                        if ($d eq "") {
                                print "  $f uses $m, which is define by an unknown\n" if $verbose > 0;
                        } elsif ($m eq "omp_lib") {
                                print "  $f uses $m, which we ignore\n" if $verbose > 0;
                        } else {
                                print "  $f uses $m, which is defined by $d, so we " if $verbose > 0;
                                if ($f eq $d) {
                                        print "ignore that\n" if $verbose > 0;
                                } else {
                                        print "add that\n" if $verbose > 0;
                                        $list{$d} = $f;
                                }
                        }
                }
        }
        $string = "";
        foreach $d (keys %list) {
                $string = "$string $d";
        }
        write if $string ne "";
}

format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<:@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$f                , $string
.
