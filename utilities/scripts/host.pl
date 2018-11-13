#!/usr/bin/perl

# debug flag (ignore all other command line arguments, except the last one)
while ($#ARGV > 0) {
	if ($ARGV[0] eq "-d") {$DEBUG="yes"; shift};
	shift;
}
if ($ARGV[0] eq "-d") {$DEBUG="yes"; shift};

# Get hostname (or a fake one as the last command line argument)
$L = `uname -n`;
if ($#ARGV == 0) {$H=$ARGV[0]} else {$H=$L};

# Remove newline
chomp($H);

# Try to translate to generic name
#######################################
$HOST=$H;

# NBI
$HOST="comp"   if $H =~ m/^comp/;
$HOST="comp"   if $H =~ m/astro.ku.dk/;
$HOST="nbi"    if $H =~ m/nbi.ku.dk/;

# Steno
$HOST="steno"  if $H =~ m/dcsc.ku.dk/;
$HOST="steno"  if $H =~ m/hpc.ku.dk/;
$HOST="steno"  if $H =~ m/^fend0[0-9]/;
$HOST="steno"  if $H =~ m/^astro0[4-9]/;
$HOST="steno"  if $H =~ m/^node[0-9][0-9][0-9].cluster/;
$HOST="astro01" if $L =~ m/^astro01/;
$HOST="astro01" if $L =~ m/^fend01/;
$HOST="fend02"  if $H =~ m/^fend02/;

# Pleiades
$HOST="pleiades" if $H =~ m/^pfe[0-9]/;
$HOST="pleiades" if $H =~ m/^mfe[0-9]/;
$HOST="pleiades" if $H =~ m/^bridge[0-9]/;
$HOST="pleiades" if $H =~ m/^r[0-9][0-9][0-9]i[0-9]n[0-9]/;

# SuperMUC
$HOST="supermuc" if $H =~ m/^login0[0-9]/;

# SuperMIG
$HOST="supermuc" if $H =~ m/^i0[0-9]r[0-9][0-9]s[0-9][0-9]/;

# Juqueen
$HOST="juqueen"  if $H =~ m/juqueen.*.zam.kfa-juelich.de/;

# Curie
$HOST="curie"    if $H =~ m/curie.tgcc.ccc.cea.fr/;

# Marconi
$HOST="marconi"  if $H =~ m/^r[0-9][0-9][0-9]u[0-9][0-9]l[0-9][0-9]/;
$HOST="marconi"  if $H =~ m/^r[0-9][0-9][0-9]c[0-9][0-9]s[0-9][0-9]/;

# Cygwin
$ARCH=`uname -s`;
chomp($ARCH);
if ($ARCH ne "Darwin") {
  $S=`uname -o`;
  $HOST="cygwin"   if $S =~ m/Cygwin/;
}

# Output result
#######################################
if ($DEBUG eq 'yes') {
  print "Fully qualified hostname is ";
  print $H; print "\n";
  print "Generic hostname is ";
}
print $HOST;
if ($DEBUG eq 'yes') {
  print "\n"
}
