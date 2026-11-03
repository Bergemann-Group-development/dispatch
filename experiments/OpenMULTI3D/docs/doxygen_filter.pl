# Doxygen input filter for the Fortran sources. Doxygen reads '%' before a word as
# its no-link marker and drops it, so self%x would render as selfx; escape every
# '%' in the !> and !< comments. Code and plain comments pass unchanged.
while (<>) {
  s{(!(?:>|<).*)}{(my $c = $1) =~ s/%/\\%/g; $c}e;
  print;
}
