__all__ = ['Params']

from bunch import Bunch

class Params (Bunch):
  """
  This class parses a Fortran namelist and converts it into a Bunch. It uses
  Wolfgang Dobler's Perl module Fortran::F90Namelist::Group, which is available
  on CPAN.
  """
  def __init__ (self, nml = "Data/params.nml"):

    from subprocess import check_output
    from sys import argv
    from re import findall
    from os.path import realpath, dirname
    from os import environ

    # Perl script that namelist using Wolfgang's library
    script = """
    use Fortran::F90Namelist::Group;
    my $nlgrp = Fortran::F90Namelist::Group->new();
    $nlgrp->parse (file => '%s');
    while (my $nl = $nlgrp->pop()) {print $nl->output (format=>'idl');}
    """ % nml

    # Make sure the Perl interpreter uses the right module
    dirname = dirname (realpath (__file__))
    env = environ.copy ()
    env['PERL5LIB'] = dirname + '/../../Perl:' + env.get ('PERL5LIB', '')

    # Execute script and get output
    output = check_output (['perl','-e',script], env = env)

    # Put quotation marks around dictionary keys
    keys = findall ('\s\S*:', output)
    for key in keys:
      output = output.replace (key, '"' + key[1:-1] + '":')

    # Replace logicals
    output = output.replace ('0L', 'False')
    output = output.replace ('-1L', 'True')

    # Join every line that ends with a continuation character
    # with the one that follows
    lines = ' '.join ([line.strip () for line in output.split ('$\n')])

    # Create a bunch of bunches of parameters
    dictionary = {}
    for line in lines.split ('\n'):
      (section, bunch) = line.split ('=')
      cmd = 'dictionary["' + section.strip () + '"] = Bunch (**' + bunch + ')'
      exec (cmd)
    Bunch.__init__ (self, **dictionary)

    # Strip string valued parameters of white spaces
    for k1, v1 in self.iteritems ():
      for k2, v2 in v1.iteritems ():
	if type (v2) == str: self[k1][k2] = v2.strip ()
