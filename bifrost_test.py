#!/usr/bin/env python
"""Summary
Attributes:
    screen_width (int): Width of the printouts
"""
# encoding: utf-8

"""
bifrost_test.py

Script to auto-test Bifrost experiments.

Created by Mikolaj Szydlarski on 2021-02-13.
Copyright (c) 2013 ITA - University of Oslo
"""

import click
import os
import re
import subprocess
import sys
import yaml

from box import Box
from enum import Enum

screen_width = 77

class ExecuteState(Enum):
  not_tested      = 0
  no_errors       = 1
  with_warnings   = 2
  with_errors     = 3

def print_state(state):
  if state == ExecuteState.not_tested:
    return "[ SKIPPED ]"
  elif state == ExecuteState.no_errors:
    return "[ ALL  OK ]"
  elif state == ExecuteState.with_warnings:
    return "[ WARNING ]"
  elif state == ExecuteState.with_errors:
    return "[ ERRORS! ]"
  else:
    return "[ NOT SET ]"

def test_passed(state):
  if state == ExecuteState.not_tested:
    return 1
  elif state == ExecuteState.with_warnings:
    return 1
  elif state == ExecuteState.no_errors:
    return 1
  else:
    return 0

def generate_TR(cfg):
  # This will generate TR dictionary from cfg
  TR = {} # dictionary to store tests results
  for eid, exp in cfg.items():
    TR[eid] = {}
    TR[eid]["tests"] = {}
    TR[eid]["comp_res"] = ExecuteState.not_tested
    for no, unit in  exp.tests.items():
      TR[eid]["tests"][no] = {}
      TR[eid]["tests"][no]["descption"] = unit.descption
      TR[eid]["tests"][no]["result"] = ExecuteState.not_tested
  return TR

def do_we_have_tests_at_our_level(cfg,lvl) -> bool:
  # Function to check if enviroment has any test at lvl we set in init
  answer = False
  for eid, exp in cfg.items():
    for no, unit in  exp.tests.items():
      answer = answer or (unit.lvl <= lvl)
  return answer

@click.command()
@click.option('--level', default=0, help="Auto test level")
@click.option('--setup_file', required=True, help="YAML config file with list of experiments to tests.")
@click.option('--only_stderr', default=True, required=False, help="If True, only stderr will be printed on the scren.")
@click.option('--save_summary', default=True, required=False, help="If True, summary will be saved to summary.txt")

def main(setup_file,level,only_stderr,save_summary):

  with open(setup_file, "r") as ymlfile:
    full_cfg = yaml.safe_load(ymlfile)
  ymlfile.close()

  cfg = Box({**full_cfg}, default_box=True, default_box_attr=None)
  SUMMARY = {}
  
  root_dir = os.getcwd() # store entry point


  print()
  commit = subprocess.run(f"git rev-parse --short HEAD", shell=True, check=True, capture_output=True,universal_newlines=True)
  print('The repo hash for this run is:', commit.stdout)
  branch = subprocess.run(f"git rev-parse --abbrev-ref HEAD", shell=True, check=True, capture_output=True,universal_newlines=True)
  print('The branch on which it was run:', branch.stdout)


  for env, s in cfg.enviroments.items(): # loop over enviroments

    SUMMARY[env] = {}
    SUMMARY[env]['do_we_bother'] = do_we_have_tests_at_our_level(s.experiments,level)

    if SUMMARY[env]['do_we_bother']:

      print()
      print(f"::: Setting up enviroment for {s.descption or '@' + env}")
      try:
        res = subprocess.run(f"{s.setup}",
          shell=True,
          check=True,
          capture_output=True,
          universal_newlines=True)
        SUMMARY[env]['setup'] = ExecuteState.no_errors
      except subprocess.CalledProcessError as e:
        print('::: Enviroment setup errors:\n')
        print(e.stderr)
        SUMMARY[env]['setup'] = ExecuteState.with_errors

      SUMMARY[env]['TR'] = generate_TR(s.experiments) # place holder for tests results for env

      if (SUMMARY[env]['setup'] == ExecuteState.no_errors): # compile only if env setup started with no errors
   
        for eid,exp in s.experiments.items():
          os.chdir(os.path.join(root_dir,str(exp.dir)))
    
          print()
          print(f' {exp.experiment} [{exp.solver}] '.center(screen_width, ':'))
          print('::: Compilation')

          try:
            res = subprocess.run(f"make clean", shell=True, check=True, capture_output=True,universal_newlines=True)
            res = subprocess.run(f"{s.setup} make COMPILER={exp.compiler or compiler} SOLVER={exp.solver} -j 4",
              shell=True,
              check=True,
              capture_output=True,
              universal_newlines=True)
          except subprocess.CalledProcessError as e:
            print('::: Compilation Errors:')
            print(e.stderr)
            compilation_ok = False
            SUMMARY[env]['TR'][eid]["comp_res"] = ExecuteState.with_errors
          else:
            if (not only_stderr): print(res.stdout)
            if (len(res.stderr)):
              print('::: Compilation Warnings:\n')
              print(res.stderr)
              SUMMARY[env]['TR'][eid]["comp_res"] = ExecuteState.with_warnings
            else:
              SUMMARY[env]['TR'][eid]["comp_res"] = ExecuteState.no_errors
            compilation_ok = True

          if compilation_ok:

            for no, unit in  exp.tests.items():
              if (unit.lvl <= level):
                if (unit.dir):
                  os.chdir(os.path.join(root_dir, unit.dir))

                if (unit.run):
                  print(f'::: Running {unit.descption}')
                  try:
                    res = subprocess.run(f"{s.setup}{unit.run}",
                      shell=True,
                      check=True,
                      capture_output=True,
                      universal_newlines=True)
                  except subprocess.CalledProcessError as e:
                    print('::: Runtime Errors:\n')
                    print(e.stderr)
                    SUMMARY[env]['TR'][eid]["tests"][no]["result"] = ExecuteState.with_errors
                  else:
                    if (not only_stderr): print(res.stdout)
                    if (len(res.stderr)):
                      print('::: Runtime Warnings:\n')
                      print(res.stderr)
                      SUMMARY[env]['TR'][eid]["tests"][no]["result"] = ExecuteState.with_warnings
                    else:
                      SUMMARY[env]['TR'][eid]["tests"][no]["result"] = ExecuteState.no_errors

        os.chdir(root_dir)

  #Save / print summary
  
  if (save_summary):
    OF = open("summary.txt", 'w')
    def dualprint(text):
      print(text)
      OF.write(text + "\n")
  else:
    def dualprint(text):
      print(text)

  print('')
  dualprint(' SUMMARY '.center(screen_width, ':'))
  dualprint('')

  gok = 0 # flag to store a final veredict for the tests
  gno = 0

  for env, s in cfg.enviroments.items():
    env_str = f":: (env) {s.descption}"
    dualprint(f"{env_str:<65}")
    if (SUMMARY[env]['do_we_bother'] and SUMMARY[env]['setup'] == ExecuteState.no_errors):
      dualprint('')
      no = 0
      ok = 0
      for eid, e in s.experiments.items():
        exp_str = f" {no+1:02} {e.experiment} [{e.solver}] - compilation" #add compiler
        exp_res = f"{print_state(SUMMARY[env]['TR'][eid]['comp_res'])}"
        dualprint(f"{exp_str:<65}{exp_res}")
        no = no + 1
        ok = ok + test_passed(SUMMARY[env]['TR'][eid]['comp_res'])

        for tid, t in e.tests.items():
          if (t.lvl <= level and SUMMARY[env]['TR'][eid]['tests'][tid]['result'] != ExecuteState.not_tested):
            test_str = f" {no+1:02} └──── {t.descption}"
            test_res = f"{print_state(SUMMARY[env]['TR'][eid]['tests'][tid]['result'])}"
            dualprint(f"{test_str:<65}{test_res}")
            no = no + 1
            ok = ok + test_passed(SUMMARY[env]['TR'][eid]['tests'][tid]['result'])   

      final_result = f"{ok:<5}/{no:>5}"
      dualprint(('-'*11).rjust(screen_width-1))
      dualprint(final_result.rjust(screen_width-1))
      dualprint('')
      gok = gok + ok
      gno = gno + no
    else:  
      if (SUMMARY[env]['do_we_bother']):
        ok = -1
        dualprint(f":: (env) setup failed ... all tests skipped")
      else:
        dualprint(f":: (env) .. no tests for auto test lvl =< {level}")
    dualprint('')
  dualprint(':'*screen_width)

  if (save_summary):
    OF.close()

  if gno == gok:
    sys.exit(0)
  else:
    sys.exit(1)

if __name__ == '__main__':
    main()
    
