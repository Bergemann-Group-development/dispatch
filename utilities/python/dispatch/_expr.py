# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 02:06:11 2019

@author: Aake
"""

import numpy as np

class Expressions(object):
    def __init__(self,task=None,verbose=0):
        """ Initialize Exxpressions class """
        self.verbose=verbose
        self.expressions=[]
        self.ops=' +-*/()[]:,'
        self.known=[]
        self.pvars=[]
        self.task=task
        
        """ If a task is given, collect all keys it accepts, except integers.
            A task is assume to have an attribute 'keys', which is a dictionry
            with categories of known variable names.  This names are added to
            the 'pvars' list of variable names.
        """
        if task:
            for k,v in task.keys.items():
                if k != 'numbers':
                    self.pvars.extend(v)
            if self.verbose:
                print('task variables:',self.pvars)
        if self.verbose>1:
            print('EXPRESSIONS:')

    def split(self,rhs):
        """ Split a RHS into words, based on operators in self.ops, returning
            a list of words that are potentially names of known variables.
        """
        words=[]
        def add(w):
            if w != '' and w not in words:
                words.append(w) 
        word=''
        for c in rhs:
            if c in self.ops:
                add(word)
                word=''
            else:
                word+=c
        add(word)
        return words

    def register(self,expression):
        """ Register expressions, checking each RHS for words that are either
            known task variables, previsouly defined, or scalars
        """
        if self.verbose>1:
            print(' ',expression)
        lhs,rhs=expression.split('=')
        self.expressions.append((lhs,rhs))
        for word in self.split(rhs):
            if word in self.known or word in self.pvars:
                continue
            else:
                if word in globals():
                    print('WARNING: {} in globals()'.format(word))
                try:
                    eval(word)
                    self.known.append(word)
                except:
                    print('ERROR: could not evaluate',word)
                    return False
        self.known.append(lhs)
        return True
        
    def evaluate(self):
        """ After all expressions registered, evaluate """
        if self.verbose:
            print('EVALUATE:')
        ll=locals()
        results={}
        p=self.task
        for lhs,rhs in self.expressions:
            for word in self.split(rhs):
                if not (word in self.known or word in ll):
                    if word in self.pvars:
                        cmd="p.var('{}')".format(word)
                        try:
                                ll[word]=eval(cmd)
                        except:
                            print('ERROR: {}'.format(cmd))
                            return False
                    else:
                        print('ERROR: {} not known'.format(word))
                    
            """ Now the whole RHS should be possible to evaluate """
            if self.verbose:
                print('  {}={}'.format(lhs,rhs))
            ll[lhs]=eval(rhs)
            results[lhs]=ll[lhs]
            self.known.append(lhs)
        return results

def expression_parser(task,iv='d',verbose=0):
    """ The expression parser takes an object that has a dictionary with
        categories of known variables (attributes), registers an expression
        of the type 'iv=expression', and then returns the value of the LHS
        if that expression (i.e. the value of 'iv'), by splitting the expression
        into words, and evaluating into the current context the values of all
        known variables.   If all the words that occur in the expression are
        known, then the expression will evaluate, and the result can be returned.
    """
    E=Expressions(task)
    E.register('iv='+iv)
    results=E.evaluate()
    return results['iv']