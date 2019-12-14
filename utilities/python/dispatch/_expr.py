# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 02:06:11 2019

@author: Aake
"""

import numpy as np

class Expressions(object):
    def __init__(self,patch=None,verbose=0):
        """ Initialize Exxpressions class """
        self.verbose=verbose
        self.expressions=[]
        self.ops=' +-*/()'
        self.known=[]
        self.pvars=[]
        self.patch=patch
        
        """ If a patch is given, collect all keys it accepts, except integers """
        if patch:
            for k,v in patch.keys.items():
                if k != 'numbers':
                    self.pvars.extend(v)
            if self.verbose:
                print('patch variables:',self.pvars)
        if self.verbose>1:
            print('EXPRESSIONS:')

    def split(self,rhs):
        """ Split a RHS into words, based on operators in self.ops """
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
            known patch variables, previsouly defined, or scalars
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
        p=self.patch
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

def expression_parser(patch,iv='d',verbose=0):
    E=Expressions(patch)
    E.register('iv='+iv)
    results=E.evaluate()
    return results['iv']