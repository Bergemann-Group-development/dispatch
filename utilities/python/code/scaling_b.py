"""
    Python code to add declaration and expression for scaling%b

    Command:
        find . -name scaling_mod.f90 | xargs grep -c 'self%b *=' | grep ':0$' | sed 's/:0$//' | xargs -L 1 python scaling.py
"""
import os
import sys

def pad(fd,n):
    for i in range(n):
        fd.write(' ')

file=sys.argv[1]
fo=open('xxxx.f90','w')
with open(file,'r+') as fd:
    for line in fd:
        words=line.split()

        # find a place to add declaration
        if len(words)>0:
            if words[0]=='contains':
                pad(fo,indent)
                fo.write('real(kind=8):: b\n')
        fo.write(line)
        
        # watch for 'self%e =' and append there
        if len(words) > 2:
            if words[0]=='self%e' and words[1]=="=":
                n1=line.find('self')
                pad(fo,n1)
                s='self%b = sqrt(self%d)*self%u'
                fo.write(s)

                # append a comment, aligned with previous, if exists
                n2=line.find('!')
                if n2<0:
                    n2=8
                else:
                    n2=n2-len(s)-n1
                pad(fo,n2)
                fo.write('! magnetic flux\n')

        # measure indent
        indent=0
        while(line[indent]==' '):
            indent+=1

# mv output file back
fo.close()
os.system('mv xxxx.f90 '+file)
