#!/usr/bin/env python
import os,sys

def main():
    try:
        xcov = sys.argv[1]
    except:
        sys.exit(sys.argv[0]+' [mosdepth.summary.txt file]')
    data = {}
    with open(xcov) as infile:
        while 1:
            line = infile.readline()
            if not line:
                break
            c = line.strip().split('\t')
            data[c[0]] = c[3]
    if float(data['chrY_region']) / float(data['total_region']) > 0.25:
        return '-y'
    else:
        return ''
if __name__ == '__main__':
    print(main())
