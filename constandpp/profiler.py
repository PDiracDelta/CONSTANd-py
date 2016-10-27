#!/usr/bin/python
# -*- coding: utf-8 -*-

import cProfile, pstats
from main import main

pr = cProfile.Profile()
pr.enable()
main(False, False)
pr.disable()
pr.dump_stats('main.stats')
outfilename='mainstats.txt'
with open(outfilename, 'wt') as output:
    stats = pstats.Stats('main.stats', stream=output)
    stats.sort_stats('cumulative', 'time')
    stats.print_stats()
with open(outfilename, 'r') as fin:
    print(fin.read())