#!/usr/bin/python

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
#
# launch #j jobs of #n events on queue #q
# output name specified by #f, output dir specified by #d 
#
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

import sys
import os
import commands
from commands import getstatusoutput
from commands import getoutput
import datetime
import argparse


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

def submitJob (jobID, nEvents, queue, fileName, fileDir, execDir):
    jobname = 'jPyt_'+queue+'_'+jobID+'.sh'
    fileName += '_'+jobID+".hepmc"
    f = open (jobname, 'w')
    f.write ('#!/bin/sh' + '\n\n')
    f.write ('source '+execDir+'minbias_setup_slc6.sh \n')
    f.write (execDir+'/../genMinBias_14TeV ' + nEvents + ' ' + fileName + ' \n\n')
    f.write ('cp ' + fileName + ' ' + fileDir + '\n')
    f.close ()
    getstatusoutput ('chmod 755 ' + jobname)
    getstatusoutput ('bsub -q ' + queue + ' ' + '-u simone.pigazzini@cern.ch ' + jobname)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == '__main__':

    parser = argparse.ArgumentParser (description = 'boh')
    parser.add_argument('-j', '--jobsNum' , default = 1, help='number of sub job')
    parser.add_argument('-n', '--eventsNum' , default = 10000, help='number of events per job')
    parser.add_argument('-f', '--fileName' , help='file name')
    parser.add_argument('-d', '--copyDir' , help='copy directory')
    parser.add_argument('-q', '--queue' , default = '1nh', help='batch queue (1nh)')
    
    args = parser.parse_args ()

    execDir = getoutput('pwd')

    print 'submitting', args.jobsNum, 'jobs to queue', args.queue
    for i in range (0, int(args.jobsNum)):
        submitJob (str (i), args.eventsNum, args.queue, args.fileName, args.copyDir, execDir) 



