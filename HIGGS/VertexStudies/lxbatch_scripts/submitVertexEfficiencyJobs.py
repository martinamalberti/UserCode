import os
import sys
import numpy
import glob

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i","--inputdir",dest="inputdir",type="string",help="Path where ntuples are located. Example: /store/cmst3/user/malberti/HIGGS/VERTEX/2012/DATA/")
parser.add_option("-c","--config",dest="config",type="string",help="Configuration file template. Example: TMVA_check_ZmumuData_TEMPLATE.cfg")
parser.add_option("-w","--workdir",dest="workdir",type="string",default="EffJobs",help="Name of the directory for jobs")
parser.add_option("-o","--outputname",dest="outputname",type="string",default="testEfficiency",help="Name of the output file")
parser.add_option("-n","--njobs",dest="njobs",type="int",help="Number of jobs")
parser.add_option("-q","--queue",dest="queue",type="string",default="1nh",help="Name of the queue on lxbatch")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true",help="Do not submit jobs")
parser.add_option("","--checkJobs",dest="checkJobs",action="store_true",default=False,help="Checks job status")
parser.add_option("","--resubmit",dest="resubmit",action="store_true",default=False,help="Resubmit job ")
parser.add_option("-j","--jobs",dest="jobs",type="int",help="Job number (for resubmission). You can resubmit one job at time for now.")

(options,args)=parser.parse_args()

eos = '/afs/cern.ch/project/eos/installation/pro/bin/eos.select'

#-----------------------------------------
# prepare the list of files to be analyzed
#-----------------------------------------
path = os.getcwd()
wdir = path+'/'+options.workdir

if not options.checkJobs and not options.resubmit:
    os.mkdir(wdir)
    listoffiles = []
    command = ('%s find -f %s | grep root > %s/list.txt' % (eos,options.inputdir,wdir))
    print command
    os.system(command)
    file = open('%s/list.txt'%wdir, 'r')
    listoffiles=[line.replace('/eos/cms/','root://eoscms//eos/cms/').replace('\n','') for line in file ]
    print 'Found %d files' %len(listoffiles)
    #print listoffiles
    
    #-----------------------------------------
    # now split the jobs
    #-----------------------------------------
    for job in range(options.njobs):
        print 'job %d' %job
        jobdir = '%s/JOB_%d'%(wdir,job)
        os.mkdir(jobdir)
        
        #--- prepare the list of files for each job 
        f = open('%s/input_%d.txt'%(jobdir,job), 'w')
        sublist = [file for i,file in enumerate(listoffiles) if (i%options.njobs==job)]
        for fname in sublist:
            f.write('%s '%fname) 
            
        #--- prepare the configuration file
        #print ('cp %s %s/%d.cfg' %(options.config,jobdir,job))
        os.system('cp %s %s/temp_%d.cfg' %(options.config,jobdir,job))
        cfgtemp = open('%s/temp_%d.cfg'%(jobdir,job),'r')
        cfg     = open('%s/%d.cfg'%(jobdir,job),'w')
        for line in cfgtemp:
            cfg.write(line.replace('INPUTLIST','%s/input_%d.txt'%(jobdir,job)).replace('OUTPUTDIR','%s/'%jobdir).replace('OUTPUTFILENAME','%s_%d.root'%(options.outputname,job)) )
        os.system('rm %s'%cfgtemp.name)
            
        #--- prepare the jobs scripts
        jobscript = open('%s/bjob_%d.sh'%(jobdir,job),'w')
        jobscript.write('cd %s \n'%jobdir)
        jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc462 \n')
        jobscript.write('eval ` scramv1 runtime -sh ` \n')
        jobscript.write('source ../../../scripts/setup.sh \n')
        jobscript.write('if ( \n')
        jobscript.write('\t touch %s/bjob_%d.run \n'%(jobdir,job))
        jobscript.write('\t MyVertexAnalysis.exe %s/%d.cfg \n'%(jobdir,job))
        jobscript.write(') then \n')
        jobscript.write('\t mv ./%s_%d.root %s \n'%(options.outputname,job,jobdir))
        jobscript.write('\t touch %s/bjob_%d.done \n'%(jobdir,job))
        jobscript.write('else \n')
        jobscript.write('\t touch %s/bjob_%d.fail \n'%(jobdir,job))
        jobscript.write('fi \n')
        os.system('chmod a+x %s/bjob_%d.sh'%(jobdir,job))
        #print 'bsub -q %s -o %s/%s_%d.log %s'%(options.queue,jobdir,options.outputname,job,jobscript.name )
        if not options.dryRun:
            os.system('bsub -q %s -o %s/%s_%d.log %s'%(options.queue,jobdir,options.outputname,job,jobscript.name ))
        
elif options.resubmit and options.jobs >-1 :
    print 'Resubmitting job %d ' %options.jobs
    resubcmd = 'bsub -q %s -o %s/JOB_%d/%s_%d.log %s/JOB_%d/bjob_%d.sh'%(options.queue,wdir,options.jobs,options.outputname,options.jobs,wdir,options.jobs,options.jobs )
    #print resubcmd
    os.system(resubcmd)

elif options.checkJobs:
    jobs = glob.glob( '%s/JOB_*/bjob*.sh'% (wdir) )
    print 'Total number of jobs: %d' %len(jobs)

    listdone = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/bjob_%d.done' % (wdir,j,j))]
    print 'Total number of DONE jobs: %s ' % len(listdone)
    print '  %s' %listdone
    for j in listdone:
        f = '%s/JOB_%d/bjob_%d.run'%(wdir,j,j)
        #print 'rm %s/JOB_%d/bjob_%d.run'%(wdir,j,j)
        if (os.path.isfile(f)):
            os.system('rm %s'%f)
        
    listrun = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/bjob_%d.run' % (wdir,j,j))]
    print 'Total number of RUNNING jobs: %d ' %len(listrun)
    print '   %s' %listrun
    
    listfailed = [j for j in range(len(jobs)) if os.path.isfile('%s/JOB_%d/bjob_%d.fail' % (wdir,j,j))]
    print 'Failed jobs: %s ' % listfailed
    print '   %s' %listfailed



    if (len(listdone) == len(jobs)):
        print "All jobs successful! Merging output files..."
        os.chdir(wdir)
#        os.system('ls JOB*/*.root')
        os.system('hadd -f testEfficiency.root JOB_*/*.root')
    else:    
        for j in listfailed:
            print 'bsub -q %s -o %s/JOB_%d/%s_%d.log %s/JOB_%d/bsub_%d.sh'%(options.queue,wdir,j,options.outputname,j,wdir,j,j )
    
