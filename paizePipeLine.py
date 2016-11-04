import os,re,sys,logging,time,urllib2
import cutadaptWrapperLib as cutLib
from optparse import OptionParser,OptionGroup
from IPython import embed
from urllib2 import HTTPError,URLError

### phase 1: download reads, call fastqc and primers filtering


#############################
# Fxns for Downloading data #
#############################

def getOutDirBaseName_downloadReads(out_prefix,read_ftp):
        """ build outdir base name to be made when downloading reads  """
        outdir = os.path.join(out_prefix,read_ftp.split('/')[-2])
        return outdir

def safeMakeDir(path):
    """ a safe way to force directory creation """
    from os import makedirs as mkdir
    import errno
    try: mkdir(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path): pass
        else: raise

def downloadReads(read_ftp,out_prefix='./'):
    """ check if the url is valid and download files to `out_prefix/ACCESSION_NAME/`
        url should be in the form of http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765601/SRR765601_*.fastq.gz"""
    read1,read2,outdir = read_ftp.replace('*','1'),read_ftp.replace('*','2'),getOutDirBaseName_downloadReads(out_prefix,read_ftp)
    safeMakeDir(outdir)
    out1,out2 = map(lambda x: '%s/%s' %(outdir,os.path.basename(x)), [read1,read2])
    downloadReads_progress(read1,out1)
    downloadReads_progress(read2,out2)
    return out1,out2
    # keep the lines below for testing, they are actually NOT USED
    try: f1 = urllib2.urlopen(read1) # if read_ftp doesn't exist, raise a HTTPError
    except URLError as e:
        e.url = read1
        raise URLError
    with open(out1,'wb') as fh:
        logging.info('downloading %s to %s' %(read1,out1))
        fh.write(f1.read())
    try: f2 = urllib2.urlopen(read2)
    except URLError as e:
        e.url = read2
        raise URLError
    with open(out2,'wb') as fh:
        logging.info('downloading %s to %s' %(read2,out2))
        fh.write(f2.read())
    return out1,out2

def downloadSRA(accession,out_prefix="./",get_accession_name=False):
    """ using aspera and the sra-toolkit it will download and extract fastq files to `out_prefix/ACCESSION_NAME/`
        get_accession_name should be set as True to extract the accession name from the ftp address """
    import re
    if get_accession_name: accession = re.search("(SRR[0-9]*)_[0-9*]\.fastq\.gz",accession).group(1)
    out_dir = os.path.join(out_prefix,accession)
    os.system("bash sraDownloader.sh %s %s" %(accession,out_dir))
    out1,out2 = map(lambda x: os.path.join(out_dir,"%s_%s.fastq.gz" %(accession,x)),[1,2])
    return out1,out2

def chunk_report(bytes_so_far, chunk_size, total_size):
    percent = float(bytes_so_far) / total_size
    percent = round(percent*100, 2)
    sys.stdout.write("Downloaded %s of %s (%0.2f%%)\r" % 
        (humansize(bytes_so_far), humansize(total_size), percent))
    if bytes_so_far >= total_size: sys.stdout.write('\n')

def chunk_read(response, chunk_size=8192, report_hook=None):
    total_size = response.info().getheader('Content-Length').strip()
    total_size = int(total_size)
    bytes_so_far = 0
    data = []
    while 1:
        chunk = response.read(chunk_size)
        bytes_so_far += len(chunk)
        if not chunk: break
        data += chunk
        if report_hook: report_hook(bytes_so_far, chunk_size, total_size)
    return "".join(data)

def humansize(nbytes):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
            nbytes /= 1024.
            i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])

def downloadReads_progress(read_ftp,out_file):
    """ wrap urllib2 download with a progress bar """
    logging.info('downloading %s to %s' %(read_ftp,out_file))
    try: response = urllib2.urlopen(read_ftp)
    except URLError as e:
        e.url = read_ftp
        raise URLError
    data = chunk_read(response, report_hook=chunk_report)
    with open(out_file,'wb') as fh: fh.write(data)
    return


###############################
# Fxn for FastQC and Cutadapt #
###############################


def getFastqcOut(outdir,f):
    """ get name of fastqc output, given fastq input name """
    import os
    outdir = os.path.abspath(outdir)
    inp_prefix = os.path.basename(f).split('.fastq')[0] + '_fastqc'
    out = os.path.join(outdir,inp_prefix,'fastqc_data.txt')
    return out

def makeOutdirName(outdir_basepath,r):
    """ produce a name for the outdir name, requires the out abs path and the read name """
    import os
    root,leaf = os.path.abspath(outdir_basepath),os.path.basename(r).split('_1.fastq')[0].split('_2.fastq')[0]
    out=os.path.join(root,leaf)    
    return out

def fastqcWrapper(outdir,*reads,**optionals):
    """ fastqc wrapping fxn """
    import os
    opts = ''
    for k,v in optionals.iteritems():
        if v == '': opts += '--%s' %k
        else: opts += '--%s %s' %(k,v)
    cl = 'fastqc %s %s -o %s' %(' '.join(reads),opts,outdir)
    #embed()
    if not os.path.exists(outdir): os.mkdir(outdir)
    logging.info('Starting fastqc analysis...')
    logging.info('Command line used: %s' %cl)
    os.system(cl)
    outnames = [outdir] + map(lambda x: getFastqcOut(outdir,x),reads)
    return outnames

def getPrimers(fastqcOuts):
    """ extracts primers and other sequences to be trimmed """
    logging.info('Parsing FastQC output to identify sequences to be removed...')
    toBeRemoved = []
    for o in fastqcOuts:
        logging.info('Parsing %s...' %o)
        trash = [t for t in cutLib.trashExtractor(o)]
        toBeRemoved.append(cutLib.getAdaptersAndPolyAs(trash))
    if len(toBeRemoved) == 1: return toBeRemoved
    adapters1,polyA1,universal1 = toBeRemoved[0]
    adapters2,polyA2,universal2 = toBeRemoved[1]
    adapters = adapters1.union(adapters2)
    logging.info('The following adapters were found:\n\t%s' %'\n\t'.join(list(adapters)))
    polyAs = polyA1 + polyA2
    logging.info('The following polyNs were found:\n\t%s' %'\n\t'.join(polyAs))
    universal = [i for i in [universal1,universal2] if i != None]
    if universal == []: universal = None
    else: universal = universal[0]
    primers = [adapters,polyAs,universal]
    return primers

def getPrimers_(fastqcOuts):
    """ extracts primers and other sequences to be trimmed (old) """
    toBeRemoved = []
    for o in fastqcOuts:
        trash = [t for t in cutLib.trashExtractor(o)]
        toBeRemoved.append(cutLib.getAdaptersAndPolyAs_old(trash)) # only difference between old and new fxn version
    if len(toBeRemoved) == 1: return toBeRemoved
    adapters1,polyA1,universal1 = toBeRemoved[0]
    adapters2,polyA2,universal2 = toBeRemoved[1]
    adapters = adapters1.union(adapters2)
    polyAs = polyA1 + polyA2
    universals = [i for i in [universal1,universal2] if i != None]
    if universals == []: universal = None
    else: universal = universals[0]
    primers = [adapters,polyAs,universal]
    return primers

def cutadaptWrapper(reads,primers,options):
    """ cutadapt wrapping fxn """
    adapters,universal = list(primers[0]) + primers[1],primers[-1]
    cutadaptCl,trimmed1,trimmed2 = cutLib.makeCutAdaptCL(reads[0],reads[1],universal,adapters)
    # embed()
    os.system(cutadaptCl)
    return [trimmed1,trimmed2]

def presetOptions(fastqcOpts=None,cutadaptOpts=None):
    """ should be used for our batch analysis """
    from multiprocessing import cpu_count
    presetFastQc = {'extract':' ','threads':cpu_count()/2}
    presetCutAdapt = {}
    if fastqcOpts == None: fastqcOptions = presetFastQc
    else: fastqcOptions = fastqcOpts
    if cutadaptOpts == None: cutadaptOptions = presetCutAdapt
    else: cutadaptOptions = cutadaptOpts
    return [fastqcOptions,cutadaptOptions]

def phase1(reads,fastqcOptions,cutadaptOptions,preset=False,outdir_basepath='.'):
    """ wrapper for phase 1 (reads processing) """
    if preset: fastqcOptions,cutadaptOptions = presetOptions()
    # fastqc
    outdir = makeOutdirName(outdir_basepath,reads[0])
    fastqcDir,fastqcOut1,fastqcOut2 = fastqcWrapper(outdir,*reads,**fastqcOptions)
    # extract primers
    # embed()
    primers = getPrimers([fastqcOut1,fastqcOut2])
    # cutadapt
    trimmedReads = cutadaptWrapper(reads,primers,cutadaptOptions)
    return trimmedReads


##################################
# Fxn for Cleaning and Profiling #
##################################


def cleanFiles(dir_list):
    """     to be called when an exception is raised or when a single iteration is over.
        Delete the input files to save disk space. """
    # TODO test and debug this
    for d in dir_list:
        if not os.path.exists(d): continue
        for f in os.listdir(d):
                f_ = os.path.join(d,f)
                if os.path.isfile(f_): os.remove(f_)
                else: cleanFiles([f_])
        os.rmdir(d)
    return

def moveFiles(accession,important_files):
        """ copy output files to SSUP server. A public key has to be set for this function to work! """
        import hashlib,os
        try: os.system('ssh -t edoardo@193.205.83.138 "mkdir -p ~/progetti_completi/Remote_Results_Paize/%s"' %accession)
        except:
                # check if dir exists
                try: os.system('ssh -t edoardo@193.205.83.138 "ls progetti_completi/Remote_Results_Paize/%s"' %accession)
                except Exception as e: raise e
                # if this is passed, continue (the directory already exists)
        for f in important_files:
                retries = 0
                while True:
                        remoteDestination = 'progetti_completi/Remote_Results_Paize/%s/%s'  %(accession,os.path.basename(f))
                        md5Start = hashlib.md5(open(f,'rb').read()).hexdigest()
                        os.system('scp %s edoardo@193.205.83.138:~/%s' %(f,remoteDestination))
                        md5End = os.popen('ssh -t edoardo@193.205.83.138 "md5sum %s"' %remoteDestination).read().split()[0]
                        if md5Start == md5End:
                                logging.info('successfully transfered %s!' %f)
                                break
                        logging.warning('failed transfering %s! Retrying again for 3 times' %f)
                        retries +=1
                        if retries >= 2:
                                logging.error('failed to transfer %s! Max tries reached' %f)
                                raise Exception('failed to transfer %s! Max tries reached' %f)
        return

def moveFiles_local(destination,accession,important_files):
        """ copy output files to a local destination """
        import os
        destination_dir = os.path.join(destination,accession)
        safeMakeDir(destination_dir)
        for f in important_files:
                new_destination = os.path.join(destination_dir,os.path.basename(f))
                os.rename(f,new_destination)
        return


def profilingPhase1(accession,times):
    time_start,time_download,time_trimming,time_final = times
    date_start,date_end = map(lambda x: time.asctime(time.localtime(x)), [time_start,time_final])
    time_to_download = time_download - time_start 
    time_to_trim = time_trimming - time_download
    time_to_tophat = time_final - time_trimming
    total_time = time_final - time_start
    s = """ Pipeline profiling for %s accession:
        started at %s \t ended at %s
        time to download reads: %s seconds (%s minutes)
        time to trim: %s (%s minutes)
        time to run Tophat2: %s (%s minutes)
        """     %(accession,date_start,date_end,
            round(time_to_download,2),round(time_to_download/60.,2),
            round(time_to_trim,2),round(time_to_trim/60.,2),
            round(time_to_tophat,2),round(time_to_tophat/60,2))
    return s

##################################
# Tophat and stuff               #
##################################


def tophat(reads,outdir='prova',reference='reference/Zea_mays.AGPv3.31.dna.genome',junctions=None,threads=1):
    """ TopHat wrapper """
    import os
    reads1,reads2 = reads
    if junctions == None: cl='tophat-2.1.1.Linux_x86_64/tophat2 -o %s -p %s -i 10 -I 60000 --library-type fr-unstranded %s %s %s' %(outdir,threads,reference,reads1,reads2)
    else: cl='tophat-2.1.1.Linux_x86_64/tophat2 -o %s -p %s -i 10 -I 600000 -G -j %s --no-novel-juncs --library-type fr-unstranded %s %s %s' %(outdir,threads,junctions,reference,reads1,reads2)
    logging.info('command line used for tophat: %s' %cl)
    os.system(cl)
    return map(lambda x: '%s/%s' %(outdir,x) ,['junctions.bed','align_summary.txt']) + list(reads) #tophat

def setJunctions(junction_file):
    """ Produce junctions for the second run of TopHat """
    from os.path import dirname,abspath,join
    from os import system
    junct_dir = dirname(abspath(junction_file))
    all_juncts = join(junct_dir,'all_junctions.bed')
    final_juncts = join(junct_dir,'all_junctions.juncts')
    cl1 = 'cat %s > %s' %(junction_file,all_juncts)
    cl2 = 'bed_to_juncs < %s | sort -k 1,4 -u | sort -k 1,1 > %s' %(all_juncts,final_juncts)
    system(cl1)
    system(cl2)

def phase2(reads,outdir):
    """ wrapper for phase 2 (tophat) """
    from os.path import join,abspath
    # first tophat run
    tophat(reads,outdir)
    # collect junctions
    junctions = join(abspath(outdir),'junctions.bed')
    # second tophat run
    tophat(reads,outdir,junctions)

##################################
# Unit testing                   #
##################################

def test_download(sra=True):
    download_reads = downloadReads
    if sra: download_reads = lambda *args,**kwargs: downloadSRA(*args,get_accession_name=True,**kwargs)
    i, reads = 'http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765378/SRR765378_*.fastq.gz',['Accessions/SRR765378/SRR765378_1.fastq.gz','Accessions/SRR765378/SRR765378_2.fastq.gz']
    accession_name = i.split('/')[-2]
    logging.info('Working on %s: url is %s' %(accession_name,i))
    out_prefix='Accessions_testing'
    outDir_download = getOutDirBaseName_downloadReads(out_prefix,i)
    print outDir_download
    try: reads = download_reads(i,out_prefix=out_prefix)
    except (HTTPError,URLError) as e: 
        logging.error('''%s for accession %s when downloading %s
                    Skipping downstream analyses for this accession...''' %(e.__repr__(),accession_name,e.url))
        print "test failed!"
    else: print "test passed!"
    
def test_trimming():
	return

def test_tophat():
	return

def test_moving(dir_= "toSend_testing",accession_name="SRR765379",local=False):
	return




##################################
# Main                           #
##################################

def main(options):
    """ wrapper for main. Quick and dirty.
    TODO:     1. add a way to restart from a stopped iteration 
        2. redirect stdout e stderr from os.system calls to the logger"""
    # 0) instanciate logger and console handler

    # set options
    urls_file,reference,threads_ = options.accession_file,options.reference,options.threads
    # if not os.path.isfile(reference): raise OSError("reference file %s was not found" %reference)
    move_files,download_reads = moveFiles,downloadReads
    if options.localOutDir: move_files = lambda *args: moveFiles_local(options.localOutDir,*args)
    if options.downloader: download_reads = lambda *args,**kwargs: downloadSRA(*args,get_accession_name=True,**kwargs)

    # set logger
    failed,successful = [],[]
    logger_name = 'paizeLogger_' + time.strftime("%d-%m-%Y_%I-%M-%S") + '.log'
    logging_format = '%(asctime)-15s %(message)s'
    logging.basicConfig(filename=logger_name,format=logging_format,level=logging.DEBUG)
    logging.getLogger().addHandler(logging.StreamHandler())

    # 1) collect initial inputs (list of directory with data for accessions)
    inputs = [i.strip() for i in open(urls_file)]
    logging.info('There are %s accessions: \n\t%s' %(len(inputs),'\n\t'.join(map(lambda x: x.split('/')[-2],inputs))))

    # 2) start iterating
    logging.info('Starting iteration:')
    for i in inputs:
        # i, reads = 'http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765378/SRR765378_*.fastq.gz',['Accessions/SRR765378/SRR765378_1.fastq.gz','Accessions/SRR765378/SRR765378_2.fastq.gz'] # debug only, delete after
        toSend,toDelete = [],[]
        time_start = time.time()
        accession_name = i.split('/')[-2]
        logging.info('Working on %s: url is %s' %(accession_name,i))
        out_prefix='Accessions'
        outDir_download = getOutDirBaseName_downloadReads(out_prefix,i)
        print outDir_download
        try: reads = download_reads(i,out_prefix='Accessions')
        except (HTTPError,URLError) as e: 
            logging.error('''%s for accession %s when downloading %s
                    Skipping downstream analyses for this accession...''' %(e.__repr__(),accession_name,e.url))
            failed.append(accession_name)
            if os.path.exists(outDir_download): cleanFiles([outDir_download])
            continue
        else: toDelete.append(outDir_download)
        time_download = time.time()
        # 2.1) go with fastq up to first tophat
        #### fastqc and trimming
        logging.info('Performing trimming')
        try: trimmed_reads = cutLib.trimmomaticWrapper(*reads)
        except Exception as e:
            logging.error('''%s
                    Skipping downstream analyses for this accession
                    Moreover, the input files will be deleted!''' %e)
            cleanFiles(toDelete)
            failed.append(accession_name)
            continue
        time_trimming = time.time()
        #### tophat
        #trimmed_reads = ["Accessions/SRR765378/trimmed_paired_SRR765378_1.fastq.gz","Accessions/SRR765378/trimmed_paired_SRR765378_2.fastq.gz"]
        logging.info('Performing tophat2 on the trimmed reads:\n%s' %('\t\n'.join(trimmed_reads)))
        # embed()
        tophat_out = 'tophat_out_' + accession_name
        try: important_outfiles = tophat(trimmed_reads,reference=reference,outdir=tophat_out,threads=threads_)
        except Exception as e:
            logging.error('''%s
                    Skipping downstream analyses for this accession
                    Moreover, the input files will be deleted!
                    The command line used for tophat was ''' %e)
            cleanFiles(toDelete)
            failed.append(accession_name)
            continue
        else:
            map(lambda x: toSend.append(x),important_outfiles)
            toDelete.append(tophat_out)
        time_final = time.time()
        logging.info('Analyses finished for accession %s!' %accession_name)
        logging.info(profilingPhase1(accession_name,[time_start,time_download,time_trimming,time_final]))
        # 2.2) copy trimmed reads to SSUP server, then delete everything but junction and stats in tophat dir
        logging.info('Copying trimmed reads to SSUP server, then deleting everything but junction and stats file...')
        # break #REMOVE THIS WHEN FINISHED TODO
        try: move_files(accession_name,toSend)
        except Exception as e:
            failed.append(accession_name)
            logging.error('''%s 
                     Failed at copying file for accession %s''' %(e,accession_name))
        else:
            logging.info('Successful iteration for accession %s!'%accession_name)
            successful.append(accession_name)
        cleanFiles(toDelete)
    # 3) end the pipeline, report number of successful and failed iterations
    # return #REMOVE THIS WHEN FINISHED TODO
    logging.info('''paizePipeline ended!
            successful iterations for %s accessions:\n%s
            --------------------------------------------
            failed iterations for %s accessions:\n%s
            ''' %(len(successful),'\n'.join(successful),len(failed),'\n'.join(failed))
            )
# test

#reads=['accession_data/aSRR765127/SRR765127_1.fastq.gz','accession_data/aSRR765127/SRR765127_2.fastq.gz']
#reads=['accession_data/SRR765128/SRR765128_1.fastq.gz','accession_data/SRR765128/SRR765128_2.fastq.gz']

#try_phase1 = phase1(reads,None,None,preset=True,outdir_basepath='.')

#trimmed1,trimmed2 = try_phase1
#tophat_outdir = 'prova_tophat'

if __name__ == "__main__":

        ###############
        # Opt Parsing #
        ###############

    usage=""" %prog [options]
        """

    parser = OptionParser(usage=usage)

    # mandatory
    group1 = OptionGroup(parser, "Mandatory Arguments")
    group1.add_option("-i", "--input", dest="accession_file",
                      help="target genome to be scaffolded", metavar="FILE")
    group1.add_option("-r","--reference",dest="reference",
                      help="genome to be used as reference, must be preindexed", metavar="FILE")
    parser.add_option_group(group1)

    # optional
    group2 = OptionGroup(parser, "Optional Arguments")
    group2.add_option("-o", "--localoutput", dest="localOutDir",
                      help="if this option is set, outputs are moved to DIR instead than to SSUP", metavar="DIR")
    group2.add_option("-t", "--threads", dest="threads",default=1,
                      help="set the number of threads. It affects the overall speed")
    group2.add_option("-s", "--sra-toolkit", dest="downloader", action="store_true",
                      help="""if this option is set, fastq files are downloaded using sra-toolkit, sensibly increasing the speed.
                                Requires the installation of ```sra-toolkit``` and ```aspera connect``` working.
                                It might be necessary to tweak the options present in ```sraDownloader.sh``` """,default=None)
    parser.add_option_group(group2)
    (options, args) = parser.parse_args()


    if not (options.accession_file and options.reference):
        parser.print_help()
        parser.error('Mandatory Arguments missing')
    main(options)

