#!/usr/bin/env python
"""
Filter QSEQ (GERALD)/FastQ reads. Store output as FastQ.
Reads are clipped at first undetermined base (. in qseq or N in fastq)
and at first base having qual below -q.
Reads (and their pairs if -p) not passing filtering are discarded. 
Orphaned reads may be store optionally (-u).

USAGE:
filterReads.py -l31 -q10 -o outDir -p -u -v _archives/wt_10b_read*

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.es

Version 0.11

Fixes:
-0.1
--wrong sep in _clipSeq

-0.11
--output always PHRED+33 quals (Sanger, CASAVA1.8+)
"""
 
import gzip, os, sys
from datetime import datetime
from Bio import SeqIO
from optparse import OptionParser

def _clipSeq( seq,quals,minLen,sep='.' ):
  """
  """
  if sep in seq:
    pos=seq.index(sep)
    seq,quals=seq[:pos],quals[:pos]
    if len(seq)<minLen or not seq: 
      #print ' short',seq
      return
  return seq
  
def _getFastQ( file,qseq,minLen=0,qualityTh=0,ASCII_offset=64 ):
  """Process each line of GERALD output and Return FastQ str. Cut seq @ first '.' position.
  Also, check for quality if qualityTh defined.
  Return None if length of the seq < minLen or an error occured during line processing.
  """
  ##GERALD (QSEQ)
  if qseq:
    qseq=file.next()[:-1]
    
    qseq_element=qseq.split('\t') #SOLEXA 90403 4 1 23 1566 0 1 ACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCG `aaaaa```aZa^`]a``a``a]a^`a\Y^`^^]V` 1
    if len(qseq_element)!=11 or qseq_element[-1]!='1': 
      return
    
    #formatting
    name      = '@%s:%s:%s:%s:%s#%s/%s' % ( qseq_element[0],qseq_element[2],qseq_element[3],qseq_element[4],qseq_element[5],qseq_element[6],qseq_element[7] )
    seq,quals = qseq_element[8],qseq_element[9]
    
    #clip seq & quals @ . ( unknown base )
    seq       = _clipSeq( seq,quals,minLen,'.' )

  ##FASTQ
  else:
    name  = file.next()[:-1]
    seq   = file.next()[:-1]
    sep   = file.next()[:-1]
    quals = file.next()[:-1]
    
    #format name - @HWI-ST227:145:C06RAACXX:7:1101:1156:2148 2:Y:0:AACT > @HWI-ST227:145:C06RAACXX:7:1101:1156:2148/2
    if len( name.split() ) > 1:
      name = '%s/%s' % ( name.split()[0],name.split()[1][0] )
    
    #clip seq & quals @ N ( unknown base )
    seq = _clipSeq( seq,quals,minLen,'N' )
  
  if not seq:
    return
  
  #return PHRED+33 quals (Sanger encoding)
  if ASCII_offset==64:
    quals=''.join( [ chr(ord(q)-31) for q in quals ] )
  
  #cut sequence & quals @ quality
  if qualityTh:
    pos=0
    for q in quals:
      phredQ=ord(q)-33 #PHRED+33 encoding
      if phredQ<qualityTh: 
        seq,quals=seq[:pos],quals[:pos]
        if len(seq)<minLen or not seq:  
          return
        break
      pos+=1
  
  #define fastQ line
  fastq='%s\n%s\n+\n%s\n' % ( name,seq,quals )
  
  return fastq

def filterSingle( qseq,fPathF,outFileF,minLen,qualityTh,ASCII_offset,gzlipEndings=('.gz','.tgz','.tar.gz') ):
  """Convert GERALD file (fPathF) to FastQ (outFileF) with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  """
  #open input file
  i=filterSkip=0
  if fPathF.endswith( gzlipEndings ):
    fileF=gzip.open( fPathF,'rb' )
  else:
    fileF=open( fPathF,'rb' )
  #parse file
  try:
    while 1:
        i+=1
        rec=_getFastQ( fileF,qseq,minLen,qualityTh,ASCII_offset )
        
        if rec: outFileF.write(rec)
        else:   filterSkip+=1
          
  except StopIteration, e: pass
  
  return i,filterSkip,0

def filterPaired( qseq,fPathF,fPathR,outFileF,outFileR,combinedOutFile,outUnpaired,minLen,qualityTh,ASCII_offset,gzlipEndings=('.gz',) ):
  """Convert GERALD files (fPathF and fPathR) to FastQ (outFileF,outFileR) with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  Write singletons FastQ to outUnpaired if exists.
  Write filtered pairs to combinedOutFile if exists.
  """
  #open input files
  i=filterSkip=single=0
  if fPathF.endswith(gzlipEndings) and fPathR.endswith(gzlipEndings):
    fileF=gzip.open( fPathF,'rb' )
    fileR=gzip.open( fPathR,'rb' )
  else:
    fileF=open( fPathF,'rb' ) #f = gzip.open('/users/tg/lpryszcz/cluster/assembly/RNAseq/lane6_CDC317_THP_1-12_read1_qseq.txt.gz', 'rb')
    fileR=open( fPathR,'rb' )
  #parse files
  try:
    while 1:
        i+=1
        rec1=_getFastQ( fileF,qseq,minLen,qualityTh,ASCII_offset )
        rec2=_getFastQ( fileR,qseq,minLen,qualityTh,ASCII_offset )
        
        #save FastQ
        if rec1 and rec2:           #of given pair is F & R pass filtering
          outFileF.write(rec1)
          outFileR.write(rec2)
          #store combinedOut if exists
          if combinedOutFile: combinedOutFile.write(rec1+rec2)
        elif outUnpaired and rec1:  #F as single read if R didn't pass filtering
          single+=1
          outUnpaired.write( rec1 )
        elif outUnpaired and rec2:  #R as single read if F didn't pass filtering
          single+=1
          outUnpaired.write( rec2 )
        else: filterSkip+=1         #nothing ff both didn't pass filtering 
          
  except StopIteration, e: pass
  
  #close in files
  fileF.close()
  fileR.close()
  return i,filterSkip,single

def processReads( fPaths,qseq,outDir,paired,storeUnpaired,minLen,qualityTh,combinedFastQ=False,ASCII_offset_33=True,replace=False,verbose=False,qseqEnds=('.txt','.gz') ):  
  """Convert GERALD/FastQ files from Fpaths list to FastQ with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  If paired, save filtered read pairs into qXX_1.fastq and qXX_2.fastq and remove reads with missing pairs or save them as singletons into qXX.unpaired.fastq (storeUnpaired).
  If combinedFastQ, save F & R into qXX.combined.fastq.
  Print stats at the end.
  """
  #define quality score type
  if ASCII_offset_33: ASCII_offset=33 #sanger
  else:               ASCII_offset=64 #solexa/illumina
  #info
  if verbose: 
    if qseq:
      print "\nConverting gerald file(s) %s into fastq and filtering..." % fPaths, datetime.now()
    else:
      print "\nFiltering FastQ file(s) %s..." % fPaths, datetime.now()
  #process geralds
  i=filterSkip=single=0
  ###
  #process paired-ends reads
  if paired: 
    #define file names for output
    outFnameF       =os.path.join( outDir,'q%s_1.fastq'         % qualityTh )
    outFnameR       =os.path.join( outDir,'q%s_2.fastq'         % qualityTh )
    unpairedFname   =os.path.join( outDir,'q%s.unpaired.fastq'  % qualityTh )
    outCombinedFname=os.path.join( outDir,'q%s.combined.fastq'  % qualityTh )
    #check if out file exists
    if replace: pass
    elif os.path.isfile( outFnameF ) or os.path.isfile( outFnameR ):
      print " At least one of out files exist: %s or %s. Exitting." % ( outFnameF,outFnameR )
      return
    #open files for writting
    outFileF=open( outFnameF,'wb' )
    outFileR=open( outFnameR,'wb' )
    #open out file for unpaired reads
    if storeUnpaired: outUnpaired     =open( unpairedFname,'wb' )
    else:             outUnpaired     =False
    #open out file for combined FastQ
    if combinedFastQ: combinedOutFile =open( outCombinedFname,'wb' )
    else:             combinedOutFile =False
    #store pairs of filenames
    fnames_pair=[]
    #process all input files
    for fname in fPaths:
      fnames_pair.append(fname)#; print fnames_pair
      if len(fnames_pair)!=2: continue
      #get F and R qseq fnames
      fPathF,fPathR=fnames_pair 
      #proces qseq files: GERALD->FASTA
      p_i,p_filterSkip,p_single = filterPaired( qseq,fPathF,fPathR,outFileF,outFileR,combinedOutFile,outUnpaired,minLen,qualityTh,ASCII_offset )
      #print info
      if verbose: 
        print '',fnames_pair,p_i,p_filterSkip,datetime.now()
      #update read counts
      i           +=p_i
      filterSkip  +=p_filterSkip
      single      +=p_single
      #reset fnames
      fnames_pair=[]
    #close outfiles
    outFileF.close(); outFileR.close()
    #store optional out files
    if storeUnpaired: outUnpaired.close()
    if combinedFastQ: combinedOutFile.close()
    #print info for run
    print '\tProcessed pairs: %s. Filtered: %s. Reads pairs included: %s [%.2f%s]. Singletons: %s [%.2f%s]' % ( i,filterSkip,(i-filterSkip),(i-filterSkip)*100.0/i,'%',single,single*100.0/i,'%' )
    
  ####
  #for single reads
  else: 
    #define out fname
    outFnameF=os.path.join( outDir,'q%s.fastq' % ( qualityTh ) )
    #check if out file exists
    if replace: pass
    elif os.path.isfile( outFnameF ):
      print " Out file exists: %s. Exitting." % ( outFnameF, )
      return
    #open files for writting
    outFileF=open( outFnameF,'wb' )
    #process all files as single reads
    for fPathF in fPaths: 
      #proces qseq file: GERALD->FASTA
      p_i,p_filterSkip,p_single = filterSingle( qseq,fPathF,outFileF,minLen,qualityTh,ASCII_offset )
      #print info
      if verbose: print '',fPathF,p_i,p_filterSkip
      #update read counts
      i           +=p_i
      filterSkip  +=p_filterSkip
    #close outfile
    outFileF.close()
    #print info for run
    print '\tProcessed reads: %s. Filtered: %s. Reads included: %s [%.2f%s].' % ( i,filterSkip,(i-filterSkip),(i-filterSkip)*100.0/i,'%' )
  
def main():
  
  parser = OptionParser() #allow_interspersed_args=True
  
  parser.add_option("-o", dest="outDir", help="define where to store output files")
  parser.add_option("-v", action="store_true", dest="verbose", default=False,
                    help="print status messages to stdout [default: %default]")
  parser.add_option("-g", dest="qseq", action="store_true",  default=False,
                    help="input QSEQ [default: FastQ]")
  parser.add_option("-l", dest="minLen", default=31, type=int,
                    help="min read lenght (shorter after quality trimming are removed) [default: %default]" )
  parser.add_option("-q", dest="qualityTh", default=0, type=int,
                    help="read is clipped @ first base having PHRED quality lower than [default: %default]" )
  parser.add_option("-t", dest="ASCII_offset_33", action="store_false", default=True,
                    help="use illumina/solexa quality encoding (ASCII offset of 64) [default: Sanger]")
  parser.add_option("-p", action="store_true", dest="paired", default=False, 
                    help="process as paired-end reads (qXX_1.fastq & qXX_2.fastq) [default: single-end]")
  parser.add_option("-u", action="store_true", dest="unpaired", default=False, 
                    help="store orphaned reads > qXX.unpaired.fastq [default: %default]")                  
  parser.add_option("-r", dest="replace", action="store_true", default=False,
                    help="overwrite output files [default: %default]")
  parser.add_option("-c", dest="combinedFastQ", action="store_true", default=False,
                    help="store combined fastQ for paired reads as well > qXX.combined.fastq [default: %default]" )
  
  ( o, fPaths ) = parser.parse_args()
  if o.verbose: 
    print o, fPaths 
  
  if not o.outDir: 
    sys.exit( "Output dir has to be specified.")
  
  ###check input parameters
  #if any file specified
  if not fPaths: sys.exit( "Error! qseq/fastq files has to be specified as input" )
  #file input are files
  for fpath in fPaths:
    if not os.path.isfile( fpath ): 
      sys.exit( "Error! No such file: %s" % fpath )
  
  #create output directore if not present already
  if not os.path.isdir( o.outDir ): 
    os.makedirs( o.outDir )
  
  #gerald2fastq
  processReads( fPaths,o.qseq,o.outDir,o.paired,o.unpaired,o.minLen,o.qualityTh,o.combinedFastQ,o.ASCII_offset_33,o.replace,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  print "#Time elapsed: %s" % dt
