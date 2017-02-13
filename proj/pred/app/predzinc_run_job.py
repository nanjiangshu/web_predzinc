#!/usr/bin/env python
# Description: run job

# ChangeLog 
#
# ChangeLog 2015-02-12 
#   submit individual sequences to the workflow, so that the result of each
#   sequence can be cached and the progress can be shown for a job with many
#   sequences
# ChangeLog 2015-09-25
#   result from cache just make a soft link, 
#   zip -rq will replace the symbolic link with the actual data when making the
#   zip file

# how to create md5
# import hashlib
# md5_key = hashlib.md5(string).hexdigest()
# subfolder = md5_key[:2]

#

import os
import sys
import subprocess
import time
import myfunc
import glob
import hashlib
import shutil
import datetime
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
rundir = os.path.dirname(os.path.realpath(__file__))

ZB_SCORE_THRESHOLD = 0.45
chde_table = {
        'C':'CYS',
        'H': 'HIS',
        'D': 'ASP',
        'E': 'GLU',
        'CYS': 'C',
        'HIS': 'H',
        'ASP': 'D',
        'GLU': 'E'
        }

blastdir = "%s/%s"%(rundir, "soft/predzinc")
os.environ['BLASTMAT'] = "%s/data"%(blastdir)
os.environ['BLASTBIN'] = "%s/bin"%(blastdir)
os.environ['BLASTDB'] = "/data/blastdb"
#blastdb = "%s/%s"%(os.environ['BLASTDB'], "uniref90.fasta" )
blastdb = "%s/%s"%(os.environ['BLASTDB'], "uniref90.fasta" )
#blastdb = "/data3/data/blastdb/swissprot"
predzinc_db = "passe"
runscript = "%s/%s"%(rundir, "soft/predzinc/predzinc.sh")

basedir = os.path.realpath("%s/.."%(rundir)) # path of the application, i.e. pred/
path_md5cache = "%s/static/md5"%(basedir)

contact_email = "nanjiang.shu@scilifelab.se"
vip_user_list = [
        "nanjiang.shu@scilifelab.se"
        ]

# note that here the url should be without http://

usage_short="""
Usage: %s seqfile_in_fasta 
       %s -jobid JOBID -outpath DIR -tmpdir DIR
       %s -email EMAIL -baseurl BASE_WWW_URL
       %s [-force]
"""%(progname, wspace, wspace, wspace)

usage_ext="""\
Description:
    run job

OPTIONS:
  -force        Do not use cahced result
  -h, --help    Print this help message and exit

Created 2015-02-05, updated 2015-02-12, Nanjiang Shu
"""
usage_exp="""
Examples:
    %s /data3/tmp/tmp_dkgSD/query.fa -outpath /data3/result/rst_mXLDGD -tmpdir /data3/tmp/tmp_dkgSD
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def CountNumPredZB(predfile, threshold=0.45):#{{{
    cntZB = 0
    cntHomo = 0
    hdl = myfunc.ReadLineByBlock(predfile)
    if not hdl.failure:
        lines = hdl.readlines()
        while lines != None:
            for line in lines:
                if not line:
                    continue
                if line[0] == "#":
                    strs = line.split()
                    if len(strs) >=6 and strs[0] == "#Homolog":
                        cntHomo += 1
                else:
                    strs = line.split("\t")
                    try:
                        score = float(strs[2])
                        if score >= threshold:
                            cntZB += 1
                    except:
                        pass
            lines = hdl.readlines()
        hdl.close()
        return (cntZB, cntHomo)
    else:
        return (cntZB, cntHomo)
#}}}
def WriteNiceResult(predfile, fpout, threshold=0.45):#{{{
    is_ZB = False
    is_has_homo = False
    li_homolog = []
    li_ZB = []
    li_other = []
    hdl = myfunc.ReadLineByBlock(predfile)
    if not hdl.failure:
        lines = hdl.readlines()
        while lines != None:
            for line in lines:
                if not line:
                    continue
                if line[0] == "#":
                    strs = line.split()
                    if len(strs) >=6 and strs[0] == "#Homolog":
                        id_homo = strs[4]
                        score_homo = strs[5]
                        li_homolog.append([id_homo, score_homo])
                else:
                    strs = line.split("\t")
                    if len(strs)>=5:
                        key = strs[0]
                        ss2 = key.split(";")
                        if len(ss2)>=5:
                            res = ss2[3]
                            series = ss2[4]
                            try:
                                score = float(strs[2])
                            except:
                                score = -100.0
                            if score >= threshold:
                                li_ZB.append([res, series, score])
                            else:
                                li_other.append([res, series, score])
            lines = hdl.readlines()
        hdl.close()
        if len(li_ZB) > 0:
            is_ZB = True
        if len(li_homolog) > 0:
            is_has_homo = True

        try:
            if len(li_ZB) > 0:
                fpout.write("The following %d residues were predicted as zinc-binding (with score >= %g, sorted by scores\n\n"%(len(li_ZB), threshold))
                fpout.write("%-3s %8s %6s\n"%("Res","SerialNo","Score"))
                li_ZB = sorted(li_ZB, key=lambda x:x[2], reverse=True)
                for item in li_ZB:
                    fpout.write("%-3s %8s %6.3f\n"%(chde_table[item[0]],
                        item[1], item[2]))
            else:
                fpout.write("No residues were predicted as zinc-binding\n\n")

            if len(li_other) >0:
                fpout.write("\n\nPrediction scores for the rest %d CHDEs, sorted by scores\n\n"%(
                    len(li_other)))
                fpout.write("%-3s %8s %6s\n"%("Res","SerialNo","Score"))
                li_other = sorted(li_other, key=lambda x:x[2], reverse=True)
                for item in li_other:
                    fpout.write("%-3s %8s %6.3f\n"%(chde_table[item[0]],
                        item[1], item[2]))

            fpout.write("//\n\n") # write finishing tag
        except:
            pass

    else:
        pass

    return (is_ZB, is_has_homo)
#}}}
def WriteTextResultFile(outfile, outpath_result, maplist, runtime_in_sec, statfile=""):#{{{
    try:
        fpout = open(outfile, "w")

        fpstat = None
        num_ZB = 0
        num_has_homo = 0

        if statfile != "":
            fpstat = open(statfile, "w")

        date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print >> fpout, "##############################################################################"
        print >> fpout, "PredZinc result file"
        print >> fpout, "Generated from %s at %s"%(g_params['base_www_url'], date)
        print >> fpout, "Total request time: %.1f seconds."%(runtime_in_sec)
        print >> fpout, "##############################################################################"
        cnt = 0
        for line in maplist:
            strs = line.split('\t')
            subfoldername = strs[0]
            length = int(strs[1])
            desp = strs[2]
            seq = strs[3]
            print >> fpout, "Sequence number: %d"%(cnt+1)
            print >> fpout, "Sequence name: %s"%(desp)
            print >> fpout, "Sequence length: %d aa."%(length)
            print >> fpout, "Sequence:\n%s\n\n"%(seq)

            is_ZB = False
            is_has_homo = False
            outpath_this_seq = "%s/%s"%(outpath_result, subfoldername)
            predfile = "%s/query.predzinc.predict"%(outpath_this_seq)
            g_params['runjob_log'].append("predfile =  %s.\n"%(predfile))
            if not os.path.exists(predfile):
                g_params['runjob_log'].append("predfile %s does not exist\n"%(predfile))
            (is_ZB, is_has_homo ) = WriteNiceResult(predfile, fpout,
                    threshold=ZB_SCORE_THRESHOLD) 

            if fpstat:
                num_ZB += is_ZB
                num_has_homo += is_has_homo

            cnt += 1

        if fpstat:
            out_str_list = []
            out_str_list.append("num_ZB %d"% num_ZB)
            out_str_list.append("num_has_homo %d"% num_has_homo)
            fpstat.write("%s"%("\n".join(out_str_list)))
            fpstat.close()
    except IOError:
        print "Failed to write to file %s"%(outfile)
#}}}
def RunJob(infile, outpath, tmpdir, email, jobid, g_params):#{{{
    all_begin_time = time.time()

    rootname = os.path.basename(os.path.splitext(infile)[0])
    starttagfile   = "%s/runjob.start"%(outpath)
    runjob_errfile = "%s/runjob.err"%(outpath)
    runjob_logfile = "%s/runjob.log"%(outpath)
    finishtagfile = "%s/runjob.finish"%(outpath)
    rmsg = ""


    resultpathname = jobid

    outpath_result = "%s/%s"%(outpath, resultpathname)
    tarball = "%s.tar.gz"%(resultpathname)
    zipfile = "%s.zip"%(resultpathname)
    tarball_fullpath = "%s.tar.gz"%(outpath_result)
    zipfile_fullpath = "%s.zip"%(outpath_result)
    outfile = "%s/%s/Topcons/topcons.top"%(outpath_result, "seq_%d"%(0))
    resultfile_text = "%s/%s"%(outpath_result, "query.result.txt")
    mapfile = "%s/seqid_index_map.txt"%(outpath_result)
    finished_seq_file = "%s/finished_seqs.txt"%(outpath_result)



    tmp_outpath_result = "%s/%s"%(tmpdir, resultpathname)
    isOK = True
    try:
        os.makedirs(tmp_outpath_result)
        isOK = True
    except OSError:
        msg = "Failed to create folder %s"%(tmp_outpath_result)
        myfunc.WriteFile(msg+"\n", runjob_errfile, "a")
        isOK = False
        pass

    try:
        os.makedirs(outpath_result)
        isOK = True
    except OSError:
        msg = "Failed to create folder %s"%(outpath_result)
        myfunc.WriteFile(msg+"\n", runjob_errfile, "a")
        isOK = False
        pass


    if isOK:
        try:
            open(finished_seq_file, 'w').close()
        except:
            pass
#first getting result from caches
# ==================================

        maplist = []
        maplist_simple = []
        toRunDict = {}
        hdl = myfunc.ReadFastaByBlock(infile, method_seqid=0, method_seq=0)
        if hdl.failure:
            isOK = False
        else:
            datetime = time.strftime("%Y-%m-%d %H:%M:%S")
            rt_msg = myfunc.WriteFile(datetime, starttagfile)

            recordList = hdl.readseq()
            cnt = 0
            origpath = os.getcwd()
            while recordList != None:
                for rd in recordList:
                    isSkip = False
                    # temp outpath for the sequence is always seq_0, and I feed
                    # only one seq a time to the workflow
                    tmp_outpath_this_seq = "%s/%s"%(tmp_outpath_result, "seq_%d"%0)
                    outpath_this_seq = "%s/%s"%(outpath_result, "seq_%d"%cnt)
                    subfoldername_this_seq = "seq_%d"%(cnt)
                    if os.path.exists(tmp_outpath_this_seq):
                        try:
                            shutil.rmtree(tmp_outpath_this_seq)
                        except OSError:
                            pass

                    maplist.append("%s\t%d\t%s\t%s"%("seq_%d"%cnt, len(rd.seq),
                        rd.description, rd.seq))
                    maplist_simple.append("%s\t%d\t%s"%("seq_%d"%cnt, len(rd.seq),
                        rd.description))
                    if not g_params['isForceRun']:
                        md5_key = hashlib.md5(rd.seq).hexdigest()
                        subfoldername = md5_key[:2]
                        md5_link = "%s/%s/%s"%(path_md5cache, subfoldername, md5_key)
                        if os.path.exists(md5_link):
                            # create a symlink to the cache
                            rela_path = os.path.relpath(md5_link, outpath_result) #relative path
                            os.chdir(outpath_result)
                            os.symlink(rela_path, subfoldername_this_seq)

                            if os.path.exists(outpath_this_seq):
                                runtime = 0.0 #in seconds
                                predfile = "%s/query.predzinc.predict"%( outpath_this_seq)
                                (numZB, cntHomo) = CountNumPredZB(predfile, threshold=ZB_SCORE_THRESHOLD)
                                # info_finish has 7 items
                                info_finish = [ "seq_%d"%cnt,
                                        str(len(rd.seq)), str(numZB),
                                        str(cntHomo), "cached", str(runtime),
                                        rd.description]
                                myfunc.WriteFile("\t".join(info_finish)+"\n",
                                        finished_seq_file, "a", isFlush=True)
                                isSkip = True

                    if not isSkip:
                        # first try to delete the outfolder if exists
                        if os.path.exists(outpath_this_seq):
                            try:
                                shutil.rmtree(outpath_this_seq)
                            except OSError:
                                pass
                        origIndex = cnt
                        numTM = 0
                        toRunDict[origIndex] = [rd.seq, numTM, rd.description] #init value for numTM is 0

                    cnt += 1
                recordList = hdl.readseq()
            hdl.close()
        myfunc.WriteFile("\n".join(maplist_simple)+"\n", mapfile)

        sortedlist = sorted(toRunDict.items(), key=lambda x:x[1][1], reverse=True)
        #format of sortedlist [(origIndex: [seq, numTM, description]), ...]

        # submit sequences one by one to the workflow according to orders in
        # sortedlist

        for item in sortedlist:
#             g_params['runjob_log'].append("tmpdir = %s"%(tmpdir))
            #cmd = [script_getseqlen, infile, "-o", tmp_outfile , "-printid"]
            origIndex = item[0]
            seq = item[1][0]
            description = item[1][2]

            outpath_this_seq = "%s/%s"%(outpath_result, "seq_%d"%origIndex)
            tmp_outpath_this_seq = "%s/%s"%(tmp_outpath_result, "seq_%d"%(0))
            if os.path.exists(tmp_outpath_this_seq):
                try:
                    shutil.rmtree(tmp_outpath_this_seq)
                except OSError:
                    pass
            try:
                os.makedirs(tmp_outpath_this_seq)
            except OSError:
                g_params['runjob_err'].append("Failed to create the tmp_outpath_this_seq %s"%(tmp_outpath_this_seq))
                continue

            seqfile_this_seq = "%s/%s"%(tmp_outpath_result, "query_%d.fa"%(origIndex))
            seqcontent = ">%d\n%s\n"%(origIndex, seq)
            myfunc.WriteFile(seqcontent, seqfile_this_seq, "w")

            if not os.path.exists(seqfile_this_seq):
                g_params['runjob_err'].append("failed to generate seq index %d"%(origIndex))
                continue

            if not os.path.exists("%s/seq.fa"%(tmp_outpath_this_seq)):
                try:
                    shutil.copyfile(seqfile_this_seq, "%s/seq.fa"%(tmp_outpath_this_seq))
                except OSError:
                    pass

            numCPU = 4
            cmd = [runscript, seqfile_this_seq,  "-outpath",
                    tmp_outpath_this_seq, "-cpu", str(numCPU), "-db",
                    predzinc_db, "-blastdb", blastdb]
            g_params['runjob_log'].append(" ".join(cmd))
            begin_time = time.time()
            try:
                rmsg = subprocess.check_output(cmd)
            except subprocess.CalledProcessError, e:
                g_params['runjob_err'].append(str(e)+"\n")
                g_params['runjob_err'].append(rmsg + "\n")
                pass
            end_time = time.time()
            runtime_in_sec = end_time - begin_time

            if os.path.exists(tmp_outpath_this_seq):
                cmd = ["mv","-f", tmp_outpath_this_seq, outpath_this_seq]
                isCmdSuccess = False
                try:
                    subprocess.check_output(cmd)
                    isCmdSuccess = True
                except subprocess.CalledProcessError, e:
                    msg =  "Failed to run prediction for sequence No. %d\n"%(origIndex)
                    g_params['runjob_err'].append(msg)
                    g_params['runjob_err'].append(str(e)+"\n")
                    pass

                if isCmdSuccess:
                    runtime = runtime_in_sec #in seconds
                    predfile = "%s/query.predzinc.predict"%( outpath_this_seq)
                    (numZB, cntHomo) = CountNumPredZB(predfile, threshold=ZB_SCORE_THRESHOLD)
                    # info_finish has 7 items
                    info_finish = [ "seq_%d"%origIndex, str(len(seq)), str(numZB),
                            str(cntHomo), "newrun", str(runtime), description]
                    myfunc.WriteFile("\t".join(info_finish)+"\n",
                            finished_seq_file, "a", isFlush=True)

                    # create or update the md5 cache
                    # create cache only on the front-end
                    if g_params['base_www_url'].find("predzinc") != -1:
                        md5_key = hashlib.md5(seq).hexdigest()
                        subfoldername = md5_key[:2]
                        md5_subfolder = "%s/%s"%(path_md5cache, subfoldername)
                        md5_link = "%s/%s/%s"%(path_md5cache, subfoldername, md5_key)
                        if os.path.exists(md5_link):
                            try:
                                os.unlink(md5_link)
                            except:
                                pass
                        subfolder_md5 = "%s/%s"%(path_md5cache, subfoldername)
                        if not os.path.exists(subfolder_md5):
                            try:
                                os.makedirs(subfolder_md5)
                            except:
                                pass

                        rela_path = os.path.relpath(outpath_this_seq, md5_subfolder) #relative path
                        try:
                            os.chdir(md5_subfolder)
                            os.symlink(rela_path,  md5_key)
                        except:
                            pass


        all_end_time = time.time()
        all_runtime_in_sec = all_end_time - all_begin_time

        if len(g_params['runjob_log']) > 0 :
            rt_msg = myfunc.WriteFile("\n".join(g_params['runjob_log'])+"\n", runjob_logfile, "a")
            if rt_msg:
                g_params['runjob_err'].append(rt_msg)

        datetime = time.strftime("%Y-%m-%d %H:%M:%S")
        if os.path.exists(finished_seq_file):
            rt_msg = myfunc.WriteFile(datetime, finishtagfile)
            if rt_msg:
                g_params['runjob_err'].append(rt_msg)

# now write the text output to a single file
        dumped_resultfile = "%s/%s"%(outpath_result, "query.predzinc.txt")
        statfile = "%s/%s"%(outpath_result, "stat.txt")
        WriteTextResultFile(dumped_resultfile, outpath_result, maplist,
                all_runtime_in_sec, statfile=statfile)


        # now making zip instead (for windows users)
        pwd = os.getcwd()
        os.chdir(outpath)
#             cmd = ["tar", "-czf", tarball, resultpathname]
        cmd = ["zip", "-rq", zipfile, resultpathname]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            g_params['runjob_err'].append(str(e))
            pass
        os.chdir(pwd)


    isSuccess = False
    if (os.path.exists(finishtagfile) and os.path.exists(zipfile_fullpath)):
        isSuccess = True
        # delete the tmpdir if succeeded
        shutil.rmtree(tmpdir) #DEBUG, keep tmpdir
    else:
        isSuccess = False
        failtagfile = "%s/runjob.failed"%(outpath)
        datetime = time.strftime("%Y-%m-%d %H:%M:%S")
        rt_msg = myfunc.WriteFile(datetime, failtagfile)
        if rt_msg:
            g_params['runjob_err'].append(rt_msg)

# send the result to email
# do not sendmail at the cloud VM
    if (g_params['base_www_url'].find("predzinc") != -1 and
            myfunc.IsValidEmailAddress(email)):
        from_email = "predzinc@bioshu.se"
        to_email = email
        subject = "Your result for PREDZINC JOBID=%s"%(jobid)
        if isSuccess:
            bodytext = """
Your result is ready at %s/pred/result/%s

Thanks for using PREDZINC

        """%(g_params['base_www_url'], jobid)
        else:
            bodytext="""
We are sorry that your job with jobid %s is failed.

Please contact %s if you have any questions.

Attached below is the error message:
%s
            """%(jobid, contact_email, "\n".join(g_params['runjob_err']))
        g_params['runjob_log'].append("Sendmail %s -> %s, %s"% (from_email, to_email, subject)) #debug
        rtValue = myfunc.Sendmail(from_email, to_email, subject, bodytext)
        if rtValue != 0:
            g_params['runjob_err'].append("Sendmail to {} failed with status {}".format(to_email, rtValue))

    if len(g_params['runjob_err']) > 0:
        rt_msg = myfunc.WriteFile("\n".join(g_params['runjob_err'])+"\n", runjob_errfile, "w")
        return 1
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    infile = ""
    tmpdir = ""
    email = ""
    jobid = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-tmpdir", "--tmpdir"] :
                (tmpdir, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-jobid", "--jobid"] :
                (jobid, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-baseurl", "--baseurl"] :
                (g_params['base_www_url'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-email", "--email"] :
                (email, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-force", "--force"]:
                g_params['isForceRun'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if jobid == "":
        print >> sys.stderr, "%s: jobid not set. exit"%(sys.argv[0])
        return 1

    if myfunc.checkfile(infile, "infile") != 0:
        return 1
    if outpath == "":
        print >> sys.stderr, "outpath not set. exit"
        return 1
    elif not os.path.exists(outpath):
        try:
            subprocess.check_output(["mkdir", "-p", outpath])
        except subprocess.CalledProcessError, e:
            print >> sys.stderr, e
            return 1
    if tmpdir == "":
        print >> sys.stderr, "tmpdir not set. exit"
        return 1
    elif not os.path.exists(tmpdir):
        try:
            subprocess.check_output(["mkdir", "-p", tmpdir])
        except subprocess.CalledProcessError, e:
            print >> sys.stderr, e
            return 1

    return RunJob(infile, outpath, tmpdir, email, jobid, g_params)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['runjob_log'] = []
    g_params['runjob_err'] = []
    g_params['isForceRun'] = False
    g_params['base_www_url'] = ""
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
