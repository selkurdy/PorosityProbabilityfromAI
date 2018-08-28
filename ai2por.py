"""
python ai2por.py ffps.csv --listcols
python ai2por.py ffps.csv --datacols 14 11
python ai2por.py ffps.csv --datacols 14 11 --segyfname ai.segy

python ai2por.py ffps.csv --datacols 14 11 --hideplot
python ai2por.py ffps.csv --datacols 14 11 --segyfname aix.sgy --hideplot --percentile 20 55 80
python ai2por.py ffps.csv --datacols 14 11 --segyfname aix.sgy --hideplot --patv --porval 0.18

python ai2por.py ffps.csv --datacols 14 11 --segyfname XAGP17_086_AI_Absolute_SK307_Baram.sgy --bins 20
"""
import  os.path
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import segyio
from shutil import copyfile




def transformp(sample,binedges,databins,por,ndata=10):
    """
    Do the transformation
    sample is one data point
    por is porosity as fraction, i.e. 0.2 value to compute its probability  
    exceedence is assumed, i.e. 1-CDF
    return value is probability of por or higher
    """
    
    j = np.digitize(sample,binedges)
    ai = stats.percentileofscore(databins[j],por)
    # returned ai is percentage from 0 to 100. Need to divide by 100
    return 1.- ai/100.0

        
def transform(sample,binedges,databins,pcnt,ndata=10):
    """
    Do the transformation
    sample is one data point
    pct is an array of percentiles needed, e.g. 10,50,90
    return value is the same number of pcnt
    """
    
    j = np.digitize(sample,binedges)
    # print('bin#',j)
    
    ai = np.percentile(databins[j],pcnt)
    # print(ai)
    return ai

        

def mystat10(a):
    return np.percentile(a,10)
def mystat50(a):
    return np.percentile(a,50)
def mystat90(a):
    return np.percentile(a,90)



def getcommandline():
    parser = argparse.ArgumentParser(description='Compute various percentiles of porosity from acoustic impedence')
    parser.add_argument('wellcsv', help='csv of well AI & well porosity')
    parser.add_argument('--listcols',action='store_true',default=False,
            help='list columns in csv file')
    parser.add_argument('--datacols', nargs= 2,type=int,default=[0,1],help='AI Por column # dfv = 0 1 ')
    parser.add_argument('--bins',type=int,default=10,help='# of AI bins. default = 10')
    parser.add_argument('--segyfname',help='segy file of AI to convert to porosity')
    parser.add_argument('--minbindata',type=int,default= 10, help='Minimum data per bin, if less use the previous bin. default=10')
    # parser.add_argument('--percentile',type=float,nargs= 3,default=[10,50,90],help='percentile. default=10 50 90')
    parser.add_argument('--percentile',type=float,default=50.0,help=' one percentile. default=50')
    parser.add_argument('--patv',action='store_true',default=False,
            help='percentile or probaility at porosity value. default= false, i.e. output for percentile')
    parser.add_argument('--porval',type =float,default=0.2,help='percentile or probaility at porosity value. Exceedence is assumed, i.e. 1-CDF. default = .2 or probaility of porosity greater than 20 percent')
    parser.add_argument('--outdir',help='output directory,default= same dir as input')
    parser.add_argument('--hideplot',action='store_true',default=False,
                        help='Only save to pdf. default =show and save')
    
    
    result=parser.parse_args()
    if not result.wellcsv:
        parser.print_help()
        exit()
    else:
        return result

        
        
def main():
    cmdl= getcommandline()
    w = pd.read_csv(cmdl.wellcsv)
    if cmdl.listcols:
        print(list(enumerate(w.columns)))
    else:
        # print(w.head())
        aiw = w.iloc[:,cmdl.datacols[0]].values
        porw = w.iloc[:,cmdl.datacols[1]].values
        # fig,ax = plt.subplots()
        # plt0 = ax.scatter(aiw,porw,alpha=.3,s=10)
        # ax.set_xlabel('Well Acoustic Impedence')
        # ax.set_ylabel('Well Porosity')
        # plt.show()
        dirsplit,fextsplit= os.path.split(cmdl.wellcsv)
        fname,fextn= os.path.splitext(fextsplit)
        
        bin_p10, bin_edges, binnumber = stats.binned_statistic(aiw,porw,statistic=mystat10,bins =cmdl.bins)
        bin_p50, bin_edges, binnumber = stats.binned_statistic(aiw,porw,statistic=mystat50,bins =cmdl.bins)
        bin_p90, bin_edges, binnumber = stats.binned_statistic(aiw,porw,statistic=mystat90,bins =cmdl.bins)
        be = np.diff(bin_edges) * .5 +bin_edges[:-1]
        fig,ax = plt.subplots()
        plt0 = ax.scatter(aiw,porw,alpha=.3,s=10)
        plt1 = ax.scatter(be,bin_p10,alpha=.7,s=40,c='r',label='P10')
        plt1 = ax.scatter(be,bin_p50,alpha=.7,s=40,c='k',label='P50')
        plt1 = ax.scatter(be,bin_p90,alpha=.7,s=40,c='m',label='P90')
        ax.set_xlabel('Well Acoustic Impedence')
        ax.set_ylabel('Well Porosity')
        ax.legend()
 
        if not cmdl.hideplot:
            plt.show()
        if cmdl.outdir:
            pdfcl = os.path.join(cmdl.outdir,fname) +"_aipor.pdf" 
        else:
            pdfcl = os.path.join(dirsplit,fname) +"_aipor.pdf" 
        fig = ax.get_figure()
        fig.savefig(pdfcl)                

        porbinned=list()
        aibinned = list()
        numporbin = list()
        numaibin = list()
        for i in range(bin_p10.shape[0]):
            data = porw[binnumber == i]
            porbinned.append(data)
            data = aiw[binnumber == i]
            aibinned.append(data)
            numporbin.append(porbinned[i].size)
            numaibin.append(aibinned[i].size)
            # print (i,porbinned[i].size,aibinned[i].size)
    
        por2 = list()
        numpor2 = list()
        
        ai2 = list()
        numai2 =list()
        
        aibe =list()
        porbe =list()
        
        for i in range(len(porbinned)):
            if porbinned[i].size > 10:
                por2.append(porbinned[i])
                ai2.append(aibinned[i])
        # if por2[0].min() > porbinned[0].min():
            # porbe.append(porbinned[0].min())
            # aibe.append(aibinned[0].min())
        # else:
            # porbe.append(por2[0].min())        
            # aibe.append(ai2[0].min())        
            
            
        for i in range(1,len(por2)-1):
            porbe.append(por2[i].max())
            aibe.append(ai2[i].max())
            
            
        # if por2[-1].max() < porbinned[-1].max():
            # porbe.append(porbinned[-1].min())
            # aibe.append(aibinned[-1].min())
            
        for i in range(len(por2)):
            numpor2.append(len(por2[i]))
            numai2.append(len(ai2[i]))
            # binnum.append(i)
        print(len(por2),len(ai2))
        print(len(porbe),len(aibe))


        
        databins = por2 #por bins
        cols = ['binnum','AIstart','AIend','PORstart','PORend','PORnum']
        aistart = []
        aiend =[]
        porstart =[]
        porend =[]
        pornum =[]
        binnum = list()
        for i in range(1,len(porbe)):
            # databins.append(por2[i])
            binnum.append(i)
            aistart.append(aibe[i-1])
            aiend.append(aibe[i])
            porstart.append(porbe[i-1])
            porend.append(porbe[i])
            pornum.append(por2[i].size)
        binsdf = pd.DataFrame({cols[0]:binnum,cols[1]:aistart,cols[2]:aiend,cols[3]:porstart,cols[4]:porend,cols[5]:pornum})
        binsdf1 = binsdf[cols].copy()
        print(binsdf1.head(cmdl.bins))
        if cmdl.outdir:
            binsdfname = os.path.join(cmdl.outdir,fname) +"_aipor.csv" 
        else:
            binsdfname = os.path.join(dirsplit,fname) +"_aipor.csv" 
        binsdf1.to_csv(binsdfname,index=False)
        print('Successfully wrote %s'%binsdfname)
        
        
        if cmdl.segyfname:
            # in_fname = cmdl.segyfname
            dirsplit,fextsplit= os.path.split(cmdl.segyfname)
            fname,fextn= os.path.splitext(fextsplit)
            if cmdl.patv:
                if cmdl.outdir:
                    outfnamep = os.path.join(cmdl.outdir,fname) +"_pprob%.0f.sgy"%(cmdl.porval *100 )
                else:
                    outfnamep = os.path.join(dirsplit,fname) +"_pprob%.0f.sgy"% (cmdl.porval * 100)
                print('Copying file, please wait ........')
                start_copy = datetime.now()
                copyfile(cmdl.segyfname, outfnamep) 
                end_copy = datetime.now()
                print('Duration of copying: {}'.format(end_copy - start_copy))
                
                start_process = datetime.now()
                with segyio.open( outfnamep, "r+" ) as srcp:
                    # Memory map file for faster reading (especially if file is big...)                    
                    srcp.mmap()
                    trnum = 0
                    for tracep in srcp.trace:
                        if trnum % 100 == 0:
                            print('Trace #: {}'.format(trnum))
                        for i in range(tracep.size):
                            # ai = transformp(tracep[i],bin_edges,databins,cmdl.porval,ndata=cmdl.minbindata)
                            ai = transformp(tracep[i],aibe,databins,cmdl.porval,ndata=cmdl.minbindata)
                            tracep[i] = ai
                            # print('%10.5f'%trace[i])
                        srcp.trace[trnum] = tracep
                        # print('src.trace',src.trace[trnum])
                        srcp.flush()
                        trnum += 1
                print('Successfully wrote %s ' % outfnamep)
                end_process = datetime.now()
                print('Duration of processing: {}'.format(end_process - start_process))
                
            
            else:
                if cmdl.outdir:
                    outfnamep = os.path.join(cmdl.outdir,fname) +"_porp%.0f.sgy"%cmdl.percentile
                else:
                    outfnamep = os.path.join(dirsplit,fname) +"_porp%.0f.sgy"%cmdl.percentile

                print('Copying files, please wait ........')
                start_copy = datetime.now()
                copyfile(cmdl.segyfname, outfnamep)            
                end_copy = datetime.now()
                print('Duration of copying: {}'.format(end_copy - start_copy))
                

                start_process = datetime.now()
                with segyio.open( outfnamep, "r+" ) as srcp:
                    # Memory map file for faster reading (especially if file is big...)                    
                    srcp.mmap()
                    trnum = 0
                    for tracep in srcp.trace:
                        if trnum % 100 == 0:
                            print('Trace #: {}'.format(trnum))
                        for i in range(tracep.size):
                            ai = transform(tracep[i],aibe,databins,cmdl.percentile,ndata=cmdl.minbindata)
                            # print('ai after:',ai)


                            tracep[i] = ai
                            # print('%10.5f'%trace[i])
                            
                        srcp.trace[trnum] = tracep  
                        # print(trnum,srcp.trace[trnum])
                        trnum += 1
                print('Successfully wrote %s ' % outfnamep)
                end_process = datetime.now()
                print('Duration of processing: {}'.format(end_process - start_process))

if __name__== '__main__':
    main() 