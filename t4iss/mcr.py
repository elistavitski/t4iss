#-*- coding: utf-8 -*-



import numpy as np 
import shutil,subprocess,os,time

from pylab import *
from matplotlib import gridspec
from matplotlib import pyplot as plt

import scipy.io as sio




def run_EFA(datafile,ncomp,defaults,plot=True):
    
    if defaults['octave_path']:
        sp_call_str = ' export OCTAVE_PATH='+defaults['mcr_path'] +'/octave_version/:$OCTAVE_PATH; octave -qW --eval '+' "efa_t4iss(\''+datafile+'\', '+str(ncomp)+')" > efa.out  '
    elif defaults['matlab_path']:
        sp_call_str = 'export MATLABPATH=' + defaults['mcr_path']+'/matlab_version/:$MATLABPATH; matlab -nodesktop  -nosplash -r '+' "efa_t4iss(\''+datafile+'\', '+str(ncomp)+')" > efa.out  '
    else:
        print('MATLAB or OCTAVE path is not defined. Exitting....')
        return
    
    # check datafile
    print(sp_call_str)
    msg = 'Error: \n datafile does not exist. Please check it !! '
    try:
        if not os.path.isfile(datafile):
            print(msg)
            return
    except Exception as exc:
        print(msg)
        print(exc)
        return 


    if os.path.isfile('efa.mat'): os.remove('efa.mat')
    subprocess.call(sp_call_str, shell=True)
    time.sleep(3)

    
    if plot:
        # EFA plots
        mat_contents = sio.loadmat('efa.mat')
        xbackward  = mat_contents['xbackward']
        ebl        = mat_contents['ebl']
        xforward   = mat_contents['xforward']  
        efl        = mat_contents['efl']  
        e  = mat_contents['e']
        
        
        fig = plt.figure(figsize=(9,6))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,2], height_ratios=[1] )
        gs.update(top=0.92, bottom=0.1, left=0.07, right=0.95, wspace=0.3, hspace=0.2)
        
        
        ax=fig.add_subplot(gs[0]) 
        for i in range(ebl.shape[1]):
            ax.plot(xbackward[0],ebl[:,i],'-r')    
        for i in range(efl.shape[1]):
            ax.plot(xforward[0],efl[:,i],'-k')     
        ax.set_title('EFA'); 
        ax.set_ylabel('log(eigenvalues)')
        
        
        import warnings; warnings.simplefilter('ignore')
        ax=fig.add_subplot(gs[1]) 
        for i in range(ebl.shape[1]):
            ax.plot(xbackward[0],np.sqrt(ebl[:,i]),'-r')    
        for i in range(efl.shape[1]):
            ax.plot(xforward[0],np.sqrt(efl[:,i]),'-k') 
        ax.set_ylabel('Singular values');
        
        
        ax=fig.add_subplot(gs[2]) 
        for i in range(e.shape[1]):
            ax.plot(e[:,i],'-o',label='comp-'+str(i+1)) 
        ax.set_title('Initial concentrations');
        ax.legend(loc='best', fontsize=7);    
        
    return
    






def do_ALS(nit=300,thr=0.1,normalization=0,save_label=None,print_log=True):
    
    if os.path.exists('als.mat'): os.remove('als.mat')

    if t4iss_defaults['octave_path']:
        sp_call_str = ' export OCTAVE_PATH='+mcr_path+'/octave_version/:$OCTAVE_PATH; octave -qW --eval '+' "als_t4iss('+str(nit)+','+str(thr)+','+str(normalization)+')" > als.out  '
    elif t4iss_defaults['matlab_path']:
        sp_call_str = ' export MATLABPATH='+mcr_path+'/matlab_version/:$MATLABPATH; matlab -nodesktop  -nosplash -r '+' "als_t4iss('+str(nit)+','+str(thr)+','+str(normalization)+')" > als.out  '
    else:
        print('MATLAB or OCTAVE path is not defined. Exitting....')
        return

    
    # check datafile
    msg = 'Error: \n efa.mat does not exist. Did you run EFA ? '
    try:
        if not os.path.isfile('efa.mat'):
            print(msg)
            return
    except Exception as exc:
        print(msg)
        print(exc)
        return     
      
   
    print('ALS run started.... \n') 
    time_start = time.time()
    subprocess.call(sp_call_str, shell=True)
    time.sleep(5)
    als_time = time.time() - time_start
    print('ALS run completed in {:.3f} seconds \n'.format(als_time))
    if print_log: print(open('als.out').read())
    
    
    if save_label is None:
        print('ALS optimization log can be found in als.out\n') 
    else:
        shutil.move('als.mat',save_label+'.mat')
        shutil.move('als.out',save_label+'.out')
        print('ALS optimization log can be found in '+save_label+'.out')
        
    return







def plot_results(read_from=None,export_data=True):
    
    if read_from is None:
        mat_contents = sio.loadmat('als.mat')
    else:
        mat_contents = sio.loadmat(read_from+'.mat')
        
    
    copt  = mat_contents['copt']
    sopt  = mat_contents['sopt']
    E  = mat_contents['Energy']

    
    fig = plt.figure(figsize=(9,5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.5,1], height_ratios=[1] )
    gs.update(top=0.95, bottom=0.1, left=0.08, right=0.95, wspace=0.1, hspace=0.05)
    
    ax = fig.add_subplot(gs[0]); ax.grid(True)
    markers = ['-ro','--gs','-.bd','-m^','-c<','--k>']
    for i in range(len(copt[0])):
        ax.plot(copt[:,i],markers[i], ms=5, alpha=0.8, label='component-'+str(i+1))

    ax.set_xlabel('Spectra #')
    ax.set_ylabel('Concentration')
    ax.legend(loc='best',fontsize=10,ncol=1)
    
    
    ax = fig.add_subplot(gs[1]); ax.grid(False)
    markers = ['-r','-g','-b','-m','-c','--k']
    for i in range(len(sopt)):
        ax.plot(E,sopt[i,:], markers[i], lw=3, alpha=0.6, label='component-'+str(i+1))
    ax.set_xlabel('Energy [eV]')
    ax.set_yticks([])
    ax.legend(loc='best',fontsize=10,ncol=1)
    
    if export_data:
        np.savetxt('concentrations.dat', copt, fmt='%1.4e')  
        sopt = np.column_stack( (E,sopt.T) )
        np.savetxt('components.dat', sopt, fmt='%1.4e') 
        print('concentrations and pure spectra in plain text were exported as concentrations.dat and components.dat \n')
    




    
