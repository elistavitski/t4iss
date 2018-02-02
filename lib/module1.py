# -*- coding: utf-8 -*-
  
import numpy as np
import os,sys,shutil,subprocess,pickle,json
from os.path import join

from pymatgen.core.periodic_table import Element
from pymatgen.analysis.xas.spectrum import XANES
from pymatgen.symmetry.analyzer import *
import pymatgen as mg
from pymatgen.io.vasp.sets import MPRelaxSet

from numpy import linalg as LA
from scipy.interpolate import InterpolatedUnivariateSpline

from pylab import *
from matplotlib import gridspec
from matplotlib import pyplot as plt

from .feff_tools import write_feffinp
from .PT import get_c



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def get_XANES(mpr,mpid,absorbing_atom,run_feff=None,dbroot=None,plot=None,n_cpu=None,export_figure=None):
    
    here = os.getcwd()
    
    if dbroot is None:
        dbroot = join(here+'/data/XANES')
        if not os.path.isdir(dbroot):
            os.makedirs(dbroot,exist_ok=True)

    if run_feff is None:
        run_feff = False
    elif run_feff:
        if os.name == 'nt':
            print('Windows OS is not supported for run_feff=True \nSorry......' )
            run_feff = False
        else:
            rFMS=8.1
            rSCF=6.1
            corehole='RPA'
            feff_cmd = here+'/lib/feff_cmd.sh'
        
    if plot is None:
        plot = True
        
    if n_cpu is None:
        n_cpu = 2   
        
    if export_figure is None:
        export_figure = False            
        
    os.chdir(dbroot)  
    
    if not os.path.isdir(mpid):
        data = mpr._make_request('/xas/spectra_for',{"material_ids": json.dumps([mpid]),
                "elements": json.dumps([absorbing_atom]) })
        if not data:
            mp_xas_avail=False
            if run_feff:
                os.makedirs(mpid,exist_ok=True)
                os.chdir(mpid)
            else:
                print('XANES is not available in MP or local database.')
                try:
                    structure = mpr.get_structure_by_material_id(mpid,final=True)
                    os.makedirs(mpid,exist_ok=True)
                    os.chdir(mpid)
                except:
                    print('structure for '+mpid+' is not available in MP.\nExitting....')
                    os.chdir(here) 
                    return      
        else:
            mp_xas_avail=True
            data = mpr.get_data(mpid, data_type="feff", prop="xas")
            os.makedirs(mpid,exist_ok=True)
            os.chdir(mpid)
    else:
        os.chdir(mpid)
        
    if not os.path.isfile('CONTCAR'):
        structure = mpr.get_structure_by_material_id(mpid,final=True)
        structure.to(fmt='poscar',filename='CONTCAR')        
    structure = mg.Structure.from_file("CONTCAR")
    
    finder = SpacegroupAnalyzer(structure)
    symmetrized_structure = finder.get_symmetrized_structure()
    [sites, indices]  = symmetrized_structure.equivalent_sites, symmetrized_structure.equivalent_indices
    
    

    feff_todos = []
    for i,s in enumerate(sites):
        if s[0].species_string is absorbing_atom:
            f = 'feff_{:03d}_{}'.format(indices[i][0]+1,absorbing_atom)
            if os.path.isdir(f):
                os.chdir(f)
                if not os.path.isfile('xanes.pkl'):                    
                    data = mpr._make_request('/xas/spectra_for',{"material_ids": json.dumps([mpid]),
                            "elements": json.dumps([absorbing_atom]) })
                    if not data:
                        mp_xas_avail=False
                    else: mp_xas_avail=True
                    if mp_xas_avail:
                        data = mpr.get_data(mpid, data_type="feff", prop="xas")
                        for xas_doc in data[0]['xas']:
                            abs_atom = xas_doc['absorbing_atom']
                            if abs_atom is indices[i][0]:                        
                                x, y = xas_doc['spectrum']
                                struct = xas_doc["structure"]
                                edge = xas_doc["edge"]
                                xanes = XANES(x, y, struct, Element(absorbing_atom), edge='K')
                                pickle.dump(xanes, open('xanes.pkl', 'wb'))
                                out = np.column_stack( (x,y) )
                                #np.savetxt('xanes.dat', out, fmt="%10.3f %6.4e") 
                    elif os.path.isfile('xmu.dat'):
                        x, y = np.loadtxt('xmu.dat', unpack=True, comments='#', usecols=(0,3), skiprows=0)
                        xanes = XANES(x, y, structure, Element(absorbing_atom), 'K')
                        pickle.dump(xanes, open('xanes.pkl', 'wb'))
                        out = np.column_stack( (x,y) )
                        np.savetxt('xanes.dat', out, fmt="%10.3f %6.4e")
                    elif run_feff:
                        write_feffinp('../CONTCAR',cai=indices[i][0],dmax=rFMS+2,rFMS=rFMS,rSCF=rSCF,corehole=corehole)
                        feff_todos.append(os.getcwd())
                else:
                    os.chdir('..')
            else:
                data = mpr._make_request('/xas/spectra_for',{"material_ids": json.dumps([mpid]),
                                        "elements": json.dumps([absorbing_atom]) })
                if data:
                    os.makedirs(f,exist_ok=True)
                    os.chdir(f)
                    data = mpr._make_request('/xas/spectra_for',{"material_ids": json.dumps([mpid]),
                            "elements": json.dumps([absorbing_atom]) })
                    if not data:
                        mp_xas_avail=False
                    else: mp_xas_avail=True
                    if mp_xas_avail:
                        data = mpr.get_data(mpid, data_type="feff", prop="xas")
                        for xas_doc in data[0]['xas']:
                            abs_atom = xas_doc['absorbing_atom']
                            if abs_atom is indices[i][0]:                        
                                x, y = xas_doc['spectrum']
                                struct = xas_doc["structure"]
                                edge = xas_doc["edge"]
                                xanes = XANES(x, y, struct, Element(absorbing_atom), edge='K')
                                pickle.dump(xanes, open('xanes.pkl', 'wb'))
                                out = np.column_stack( (x,y) )
                                #np.savetxt('xanes.dat', out, fmt="%10.3f %6.4e") 
                    os.chdir('..')
                elif run_feff:
                    os.makedirs(f,exist_ok=True)
                    os.chdir(f)
                    write_feffinp('../CONTCAR',cai=indices[i][0],dmax=rFMS+2,rFMS=rFMS,rSCF=rSCF,corehole=corehole)
                    feff_todos.append(os.getcwd())
                    os.chdir('..')
                else:
                    break

    if run_feff:
        chunks = [feff_todos [i:i+n_cpu] for i in range(0, len(feff_todos), n_cpu) ]    
        for c in chunks:
            f=open('feff.sh',"w+")
            for i in c:
                f.write('cd '+i+'\n')
                f.write(feff_cmd+' &\n')
                print('running feff at '+i)        
            f.write('wait')
            f.close()
            subprocess.call(' chmod +x feff.sh ', shell=True)
            subprocess.call(' ./feff.sh ', shell=True)
            if os.path.exists('feff.sh'): os.remove('feff.sh')

    spectra = []
    for i,s in enumerate(sites):
        if s[0].species_string is absorbing_atom:
            f = 'feff_{:03d}_{}'.format(indices[i][0]+1,absorbing_atom)   
            if os.path.isdir(f):
                os.chdir(f)
                if os.path.isfile('xanes.pkl'):
                    xanes = pickle.load(open('xanes.pkl', 'rb'))
                elif os.path.isfile('xmu.dat'):
                    abs_specie = absorbing_atom
                    x, y = np.loadtxt('xmu.dat', unpack=True, comments='#', usecols=(0,3), skiprows=0)
                    xanes = XANES(x, y, structure, Element(abs_specie), 'K')
                    pickle.dump(xanes, open('xanes.pkl', 'wb'))
                    out = np.column_stack( (x,y) )
                    np.savetxt('xanes.dat', out, fmt="%10.3f %6.4f")                    
                else: 
                    xanes = []                
                s_weight = len(indices[i])  
                local_environment = symmetrized_structure.get_sites_in_sphere(symmetrized_structure[indices[i][0]].coords,10.1)
                spectra.append([s_weight,xanes,local_environment])
                os.chdir('..')                
            else:
                xanes = []
                s_weight = len(indices[i])  
                local_environment = symmetrized_structure.get_sites_in_sphere(symmetrized_structure[indices[i][0]].coords,10.1)
                spectra.append([s_weight,xanes,local_environment])                            
    if plot:
        try:
            os.chdir(here)
            plot_XANES(spectra,mpid,export_figure)
        except:
            print('Error: \n Unable to plot. Something is wrong...')
            os.chdir(here)
            return
                
    os.chdir(here)
    return


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_XANES(sp,mpid=None,export_figure=None):
    
   
    fig = plt.figure(figsize=(9,5+len(sp)/4))    
    gs1 = gridspec.GridSpec(1, 2, width_ratios=[2,2] )
    gs1.update(top=0.90, bottom=0.1, left=0.07, right=0.97, wspace=0.15, hspace=0.05)
    gs2 = gridspec.GridSpec(1, 2, width_ratios=[2,2] )
    gs2.update(top=0.95, bottom=0.1, left=0.02, right=0.97, wspace=0.15, hspace=0.05)
    
    ax=fig.add_subplot(gs1[0])    
    left, bottom, width, height = [0.33, 0.15, 0.15, 0.3]
    inset_ax = fig.add_axes([left, bottom, width, height])        
    # for labels
    env0 = sp[0][2]
    def getKey1(item): return item[1]
    env0 = sorted(env0, key=getKey1)
    # get species
    ss = []
    for i in env0:
        ss.append(i[0].specie.name)    
    labels = list(set(ss))    
    for i in labels:
        c = get_c(i)
        ax.plot(0,0,'o',color=c,label=i,ms=8)  
    ax.legend(loc='upper left',fontsize=12,ncol=1)    
    if mpid:
        ax.set_title(mpid)
    else:
        mpid = 'mp'
        
    s = 0
    for i in sp:
        multip = str(i[0])
        site_text = env0[0][0].specie.name+'-'+str(s+1)+'\n(x'+multip+')'
        ax.annotate(site_text,(-2,s), fontsize=8)        
        env = i[2]
        env = sorted(env, key=getKey1)    
        ss = []
        ds = []
        for i in env:
            ss.append(i[0].specie.name)  
            ds.append(i[1])
        ds = np.array(ds,np.float)
            
        ax.plot(s+ds[0:21],'k-')
        for i,d in enumerate(ds[0:21]):
            c = get_c(ss[i])
            ax.plot( i, s+d, 'o', color=c, ms=9, alpha=0.8 ) 
            inset_ax.plot( i, d, 'o', color=c, ms=6, alpha=0.8 )
        s += 1         
    
    ax.set_xticks(list(range(1,21)))
    ax.set_xlim([-2.5,21])
    ax.set_xlabel('Neighbour index #')
    ax.set_ylabel('Distance to absorbing atom ($\AA$)')
    
    if sp[0][1]:     
        ax=fig.add_subplot(gs2[1])
        minmaxs = []
        yshift=0
        for s,i in enumerate(sp):
            multip = str(i[0])
            xas_text = env0[0][0].specie.name+'-'+str(s+1)+'\n(x'+multip+')'
            ax.plot(i[1].energy,yshift+i[1].intensity,'-')
            ax.annotate(xas_text,(i[1].energy[-1],yshift+i[1].intensity[-1]), fontsize=8)
            yshift = yshift + 0.5
            minmaxs.append([i[1].energy[0],i[1].energy[-1]])
        minmaxs = np.array(minmaxs)  
          
        e_int = np.linspace(max(minmaxs[:,0]),min(minmaxs[:,1]),300)
        multips = 0
        ts = e_int*0
        for i in sp:
            multips = multips + i[0]
            f = InterpolatedUnivariateSpline(i[1].energy,i[1].intensity)
            i_int = f(e_int)
            ts = ts + i_int
        ts = ts/ts[-1]
        if len(sp) > 1:
            ax.plot(e_int,yshift+ts,'k-')
            ax.annotate('total\n',(e_int[-1],yshift+ts[-1]), fontsize=8)
        
        ax.set_yticks([])    
        ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('Normalized $\mu$(E)')
        ax.set_xlim(i[1].energy[0]-2,i[1].energy[-1]+6)  
        
        if export_figure:
            savefig('./plots/'+mpid+'_'+env0[0][0].specie.name+'.png', format='png', dpi=150)  
        
    else: print('XANES is not available. Try \"run_feff=True\"')
    
    return


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# This searches for structures in MP database.
def search_MP(mpr,search_pattern,nmax=None):
    
    if nmax is None: nmax=11
    
    mpid_list = []
    data = mpr.get_data(search_pattern, data_type="vasp", prop="nsites")
    for i in range(len(data)):
        if data[i]['nsites'] <= nmax: mpid_list.append(data[i]['material_id'])
        
    l = len(mpid_list)
    print("""Found %(l)i structures""" % vars())
    
    return mpid_list

# ICSD 0nly
def search_MP_icsd(mpr,search_pattern,nmax=None):
    
    if nmax is None: nmax=11
    
    found0 = []
    data = mpr.get_data(search_pattern, data_type="vasp", prop="nsites")
    for i in range(len(data)):
        if data[i]['nsites'] <= nmax: found0.append(data[i]['material_id'])
    
    mpid_list_icsd = []
    for i in found0:
        data = mpr.get_data(i, data_type="vasp")
        if data[0]['icsd_ids']:
            mpid_list_icsd.append(i)
            
    l = len(mpid_list_icsd)
    print("""Found %(l)i structures""" % vars())            
    
    return mpid_list_icsd


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






