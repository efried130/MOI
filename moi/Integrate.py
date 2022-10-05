#Standard imports
import warnings

# Third-party imports
import numpy as np
from scipy import optimize
from numpy import random

class Integrate:
     """Integrates reach-level FLPE algorithm data.
     Attributes
     ----------
     alg_dict: dict
         dictionary of algorithm data stored by algorithm name as numpy arrays
     basin_dict: dict
         dict of reach_ids and SoS file needed to process entire basin of data
     integ_dict: dict
         dict of integrator estimate data
     moi_params: ??
         ??
     sos_dict: dict
         dictionary of SoS data
     Methods
     -------
     get_pre_mean_q()
         calculate the mean discharge for each reach-level FLPE algorithm
     integrate()
         integrate and store reach-level data
     """

     def __init__(self, alg_dict, basin_dict, sos_dict, sword_dict, obs_dict,VerboseFlag):
          """
          Parameters
          ----------
          alg_dict: dict
               dictionary of algorithm data stored by algorithm name as numpy arrays
          basin_dict: dict
               dict of reach_ids and SoS file needed to process entire basin of data
          sos_dict: dict
               dictionary of SoS data
          """

          self.alg_dict = alg_dict
          self.basin_dict = basin_dict
          self.obs_dict = obs_dict
          self.sword_dict = sword_dict
          self.integ_dict = {
               "pre_q_mean": np.array([]),
               "q_mean": np.array([]),
               "flpe": {
                    "geobam" : np.array([]),
                    "hivdi" : np.array([]),
                    "metroman" : np.array([]),
                    "momma" : np.array([]),
                    "sad" : np.array([])
               }
          }
          self.moi_params = None
          self.sos_dict = sos_dict
          self.VerboseFlag = VerboseFlag

          self.get_pre_mean_q()

     def get_pre_mean_q(self):
          """Calculate the mean discharge for each reach-level FLPE algorithm.
          This represents the mean in time for each algorithm and each reach.
          This should be done prior to integration operations.
          """

          for alg in self.alg_dict:
               for reach in self.alg_dict[alg]:
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         if self.alg_dict[alg][reach]['s1-flpe-exists']:
                              self.alg_dict[alg][reach]['qbar']=np.nanmean(self.alg_dict[alg][reach]['q'])
                              if np.isnan(self.alg_dict[alg][reach]['qbar']):
                                   #print('Calculated a mean flow nan for ',alg,'reachid=',reach,'. Using prior instead.')
                                   self.alg_dict[alg][reach]['qbar']=self.sos_dict[reach]['Qbar']
          


     def pull_sword_attributes_for_reach(self,k):
         """
         Pull out needed SWORD data from the continent dataset arrays for a particular reach    
         """
    
         sword_data_reach={}
         # extract all single-dimension variables, including number of orbits and reach ids needed for multi-dim vars
         for key in self.sword_dict:
             if np.shape(self.sword_dict[key]) == (self.sword_dict['num_reaches'],):
                 sword_data_reach[key]=self.sword_dict[key][k]

         # extract multi-dim vars
         for key in self.sword_dict:            
             if key == 'rch_id_up':
                 sword_data_reach[key]=self.sword_dict[key][0:sword_data_reach['n_rch_up'],k]
             elif key == 'rch_id_dn':
                 sword_data_reach[key]=self.sword_dict[key][0:sword_data_reach['n_rch_down'],k]
             elif key == 'swot_orbits':
                 sword_data_reach[key]=self.sword_dict[key][0:sword_data_reach['swot_obs'],k]
    
         return sword_data_reach

     def ChecksPriorToAddingJunction(self,junction_to_check):
        #check to see if this one already exists
         AlreadyExists=False
         for junction in self.junctions:
             if junction['upflows']==junction_to_check['upflows'] and junction['downflows']==junction_to_check['downflows']: 
                 AlreadyExists=True    
    
         #check to see if all reaches we've identified area in the basin 
         AllReachesInReachFile=True
         for r in junction_to_check['upflows']:
             if str(r) not in self.basin_dict['reach_ids']:
                 AllReachesInReachFile=False
         for r in junction_to_check['downflows']:
             if str(r) not in self.basin_dict['reach_ids']:
                 AllReachesInReachFile=False            
            
         return AlreadyExists,AllReachesInReachFile

     def CreateJunctionList(self):
         # create list of junctions
         self.junctions=list()

         for reach in self.basin_dict['reach_ids']:
             reach=np.int64(reach)
             k=np.argwhere(self.sword_dict['reach_id'] == reach)
             k=k[0,0]
    
             # extract reach dictionary for reach k
             sword_data_reach=self.pull_sword_attributes_for_reach(k)    
    
             #1 try adding the upstream junction
             junction_up=dict()    
             junction_up['originating_reach_id']=reach
    
             #1.1 add the reaches upstream of this junction
             junction_up['upflows']=list()
             for i in range(sword_data_reach['n_rch_up']):
                 junction_up['upflows'].append(sword_data_reach['rch_id_up'][i] )
       
             #1.2 for one of the reaches upstream of this junction, add all their downstream reaches
             junction_up['downflows']=list()
             kup=np.argwhere(self.sword_dict['reach_id'] == junction_up['upflows'][0])
             kup=kup[0,0]
             sword_data_reach_up=self.pull_sword_attributes_for_reach(kup)
             for j in range(sword_data_reach_up['n_rch_down']):
                 junction_up['downflows'].append(sword_data_reach_up['rch_id_dn'][j] )

             AlreadyExists,AllReachesInReachFile=self.ChecksPriorToAddingJunction(junction_up)
    
             if not AlreadyExists and AllReachesInReachFile:
                 self.junctions.append(junction_up)

             #2 try adding the downstream junction
             junction_dn=dict()
             junction_dn['originating_reach_id']=reach #just adding this for bookkeeping/debugging purposes
    
             #2.1 add the reaches downstream of this junction
             junction_dn['downflows']=list()
             for i in range(sword_data_reach['n_rch_down']):
                 junction_dn['downflows'].append(sword_data_reach['rch_id_dn'][i] )            

             #2.2 for one of the reaches downstream of the junction, add all their upstream reaches
             junction_dn['upflows']=list()
             kdn=np.argwhere(self.sword_dict['reach_id'] == junction_dn['downflows'][0])
             kdn=kdn[0,0]
             sword_data_reach_dn=self.pull_sword_attributes_for_reach(kdn)
             for j in range(sword_data_reach_dn['n_rch_up']):
                 junction_dn['upflows'].append(sword_data_reach_dn['rch_id_up'][j] )

             AlreadyExists,AllReachesInReachFile=self.ChecksPriorToAddingJunction(junction_dn)

             if not AlreadyExists and AllReachesInReachFile:
                 self.junctions.append(junction_dn) 

     def RemoveDamReaches(self):
         for reachid in self.basin_dict['reach_ids']:
             k=np.argwhere(self.sword_dict['reach_id']==int(reachid))
             k=k[0,0]
    
             if self.sword_dict['n_rch_down'][k] == 1:
                 # 0. find the reach downstream of the target reach
                 rid_down=str(self.sword_dict['rch_id_dn'][0,k])
        
                 # if the reach downstream of the target reach is a dam...
                 if rid_down[-1] == '4':
            
                     k_down=np.argwhere(self.sword_dict['reach_id']==int(rid_down))
                     k_down=k_down[0,0]
            
                     # 1. find the reach downstream of the reach downstream of the target reach
                     if self.sword_dict['n_rch_down'][k_down] == 1:
                         rid_down_down=str(self.sword_dict['rch_id_dn'][0,k_down])
                         k_down_down=np.argwhere(self.sword_dict['reach_id']==int(rid_down_down))
                         k_down_down=k_down_down[0,0]             
                                
                         if rid_down_down[-1] == '1':
#                             if self.VerboseFlag:
#                                 print('Removing reach ',rid_down)
                             #2. point the target reach at the reach downstream of the reach downstream of the target reach
                             self.sword_dict['rch_id_dn'][0,k] = int(rid_down_down)
                             #3. point the reach downstream of the reach downstream of the target reach, back to the target reach
                             self.sword_dict['rch_id_up'][0,k_down_down]=int(reachid)
                    
                         elif rid_down_down[-1] == '4':
                             # this happens if there are two type 4 reaches in a row, downstream of the target reach
                             if self.sword_dict['n_rch_down'][k_down_down] == 1:
#                                 if self.VerboseFlag:
#                                     print('Removing reach ',rid_down, ' and',rid_down_down)
                    
                                 #1b. find the reach downstream of the reach downstream of the reach downstream of the target reach
                                 rid_down_down_down=str(self.sword_dict['rch_id_dn'][0,k_down_down])
                                 k_down_down_down=np.argwhere(self.sword_dict['reach_id']==int(rid_down_down_down))
                                 k_down_down_down=k_down_down_down[0,0]        
                        
                                 if rid_down_down_down[-1]=='1':
                                     #2b. point the target reach at the reach downstream of the reach downstream of the reach downstream of the target reach
                                     self.sword_dict['rch_id_dn'][0,k] = int(rid_down_down_down)
                                     #3b. point the reach downstream of the reach downstream of the reach downstream of the target reach, back to the target reach
                                     self.sword_dict['rch_id_up'][0,k_down_down_down]=int(reachid)
                    

     def MOI_ObjectiveFunc(self,Q,Qbar,sigmaQ):
         # Q - value of discharge vector at which to evaluate objective function
         # Qbar - prior value of Q e.g. stage 1 McFLI
         # sigmaQ - vector of Q uncertainty
         n=np.size(Q,0)
         C=np.diag( np.reshape(sigmaQ**-1,[n,]) )
         d=Qbar*sigmaQ**(-1)
         res=C@Q-d
         y=np.linalg.norm(res,2)
         return y


     def bam_objfun(self,params,obs,qbar_target): 
          qbam=self.bam_flowlaw(params,obs)
          qbam_bar=np.nanmean(qbam)
          y=abs(qbam_bar-qbar_target)
          return y

     def bam_flowlaw(self,params,obs):
          d_x_area=obs['dA']
          reach_width=obs['w']
          reach_slope=obs['S']
          bam_n=params[0]
          bam_Abar=params[1]
          #keep this equation exactly how this function is written in riverobs 
          qbam = ((d_x_area+bam_Abar)**(5/3) * reach_width**(-2/3) * \
                    (reach_slope)**(1/2)) / bam_n
          qbam=np.reshape(qbam,(1,len(d_x_area)))
          return qbam

     def hivdi_objfun(self,params,obs,qbar_target): 
          q=self.hivdi_flowlaw(params,obs)
          qbar=np.nanmean(q)
          y=abs(qbar-qbar_target)
          return y

     def hivdi_flowlaw(self,params,obs):
          d_x_area=obs['dA']
          reach_width=obs['w']
          reach_slope=obs['S']
          hivdi_alpha=params[0]
          hivdi_beta=params[1]
          hivdi_Abar=params[2]
          #keep this equation exactly how this function is written in riverobs 
          hivdi_n_inv = hivdi_alpha * (
               (d_x_area+hivdi_Abar)/reach_width)**hivdi_beta
          qhivdi = (
               (d_x_area+hivdi_Abar)**(5/3) * reach_width**(-2/3) *
               (reach_slope)**(1/2)) * hivdi_n_inv
          qhivdi=np.reshape(qhivdi,(1,len(d_x_area)))
          return qhivdi

     def metroman_objfun(self,params,obs,qbar_target): 
          q=self.metroman_flowlaw(params,obs)
          qbar=np.nanmean(q)
          y=abs(qbar-qbar_target)
          return y

     def metroman_flowlaw(self,params,obs):
          d_x_area=obs['dA']
          reach_width=obs['w']
          reach_slope=obs['S']
          metro_ninf=params[0]
          metro_p=params[1]
          metro_Abar=params[2]
          #keep this equation exactly how this function is written in riverobs 
          metro_n = metro_ninf * (
               (d_x_area+metro_Abar) / reach_width)**metro_p
          metro_q = (
               (d_x_area+metro_Abar)**(5/3) * reach_width**(-2/3) *
               (reach_slope)**(1/2)) / metro_n

          metro_q=np.reshape(metro_q,(1,len(d_x_area)))
          return metro_q

     def momma_objfun(self,params,obs,qbar_target,aux_var): 
          q=self.momma_flowlaw(params,obs,aux_var)
          qbar=np.nanmean(q)
          y=abs(qbar-qbar_target)

          if np.any(np.isinf(q)):
               # if the search sets H<B then the flow law returns inf
               # in that case return an arbitrary very large number
               y=3e38

          return y

     def momma_flowlaw(self,params,obs,aux_var):
          reach_height=obs['h']
          reach_width=obs['w']
          reach_slope=obs['S']
          momma_B=params[0]
          momma_H=params[1]
          momma_Save=aux_var
          #keep this part exactly how this function is written in riverobs 
          momma_r = 2
          momma_nb = 0.11 * momma_Save**0.18

          momma_q=np.empty( (obs['nt'],)) 

          if momma_H <= momma_B:
               momma_q=np.inf
          else:
               for t in range(obs['nt']):
                    log_factor = np.log10((momma_H-momma_B)/(reach_height[t]-momma_B))

                    if reach_height[t] <= momma_H:
                         momma_n = momma_nb*(1+log_factor)
                         log_check = log_factor > -1
                    else:
                         momma_n = momma_nb*(1-log_factor)
                         log_check = log_factor < 1

                    momma_q[t] = (
                         ((reach_height[t] - momma_B)*(momma_r/(1+momma_r)))**(5/3) *
                         reach_width[t] * reach_slope[t]**(1/2)) / momma_n

               momma_q=np.reshape(momma_q,(1,len(reach_height)))
          return momma_q

     def sad_objfun(self,params,obs,qbar_target): 
          qsad=self.sad_flowlaw(params,obs)
          qsad_bar=np.nanmean(qsad)
          y=abs(qsad_bar-qbar_target)
          return y

     def sad_flowlaw(self,params,obs):
          d_x_area=obs['dA']
          reach_width=obs['w']
          reach_slope=obs['S']
          sad_n=params[0]
          sad_Abar=params[1]
          #keep this equation exactly how this function is written in riverobs 
          qsad = ((d_x_area+sad_Abar)**(5/3) * reach_width**(-2/3) * \
                    (reach_slope)**(1/2)) / sad_n
          qsad=np.reshape(qsad,(1,len(d_x_area)))
          return qsad

     def sic4dvar_objfun(self,params,obs,qbar_target): 
          qsic4dvar=self.sic4dvar_flowlaw(params,obs)
          qsic4dvar_bar=np.nanmean(qsic4dvar)
          y=abs(qsic4dvar_bar-qbar_target)
          return y

     def sic4dvar_flowlaw(self,params,obs):
          d_x_area=obs['dA']
          reach_width=obs['w']
          reach_slope=obs['S']
          sic4dvar_n=params[0]
          sic4dvar_Abar=params[1]
          #keep this equation exactly how this function is written in riverobs 
          qsic4dvar = ((d_x_area+sic4dvar_Abar)**(5/3) * reach_width**(-2/3) * \
                    (reach_slope)**(1/2)) / sic4dvar_n
          qsic4dvar=np.reshape(qsic4dvar,(1,len(d_x_area)))
          return qsic4dvar

     def integrate(self):
          """Integrate reach-level FLPE data."""

          #0 create list of junctions, and figure out problem dimensions
          #0.1 remove type 4 reaches from topology
          self.RemoveDamReaches()
          #0.2 create junction list
          self.CreateJunctionList()       

#          if self.VerboseFlag:
#              print(self.junctions)

          #get sizes of the matrix sizes m & n
          m=0 #number of junctions
          for junction in self.junctions:
              m+=1
              junction['row_num']=m-1    
    
          n=0 #number of reaches
          #for reach in reaches:
          for reach in self.basin_dict['reach_ids']:
              n+=1

          if self.VerboseFlag:
              print('Number of junctions = ',m)
              print('Number of reaches= ',n)

          #1 compute "integrated" discharge. 
          for alg in self.alg_dict:
               print('RUNNING MOI for ',alg)
               Qbar=np.empty([n,])
               i=0
               for reach in self.alg_dict[alg]:
                    Qbar[i]=self.alg_dict[alg][reach]['qbar']
#                    if self.VerboseFlag:
#                        print('Qbar for ',reach,'=',Qbar[i])
                    i+=1

               #define G matrix
               G=np.zeros((m,n))

               for junction in self.junctions:
                   row=junction['row_num']
                   upcols=list()
                   for upflow in junction['upflows']:
                       try:
                            kup=self.basin_dict['reach_ids'].index(str(upflow))
                            upcols.append(kup)
                       except: 
                            print('did not find that one')

                   downcols=list()
                   for downflow in junction['downflows']:
                       try:
                            kdn=self.basin_dict['reach_ids'].index(str(downflow))
                            downcols.append(kdn)
                       except:
                            print('did not find that one')
    
                   for upcol in upcols:
                       G[row,upcol]=1
                   for downcol in downcols:
                       G[row,downcol]=-1

#               if self.VerboseFlag:
#                   print(G)
 
               # solve integrator problem
               bignumber=1e9
     
               sigQ=Qbar*0.4

               cons_massbalance=optimize.LinearConstraint(G,np.zeros(m,),np.zeros(m,))
               cons_positive=optimize.LinearConstraint(np.eye(n),np.zeros(n,),np.ones(n,)*bignumber)

               res=optimize.minimize(fun=self.MOI_ObjectiveFunc,x0=np.reshape(Qbar,[n,]),args=(Qbar,sigQ),method='SLSQP',                      
                      constraints=(cons_massbalance,cons_positive))

               Qintegrator=res.x

               # try computing uncertainty
               if alg == 'metroman_ignore':
                   nEnsemble=20
                   #covQind=sigQ**2*np.eye(n)
                   sigQv=np.reshape(sigQ,(n,1))
                   rho=0.7
                   covQ = np.matmul(sigQv,  sigQv.transpose()) * (rho* np.ones((n,n)) + (np.eye(n)-rho*np.eye(n) )   )  
                   Qens=random.multivariate_normal(Qbar,covQ,nEnsemble)
                   Qmin=10.
                   Qens[Qens<Qmin]=Qmin
               
                   Qensc=np.empty((nEnsemble,n))
                   for i in range(nEnsemble): 
                        res=optimize.minimize(fun=self.MOI_ObjectiveFunc,x0=np.reshape(Qens[i,:],[n,]),args=(Qens[i,:],sigQ),method='SLSQP',                      
                            constraints=(cons_massbalance,cons_positive))
                        Qensc[i,:]=res.x

                   stdQc=Qensc.std(axis=0)
                   stdQc_rel=stdQc/Qintegrator 
               else:
                   stdQc_rel=np.full(n,0.4)

#               if self.VerboseFlag:
#                    print('Integrator Q=',Qintegrator)

               if not res.success:
                    print('Optimization failed for ', alg)
                    if self.VerboseFlag: 
                        print('Qbar=',Qbar)
                    Qintegrator=Qbar

               i=0
               for reach in self.alg_dict[alg]:
                    self.alg_dict[alg][reach]['integrator']={}
                    self.alg_dict[alg][reach]['integrator']['qbar']=Qintegrator[i]
                    self.alg_dict[alg][reach]['integrator']['sbQrel']=stdQc_rel[i]
                    i+=1

          #2 compute optimal parameters for each algorithm's flow law
          #2.1 geobam   
          for reach in self.alg_dict['geobam']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    nflpe=np.nanmean(self.alg_dict['geobam'][reach]['n'])
               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(nflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['geobam'][reach]['n']),np.nanmean(self.alg_dict['geobam'][reach]['a0']))
                    param_bounds=( (0.001,np.inf),(-min(self.obs_dict[reach]['dA'])+1,np.inf))
                    res = optimize.minimize(fun=self.bam_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['geobam'][reach]['integrator']['qbar'] ),
                                        bounds=param_bounds )
                    param_est=res.x

                    #store output
                    self.alg_dict['geobam'][reach]['integrator']['n']=param_est[0]
                    self.alg_dict['geobam'][reach]['integrator']['a0']=param_est[1]
                    self.alg_dict['geobam'][reach]['integrator']['q']=self.bam_flowlaw(param_est,self.obs_dict[reach])
               else: 
                    self.alg_dict['geobam'][reach]['integrator']['n']=np.nan
                    self.alg_dict['geobam'][reach]['integrator']['a0']=np.nan
                    self.alg_dict['geobam'][reach]['integrator']['q']=np.full( (1,self.obs_dict[reach]['nt']),np.nan)

          #2.2 hivdi
          for reach in self.alg_dict['hivdi']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    alphaflpe=np.nanmean(self.alg_dict['hivdi'][reach]['alpha'])

               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(alphaflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['hivdi'][reach]['alpha']), \
                              np.nanmean(self.alg_dict['hivdi'][reach]['beta']),\
                              np.nanmean(self.alg_dict['hivdi'][reach]['a0']))
                    param_bounds=( (0.001,np.inf),(-1e2,1e2),(-min(self.obs_dict[reach]['dA'])+1,np.inf))
                    res = optimize.minimize(fun=self.hivdi_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['hivdi'][reach]['integrator']['qbar'] ),
                                        bounds=param_bounds )
                    param_est=res.x

                    #store output
                    self.alg_dict['hivdi'][reach]['integrator']['alpha']=param_est[0]
                    self.alg_dict['hivdi'][reach]['integrator']['beta']=param_est[1]
                    self.alg_dict['hivdi'][reach]['integrator']['Abar']=param_est[2]
                    self.alg_dict['hivdi'][reach]['integrator']['q']=self.hivdi_flowlaw(param_est,self.obs_dict[reach])
               else: 
                    self.alg_dict['hivdi'][reach]['integrator']['alpha']=np.nan
                    self.alg_dict['hivdi'][reach]['integrator']['beta']=np.nan
                    self.alg_dict['hivdi'][reach]['integrator']['Abar']=np.nan
                    self.alg_dict['hivdi'][reach]['integrator']['q']=np.full( (1,self.obs_dict[reach]['nt']),np.nan)


          #2.3 MetroMan
          for reach in self.alg_dict['metroman']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    naflpe=np.nanmean(self.alg_dict['metroman'][reach]['na'])
               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(naflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['metroman'][reach]['na']), \
                              np.nanmean(self.alg_dict['metroman'][reach]['x1']),\
                              np.nanmean(self.alg_dict['metroman'][reach]['a0']))
                    param_bounds=( (0.001,np.inf),(-1e2,1e2),(-min(self.obs_dict[reach]['dA'])+1,np.inf))
                    res = optimize.minimize(fun=self.metroman_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['metroman'][reach]['integrator']['qbar'] ),
                                        bounds=param_bounds )
                    param_est=res.x

                    #store output
                    self.alg_dict['metroman'][reach]['integrator']['na']=param_est[0]
                    self.alg_dict['metroman'][reach]['integrator']['x1']=param_est[1]
                    self.alg_dict['metroman'][reach]['integrator']['a0']=param_est[2]
                    self.alg_dict['metroman'][reach]['integrator']['q']=self.metroman_flowlaw(param_est,self.obs_dict[reach])
               else: 
                    self.alg_dict['metroman'][reach]['integrator']['na']=np.nan
                    self.alg_dict['metroman'][reach]['integrator']['x1']=np.nan
                    self.alg_dict['metroman'][reach]['integrator']['a0']=np.nan
                    self.alg_dict['metroman'][reach]['integrator']['q']=np.full((1,self.obs_dict[reach]['nt']),np.nan)

          #2.4 MOMMA
          for reach in self.alg_dict['momma']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    Bflpe=np.nanmean(self.alg_dict['momma'][reach]['B'])
               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(Bflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['momma'][reach]['B']), \
                              np.nanmean(self.alg_dict['momma'][reach]['H']))
                              
                    param_bounds=( (0.001,np.min(self.obs_dict[reach]['h'])-0.001),(0.001,np.inf))

                    aux_var=self.alg_dict['momma'][reach]['Save']

                    res = optimize.minimize(fun=self.momma_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['momma'][reach]['integrator']['qbar'],aux_var ),
                                        bounds=param_bounds )
                    param_est=res.x

                    #store output
                    self.alg_dict['momma'][reach]['integrator']['B']=param_est[0]
                    self.alg_dict['momma'][reach]['integrator']['H']=param_est[1]
                    self.alg_dict['momma'][reach]['integrator']['Save']=aux_var
                    self.alg_dict['momma'][reach]['integrator']['q']=self.momma_flowlaw(param_est,self.obs_dict[reach],aux_var)
               else: 
                    self.alg_dict['momma'][reach]['integrator']['B']=np.nan
                    self.alg_dict['momma'][reach]['integrator']['H']=np.nan
                    self.alg_dict['momma'][reach]['integrator']['Save']=np.nan
                    self.alg_dict['momma'][reach]['integrator']['q']=np.full( (1,self.obs_dict[reach]['nt']),np.nan)

          #2.5 SAD
          for reach in self.alg_dict['sad']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    nflpe=np.nanmean(self.alg_dict['sad'][reach]['n'])
               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(nflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['sad'][reach]['n']), \
                              np.nanmean(self.alg_dict['sad'][reach]['a0']))

                    param_bounds=( (0.001,np.inf),(-min(self.obs_dict[reach]['dA'])+1,np.inf))

                    res = optimize.minimize(fun=self.sad_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['sad'][reach]['integrator']['qbar'] ),
                                        bounds=param_bounds )

                    param_est=res.x

                    #store output
                    self.alg_dict['sad'][reach]['integrator']['n']=param_est[0]
                    self.alg_dict['sad'][reach]['integrator']['a0']=param_est[1]
                    self.alg_dict['sad'][reach]['integrator']['q']=self.sad_flowlaw(param_est,self.obs_dict[reach])
               else: 
                    self.alg_dict['sad'][reach]['integrator']['n']=np.nan
                    self.alg_dict['sad'][reach]['integrator']['a0']=np.nan
                    self.alg_dict['sad'][reach]['integrator']['q']=np.full( (1,self.obs_dict[reach]['nt']),np.nan)

          #2.6 SIC4DVar
          for reach in self.alg_dict['sic4dvar']:
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    nflpe=np.nanmean(self.alg_dict['sic4dvar'][reach]['n'])
               if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(nflpe)):
                    with warnings.catch_warnings():
                         warnings.simplefilter("ignore", category=RuntimeWarning)
                         init_params=(np.nanmean(self.alg_dict['sic4dvar'][reach]['n']), \
                              np.nanmean(self.alg_dict['sic4dvar'][reach]['a0']))

                    param_bounds=( (0.001,np.inf),(-min(self.obs_dict[reach]['dA'])+1,np.inf))

                    res = optimize.minimize(fun=self.sic4dvar_objfun,
                                        x0=init_params,
                                        args=(self.obs_dict[reach],self.alg_dict['sic4dvar'][reach]['integrator']['qbar'] ),
                                        bounds=param_bounds )

                    param_est=res.x

                    #store output
                    self.alg_dict['sic4dvar'][reach]['integrator']['n']=param_est[0]
                    self.alg_dict['sic4dvar'][reach]['integrator']['a0']=param_est[1]
                    self.alg_dict['sic4dvar'][reach]['integrator']['q']=self.sic4dvar_flowlaw(param_est,self.obs_dict[reach])
               else: 
                    self.alg_dict['sic4dvar'][reach]['integrator']['n']=np.nan
                    self.alg_dict['sic4dvar'][reach]['integrator']['a0']=np.nan
                    self.alg_dict['sic4dvar'][reach]['integrator']['q']=np.full( (1,self.obs_dict[reach]['nt']),np.nan)
