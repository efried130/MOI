#Standard imports
import warnings

# Third-party imports
import numpy as np
from scipy import optimize

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

     def __init__(self, alg_dict, basin_dict, sos_dict, obs_dict):
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
                         self.alg_dict[alg][reach]['qbar']=np.nanmean(self.alg_dict[alg][reach]['q'])

     def integrate(self):
          """Integrate reach-level FLPE data."""

          nReach=len(self.basin_dict['reach_ids'])

          #1 compute "integrated" discharge. for a stream network this is somewhat complicated and
          #   has not yet been implemented. for a set of several sequential reaches it is just the 
          #   discharge averaged over the reaches. 
          for alg in self.alg_dict:
               Qbar=np.empty([nReach,])
               i=0
               for reach in self.alg_dict[alg]:
                    Qbar[i]=self.alg_dict[alg][reach]['qbar']
                    i+=1
               with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    QbarBar=np.nanmean(Qbar)
               for reach in self.alg_dict[alg]:
                    self.alg_dict[alg][reach]['integrator']={}
                    self.alg_dict[alg][reach]['integrator']['qbar']=QbarBar

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