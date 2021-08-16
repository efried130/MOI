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
    stage1_estimates: ??
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
                  self.alg_dict[alg][reach]['qbar']=np.nanmean(self.alg_dict[alg][reach]['q'])

    def integrate(self):
        """Integrate reach-level FLPE data.
        
        TODO: Implement

        Parameters
        ----------
        """

        nReach=len(self.basin_dict['reach_ids'])

        for alg in self.alg_dict:
             #1 compute "integrated" discharge. for a stream network this is somewhat complicated. for 
              #  a set of several sequential reaches as in the Sacramento, it is just the average value
             Qbar=np.empty([nReach,])
             i=0
             for reach in self.alg_dict[alg]:
                   Qbar[i]=self.alg_dict[alg][reach]['qbar']
                   i+=1
             QbarBar=np.nanmean(Qbar)
             for reach in self.alg_dict[alg]:
                   self.alg_dict[alg][reach]['integrator']={}
                   self.alg_dict[alg][reach]['integrator']['qbar']=QbarBar

        #2 compute optimal parameters for each algorithm's flow law
        #   just working with BAM to start with
        for reach in self.alg_dict['geobam']:
             nflpe=np.nanmean(self.alg_dict['geobam'][reach]['n'])
             if self.obs_dict[reach]['nt'] > 0 and (not np.isnan(nflpe)):
                  #print('adjusting flow law parameters for reach',reach)
                  #print(self.obs_dict[reach]['nt'])
                  #print(self.alg_dict['geobam'][reach]['a0'])
                  init_params=(np.nanmean(self.alg_dict['geobam'][reach]['n']),np.nanmean(self.alg_dict['geobam'][reach]['a0']))
                  param_bounds=( (0.001,np.inf),(-min(self.obs_dict[reach]['dA'])+1,np.inf))
                  #print(init_params)
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
