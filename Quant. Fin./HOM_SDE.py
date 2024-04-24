"""
Created on Tue Oct 31 16:24:10 2023

Topic: HOM SDE simulation
Basic outline for Simulate & Fit the following SDE (analytical solution unknown) 
    dX(t) = [a(b-X(t-tau))]dt + [s2.x(t)]dWt 
(with modification in self.drift & self.diffusion of My_SDE , can be generalized for any SDE )

@author: R.Nandi
"""

import numpy as np
import pandas as pd
from bayes_opt import BayesianOptimization
from bayes_opt.util import UtilityFunction
from scipy.stats import norm as normal_dist


def plot_simulationpath(X,Med_line=False,append_title=""):
    if Med_line:
        X.median(axis=1).plot(xlabel='t',ylabel='avg. X(t)',title="simulation path"+append_title)
    else:
        X.plot(legend=False,xlabel='t',ylabel='X(t)',title="simulation path"+append_title)


class My_SDE:
    """
    Given some initial sequence, simulate the stochastic process for some future time points
    size of the time step = t_scale/#time_points
    The resulting stoc. proc. is the solution of the following SDE-
            dX(t) = [a(b-X(t-tau))]dt + [s2.x(t)]dWt

    """
    def __init__(self, a, b, tau, s2, X0):
        self.a = a # drift
        self.b = b # mean-reversion level
        self.tau = round(tau) # time lag
        self.s2 = s2 # volatility
        self.X0 = X0 # initial value sequence
        self.drift = np.vectorize(lambda x : (self.a)*((self.b) - x))
        self.diffusion = np.vectorize(lambda x: (self.s2)*x)

    def simulate(self, t, n_simu=100,t_scale=1):
        initial = len(self.X0)
        t_spread = (max(t)-min(t))
        t_len = len(t)
        X = [pd.DataFrame(np.matmul(np.ones((n_simu, 1)), np.array(self.X0).reshape(1, -1)))]
        X += [pd.DataFrame(np.zeros((n_simu,t_len)))]
        X = pd.concat(X,axis=1,ignore_index=True)
        #   initiates a DataFrame with n_simu rows ,
        #   each row is just a copy of initial value sequence and additionaly contains placeholder for new simulated values
        dWt_unscaled = np.random.normal(size=(n_simu,t_len))
        #   generates n_simu*length(time) obs. from standard normal rv,
        #   these will be scaled later and will be used as sample increment of Brownian motion
        t = np.append(t,np.nan)
        dt = (t[0] - (self.X0.index)[-1])*t_scale/t_spread
        for i in range(initial,initial+t_len):
            i_ = i-initial
            X.iloc[:,i] = X.iloc[:,i-1] + (self.drift(X.iloc[:,i-1-self.tau]))*dt + self.diffusion(X.iloc[:,i-1])*np.sqrt(dt)*dWt_unscaled[:,i_]
            dt = (t[i_+1]-t[i_])*t_scale/t_spread
        X = X.iloc[:,initial-1:]
        X = X.T
        X.index = ['']+list(t[:-1])
        return X

    def obs_vs_fitted(self,obs,t_scale=1):
        initial = len(self.X0)
        t = obs.index
        t_spread = (max(t)-min(t))
        t_len = len(t)
        X = pd.concat([self.X0,obs])
        Fitted = X.iloc[initial-1:].copy()
        t = np.append(t,np.nan)
        dt = (t[0] - (self.X0.index)[-1])*t_scale/t_spread
        for i in range(initial,initial+t_len):
            i_ = i-initial
            Fitted.iloc[i_+1] = X.iloc[i-1] + (self.drift(X.iloc[i-1-self.tau]))*dt
            dt = (t[i_+1]-t[i_])*t_scale/t_spread
        return pd.DataFrame({"observed":X.iloc[initial-1:],
                             "fitted":Fitted,
                             "fit-obs":Fitted-X.iloc[initial-1:]})

    def forecast(self,recent_train_data,test_time,t_scale=1):
        initial = len(recent_train_data)
        t_spread = (max(test_time)-min(test_time))
        t_len = len(test_time)
        X = pd.concat([recent_train_data,pd.Series(index=test_time,dtype=float)])
        test_time = np.append(test_time,np.nan)
        dt = (test_time[0] - (recent_train_data.index)[-1])*t_scale/t_spread
        for i in range(initial,initial+t_len):
            i_ = i-initial
            X.iloc[i] = X.iloc[i-1] + (self.drift(X.iloc[i-1-self.tau]))*dt
            dt = (test_time[i_+1]-test_time[i_])*t_scale/t_spread
        return X.iloc[initial:]




def estimate_parameters(train_data,past_records,bound_a_b_s2=[(-1e+1,1e+1),(-1e+1,1e+1),(1e-10,1e+1)],use_lag=0,n_simu=int(1e+4),init_points=10,n_iter=(20,50)):
    def scaled_loglikelihood_(a,b,s2):
        sde = My_SDE(a, b, use_lag, s2, past_records)
        simulated_data = (sde.simulate(train_data.index,n_simu)).iloc[1:]
        mu_ = np.median(simulated_data,axis=1)
        sigma_ = 1.4826*np.median(np.abs(simulated_data-mu_.reshape(-1,1)),axis=1)
        return np.mean(normal_dist.logpdf(np.array(train_data),loc=mu_,scale=sigma_))
    optimizer = BayesianOptimization(scaled_loglikelihood_,
                                     pbounds={"a":bound_a_b_s2[0],"b":bound_a_b_s2[1],"s2":bound_a_b_s2[2]},verbose=0)
    optimizer.acquisition_function = UtilityFunction(kind="ei") # focus exploration , slower
    optimizer.maximize(init_points=init_points,n_iter=n_iter[0])
    optimizer.acquisition_function = UtilityFunction(kind="poi") # focus refinement , faster
    optimizer.maximize(init_points=init_points,n_iter=n_iter[1])
    return optimizer.max['params']

