# Git Repository for ``Computing Longitudinal Moments for Heterogeneous Agent Models''
#### **Sergio Ocampo and Baxter Robinson**

#### **University of Western Ontario**

#### **Contact:** socampod@uwo.ca; brobin63@uwo.ca

<br/>
This repository contains code that implements the histogram iteration method for computing longidudinal moments in heterogeneous agent models.<br/>
The paper describing the method is located in the main folder.<br/>
The method is implemented in Julia and Matlab.<br/>
<br/>
<br/>

---
### **Computing Moments**


We take as given the model's solution in the form of policy functions for agents that, together with the stochastic processes of exogenous states, implies an evolution for the agents in the economy. 
This evolution is captured by a Markov kernel, 
  <img src="https://render.githubusercontent.com/render/math?math={\large{T\left(s^{'}|s\right)}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}T\left(s^{'}|s\right)}#gh-dark-mode-only">, 
that maps the transition of agents from a current state $s$ into a future state 
  <img src="https://render.githubusercontent.com/render/math?math={\large{s^{'}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}s^{'}}#gh-dark-mode-only">
in the state space 
  <img src="https://render.githubusercontent.com/render/math?math={\large{\cal S}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}{\cal S}}#gh-dark-mode-only">.
The stationary distribution, 
  <img src="https://render.githubusercontent.com/render/math?math={\large\lambda}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda}#gh-dark-mode-only">, 
is the solution to<br/>
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\lambda\left(s^{'}\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\lambda\left(s\right)ds}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}\lambda\left(s^{'}\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\lambda\left(s\right)ds}#gh-dark-mode-only">
</p>

We describe how to use 
  <img src="https://render.githubusercontent.com/render/math?math={\large\lambda}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda}#gh-dark-mode-only"> 
and 
  <img src="https://render.githubusercontent.com/render/math?math={\large{T}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}T}#gh-dark-mode-only"> 
to directly compute cross-sectional and longitudinal moments, focusing on the distribution of agents rather than a simulated sample of them. 
<br/>
<br/>

---
### **Cross-sectional moments.**


These moments involve taking expectations over some variable of interest (x) for some sub-population characterized by states 
  <img src="https://render.githubusercontent.com/render/math?math={\large{s\in S\subseteq{\cal S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}s\in S\subseteq{\cal S}}#gh-dark-mode-only">,
<p align="center">
    <img src="https://render.githubusercontent.com/render/math?math={\Large{E\left[x|s \in S\right] = \int_{s \in S} x\left(s\right) \lambda_S\left(s\right) ds,}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}E\left[x|s \in S\right] = \int_{s \in S} x\left(s\right) \lambda_S\left(s\right) ds,}#gh-dark-mode-only">
</p>
where
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda_S \equiv \frac{\mathbb{I}_{s \in S}\lambda\left(s\right)}{\int\mathbb{I}_{s \in S}\lambda\left(s\right)ds}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda_S \equiv \frac{\mathbb{I}_{s \in S}\lambda\left(s\right)}{\int\mathbb{I}_{s \in S}\lambda\left(s\right)ds}}#gh-dark-mode-only"> 
is the marginal distribution of the sub-population in S, with 
  <img src="https://render.githubusercontent.com/render/math?math={\large{\mathbb{I}_{s \in S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\mathbb{I}_{s \in S}}#gh-dark-mode-only"> 
and indicator variable for whether or not 
  <img src="https://render.githubusercontent.com/render/math?math={\large{s \in S}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}s \in S}#gh-dark-mode-only">.


The equation above applies to a wide range of moments. 
For example, the moments of the wealth distribution (an endogenous state) for the whole population or a subgroup (say among the top income earners), or to define percentiles or other descriptors of the distribution. 
It also applies to moments that depend on both endogenous and exogenous states, for example by making x total income. 
    

As is well understood, these moments can be computed immediately from the solution of the model's stationary distribution 
$\left(\lambda\right)$, 
either by approximating the integral or by calculating the moment from the discrete approximation of the distribution itself.
<br/>
<br/>

---
### **Longitudinal moments.**


Consider an outcome of interest 
  <img src="https://render.githubusercontent.com/render/math?math={\large{x\left(s,s^{'}\right)}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}x\left(s,s^{'}\right)}#gh-dark-mode-only"> 
that depends on the initial and final state of an agent.
This outcome could be any function of the initial, final, or intervening states of the agent. 
For example, it could be an indicator function for whether the agent satisfies some condition in the future as being a top earner or having a certain occupation, or the agent's income. 
The expectation of interest depends on whether we focus on the behavior of the group of agents (as in transition rates) or of individual agents (as in the auto-correlation of income). 
In the first case, we must follow the group 
  <img src="https://render.githubusercontent.com/render/math?math={\large{\left(S\right)}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\left(S\right)}#gh-dark-mode-only"> 
as a whole and compute
<p align="center">
    <img src="https://render.githubusercontent.com/render/math?math={\Large{E\left[x|s\in S\right]=\int_{s\in S}\,\int_{s^{'}\in{\cal S}}x\left(s,s^{'}\right)\lambda^{'}_{S}\left(s^{'}\right)ds^{'}\,\lambda_S\left(s\right)ds,}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}E\left[x|s\in S\right]=\int_{s\in S}\,\int_{s^{'}\in{\cal S}}x\left(s,s^{'}\right)\lambda^{'}_{S}\left(s^{'}\right)ds^{'}\,\lambda_S\left(s\right)ds,}#gh-dark-mode-only">
</p>
where
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda^{'}_{S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda^{'}_{S}}#gh-dark-mode-only"> 
is the future distribution of agents conditional on the initial distribution
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda_{S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda_{S}}#gh-dark-mode-only">.
In the second case, we must follow the possible paths of each individual and compute
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\Large{E\left[x|s\in S\right]=\int_{s\in S}\,\int_{s^{'}\in{\cal S}}x\left(s,s^{'}\right)\lambda^{'}_{\{s\}}\left(s^{'}\right)ds^{'}\,\lambda_S\left(s\right)ds,}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}E\left[x|s\in S\right]=\int_{s\in S}\,\int_{s^{'}\in{\cal S}}x\left(s,s^{'}\right)\lambda^{'}_{\{s\}}\left(s^{'}\right)ds^{'}\,\lambda_S\left(s\right)ds,}#gh-dark-mode-only">
</p>
where
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda^{'}_{\{s\}}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda^{'}_{\{s\}}}#gh-dark-mode-only"> 
is the future distribution of the mass of agents that starts in state
  <img src="https://render.githubusercontent.com/render/math?math={\large{s \in S}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}s \in S}#gh-dark-mode-only">  
(i.e., given an initial distribution 
  <img src="https://render.githubusercontent.com/render/math?math={{\delta_{\{s\}}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}\delta_{\{s\}}}#gh-dark-mode-only">).<br/>
<br/>


The difficulty in evaluating the expectations needed for longitudinal moments resides in obtaining the future distributions
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda^{'}_{S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda^{'}_{S}}#gh-dark-mode-only">
and 
  <img src="https://render.githubusercontent.com/render/math?math={\large{\lambda^{'}_{\{s\}}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\lambda^{'}_{\{s\}}}#gh-dark-mode-only">
because this requires accounting for the variation in individual paths between the initial and final period. 


We directly compute the future distributions by iterating forward the initial distribution of agents using the Markov kernel
  <img src="https://render.githubusercontent.com/render/math?math={\large{\left(T\right)}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}\left(T\right)}#gh-dark-mode-only">,
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\Large{\lambda^{'}_{S}\left(s^'\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\lambda_{S}\left(s\right)ds\,\quad\lambda^{'}_{\{s\}}\left(s^{'}\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\delta_{\{s\}}\left(s\right)ds.}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}\lambda^{'}_{S}\left(s^'\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\lambda_{S}\left(s\right)ds\,\quad\lambda^{'}_{\{s\}}\left(s^{'}\right)=\int_{s\in{\cal S}}T\left(s^{'}|s\right)\delta_{\{s\}}\left(s\right)ds.}#gh-dark-mode-only">
</p>
Performing this iteration is relatively costless, as similar iterations are involved in finding the stationary distribution and the number of iterations required to compute the required future distributions are finite (and known).
Once the future distributions are obtained, the moments can be computed directly from them.