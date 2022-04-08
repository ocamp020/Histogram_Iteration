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
#### **Computing Moments**


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
#### **Cross-sectional moments.**


These moments involve taking expectations over some variable of interest (x) for some sub-population characterized by states  
  <img src="https://render.githubusercontent.com/render/math?math={\large{s\in S\subseteq{\cal S}}}#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\large\color{white}s\in S\subseteq{\cal S}}#gh-dark-mode-only">,
\begin{align}
    E\left[x|s \in S\right] = \int_{s \in S} x\left(s\right) \lambda_S\left(s\right) ds,
    \label{eq: Cross-sectional moment}
\end{align}
where 
$\lambda_S \equiv \nicefrac{\mathbb{I}_{s \in S}\lambda\left(s\right)}{\int\mathbb{I}_{s \in S}\lambda\left(s\right)ds}$ 
is the marginal distribution of the sub-population in S, with 
$\mathbb{I}_{s \in S}$ 
and indicator variable for whether or not 
$s \in S$.


The equation above applies to a wide range of moments. 
For example, the moments of the wealth distribution (an endogenous state) for the whole population or a subgroup (say among the top income earners), or to define percentiles or other descriptors of the distribution. 
It also applies to moments that depend on both endogenous and exogenous states, for example by making x total income. 
    

As is well understood, these moments can be computed immediately from the solution of the model's stationary distribution 
$\left(\lambda\right)$, 
either by approximating the integral or by calculating the moment from the discrete approximation of the distribution itself.
<br/>
<br/>

---
#### **Longitudinal moments.**
<br/>

Consider an outcome of interest $x\left(s,s^'\right)$ that depends on the initial and final state of an agent.
This outcome could be any function of the initial, final, or intervening states of the agent. 
For example, it could be an indicator function for whether the agent satisfies some condition in the future as being a top earner or having a certain occupation, or the agent's income. 
The expectation of interest depends on whether we focus on the behavior of the group of agents (as in transition rates) or of individual agents (as in the auto-correlation of income). 
In the first case, we must follow the group $(S)$ as a whole and compute
\begin{align}
    E\left[x|s \in S\right] = \int_{s \in S} \, \int_{s^'\in{\cal S}} x\left(s,s^'\right) \lambda^{'}_{S}\left(s^'\right) ds^' \, \lambda_S\left(s\right)ds,
    \label{eq: Longitudinal moment pop}
\end{align}
where $\lambda^{'}_{S}$ is the future distribution of agents conditional on the initial distribution $\lambda_{S}$.
In the second case, we must follow the possible paths of each individual and compute
\begin{align}
    E\left[x|s \in S\right] = \int_{s \in S} \, \int_{s^'\in{\cal S}} x\left(s,s^'\right) \lambda^{'}_{\{s\}}\left(s^'\right) ds^' \, \lambda_S\left(s\right)ds,
    \label{eq: Longitudinal moment ind}
\end{align}
where $\lambda^{'}_{\{s\}}$ is the future distribution of the mass of agents that starts in state $s \in S$ (i.e., given an initial distribution $\delta_{\{s\}}$).

% How to solve integrals 
The difficulty in evaluating the expectations in \eqref{eq: Longitudinal moment pop} and \eqref{eq: Longitudinal moment ind} resides in obtaining $\lambda^{'}_{S}$ and $\lambda^{'}_{\{s\}}$ because this requires accounting for the variation in individual paths between the initial and final period. 
We directly compute $\lambda^{'}_{S}$ and $\lambda^{'}_{\{s\}}$ by iterating forward the initial distribution of agents using the Markov kernel $T$,
\begin{align}
    \lambda^{'}_{S}\left(s^'\right) = \int_{s\in{\cal S}} T\left(s^'|s\right) \lambda_{S}\left(s\right) ds; \quad \lambda^{'}_{\{s\}}\left(s^'\right) = \int_{s\in{\cal S}} T\left(s^'|s\right) \delta_{\{s\}}\left(s\right) ds. 
    \label{eq: Iterating_Distribution} 
\end{align}
Performing the iteration in \eqref{eq: Iterating_Distribution} is relatively costless, as similar iterations are involved in finding the stationary distribution (see \ref{eq: Stationary_Distribution}) and the number of iterations required to compute $\lambda^{'}_{S}$ are finite (and known).\footnote{ 
    The integral in the computation of $\lambda^{'}_{\{s\}}$ is of course superfluous.
    Nevertheless, integration becomes necessary when iterating more than one period into the future, as the initial (degenerate) distribution $\delta_{\{s\}}$ generically distributes mass across the state space ${\cal S}$.
    }
Once $\lambda^{'}_{S}$ and $\lambda^{'}_{\{s\}}$ are obtained, the moments can be computed. % from \eqref{eq: Longitudinal moment pop} and \eqref{eq: Longitudinal moment ind}.