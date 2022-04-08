# Git Repository for ``Computing Longitudinal Moments for Heterogeneous Agent Models''
#### Sergio Ocampo and Baxter Robinson

#### University of Western Ontario

#### **Contact:** socampod@uwo.ca; brobin63@uwo.ca

This repository contains code that implements the histogram iteration method for computing longidudinal moments in heterogeneous agent models.<br/>
The paper describing the method is located in the main folder.<br/>
The method is implemented in Julia and Matlab.<br/>


#### Computing Moments

We take as given the model's solution in the form of policy functions for agents that, together with the stochastic processes of exogenous states, implies an evolution for the agents in the economy. 
This evolution is captured by a Markov kernel, <img src="https://render.githubusercontent.com/render/math?math=T\left(s^'|s\right)>, that maps the transition of agents from a current state $s$ into a future state $s'$ in the state space <img src="https://render.githubusercontent.com/render/math?math={\cal S}>. 
The stationary distribution, $\lambda$, is the solution to<br/>

<img src="https://render.githubusercontent.com/render/math?math=\lambda\left(s^'\right) = \int_{s\in{\cal S}} T\left(s^'|s\right) \lambda\left(s\right) ds,">

We describe how to use $\lambda$ and $T$ to directly compute cross-sectional and longitudinal moments, focusing on the distribution of agents rather than a simulated sample of them. 

# <img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">