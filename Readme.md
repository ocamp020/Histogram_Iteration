# Git Repository for ``Computing Longitudinal Moments for Heterogeneous Agent Models''
#### **Sergio Ocampo and Baxter Robinson**

#### **University of Western Ontario**

#### **Contact:** socampod@uwo.ca; brobin63@uwo.ca

<br/>
This repository contains code that implements the histogram iteration method for computing longidudinal moments in heterogeneous agent models.<br/>
The paper describing the method is located in the main folder.<br/>
The method is implemented in Julia and Matlab.<br/>
<br/>

#### **Computing Moments**

We take as given the model's solution in the form of policy functions for agents that, together with the stochastic processes of exogenous states, implies an evolution for the agents in the economy. 
This evolution is captured by a Markov kernel, 
  <img src="https://render.githubusercontent.com/render/math?math={
    T\left(s^{'} |s\right)
    }#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}
    T\left(s^{'} |s\right)
    }#gh-dark-mode-only">, 
that maps the transition of agents from a current state $s$ into a future state 
<img src="https://render.githubusercontent.com/render/math?math={
    s^{'}
    }#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}
    s^{'}
    }#gh-dark-mode-only">
in the state space 
  <img src="https://render.githubusercontent.com/render/math?math={
    {\cal S}
    }#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\color{white}
    {\cal S}
    }#gh-dark-mode-only">.
The stationary distribution, &lambda, is the solution to<br/>
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math={\Large
    \lambda\left(s^{'}\right) = \int_{s\in{\cal S}} T\left(s^{'}|s\right) \lambda\left(s\right) ds
    }#gh-light-mode-only">
  <img src="https://render.githubusercontent.com/render/math?math={\Large\color{white}
    \lambda\left(s^{'}\right) = \int_{s\in{\cal S}} T\left(s^{'}|s\right) \lambda\left(s\right) ds
    }#gh-dark-mode-only">
</p>
We describe how to use $\lambda$ and $T$ to directly compute cross-sectional and longitudinal moments, focusing on the distribution of agents rather than a simulated sample of them. 

