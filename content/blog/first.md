---
title: "Weekly Project I"
hero_image: "hero.jpg"
date: 2018-07-24T17:44:36-07:00
description: "Predator Prey models"
---

<h2>Background</h2>
The purpose of these weekly  projects is to act as a way for me to learn new skills, refine skills or just work on something small that I am interested in. As I am not spending a great deal of time on them they may not be the cleanest or most effcient, but I hope that through them I can grow my skills in mathematics, computer science and robotics. 

## The fox and the lynx 
Modelling the population of predators and prey is not a new problem, it can trace its lineage back to 1925 where the model was applied to the populations  of lynx used for fur and the foxes that preyed upon them. The model to describe this was indepenantly developed by Alfred J. Lotka in 1925 and Vito Volterra in 1926. The model is a system of 1st order non-linear ODEs. 
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\LARGE&space;\frac{du}{dt}&space;=&space;\alpha&space;u&space;-&space;\beta&space;uv&space;\\&space;\frac{dv}{dt}&space;=&space;-\gamma&space;v&space;+&space;\sigma&space;\beta&space;uv" title="\LARGE \frac{du}{dt} = \alpha u - \beta uv \\ \frac{dv}{dt} = -\gamma v + \sigma \beta uv" />
</p>
where: 

- *u* is the population of the prey 
- *v* is the population of the predators 

- alpha is the growth of the prey without predators 
- beta is the death rate of the prey due to being eaten
- gamma is the natural death rate of the predators 
- sigma is the amount of prey that must be eaten for a predator to reproduce 

##  Solving the equillibrium points
Before we go on to solve the system of equations numerically, lets do some quick maths to investigate the equillibrium points. At the equillibrium points the rate of growth in population will be 0.
<p align = "center">
	<img src="https://latex.codecogs.com/svg.latex?\large&space;0&space;=&space;\alpha&space;u&space;-&space;\beta&space;uv&space;\\&space;0&space;=&space;-\gamma&space;v&space;&plus;&space;\sigma&space;\beta&space;uv&space;\\&space;\text{Solution&space;1:}\\&space;u&space;=0&space;,&space;v=&space;0&space;\\&space;\text{Solution&space;2:}\\&space;0&space;=&space;\alpha&space;u&space;-&space;\beta&space;uv\\&space;\alpha&space;u&space;=&space;\beta&space;uv\\&space;\alpha&space;=&space;\beta&space;v\\&space;v&space;=&space;\frac{\alpha}{\beta}&space;\\&space;0&space;=&space;-\gamma&space;v&space;&plus;&space;\sigma&space;\beta&space;uv&space;\\&space;\gamma&space;v=\beta&space;\sigma&space;uv&space;\\&space;\gamma&space;=&space;\beta&space;\sigma&space;u&space;\\u&space;=&space;\frac{\gamma}{\beta&space;\sigma}" title="\large 0 = \alpha u - \beta uv \\ 0 = -\gamma v + \sigma \beta uv \\ \text{Solution 1:}\\ u =0 , v= 0 \\ \text{Solution 2:}\\ 0 = \alpha u - \beta uv\\ \alpha u = \beta uv\\ \alpha = \beta v\\ v = \frac{\alpha}{\beta} \\ 0 = -\gamma v + \sigma \beta uv \\ \gamma v=\beta \sigma uv \\ \gamma = \beta \sigma u u = \frac{\gamma}{\beta \sigma}" />
</p>

Now that we have our equillibrium points we are able to linearise the system at these points to create a state space model, this will allow us to asses the stability of the predator prey ecosystem at these two points. Then once we have done this we can solve numerically for the not so linear parts. 

Using my knowledge of control systems the linear system can be represented as a state space. We want an equation that looks like 

X' = AX+B

where X = [u,v]

A is the Jacobian of the system of equations when the values of the equillibrium are substituted in, once this has been done the stability at equillibrium can be assesed. 

<p align="center">
	<img src="https://latex.codecogs.com/svg.latex?\large&space;A&space;=&space;\begin{bmatrix}&space;\frac{\partial&space;f_1}{\partial&space;x_1}&space;&&space;\frac{\partial&space;f_1}{\partial&space;x_2}&space;\\&space;\frac{\partial&space;f_2}{\partial&space;x_1}&space;&&space;\frac{\partial&space;f_2}{\partial&space;x_2}&space;\end{bmatrix}&space;\\&space;A&space;=&space;\begin{bmatrix}&space;\alpha&space;-&space;\beta&space;v&space;&&space;-\beta&space;u&space;\\&space;\beta&space;\sigma&space;v&space;&&space;-\gamma&space;&plus;&space;\beta&space;\sigma&space;u&space;\end{bmatrix}&space;\\&space;\text{at&space;}&space;u,v&space;=&space;0,0&space;\\&space;A&space;=&space;\begin{bmatrix}&space;\alpha&0\\0&-\gamma&space;\end{bmatrix}&space;\\&space;\text{at&space;}u,v&space;=&space;\frac{\gamma}{\beta&space;\sigma},\frac{\alpha}{\beta}&space;\\&space;A&space;=&space;\begin{bmatrix}&space;0&space;&&space;\frac{-\gamma}{\sigma}\\&space;\sigma&space;\alpha&space;&&space;0&space;\end{bmatrix}" title="\large A = \begin{bmatrix} \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} \\ \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} \end{bmatrix} \\ A = \begin{bmatrix} \alpha - \beta v & -\beta u \\ \beta \sigma v & -\gamma + \beta \sigma u \end{bmatrix} \\ \text{at } u,v = 0,0 \\ A = \begin{bmatrix} \alpha&0\\0&-\gamma \end{bmatrix} \\ \text{at }u,v = \frac{\gamma}{\beta \sigma},\frac{\alpha}{\beta} \\ A = \begin{bmatrix} 0 & \frac{-\gamma}{\sigma}\\ \sigma \alpha & 0 \end{bmatrix}" />
</p>




## Solving the system numerically
We can solve the system of equations numerically using numerical methods in python, typically this would be done using scipy and the built in ODE solvers, but I have opted to implement Runge-Kutta 45 with an optional step size optimiser. This is more as a learning opportunity and is probably not as effcient as doing it like a normal person would.

The RKF45 method is implemented in the below class
        if self.t >= self.stop:
            return(False)
        self.calculate_Ks()
        u_next,v_next = self.step4()
        self.u_data.append(self.u)
        self.v_data.append(self.v)
        self.u = u_next
        self.v = v_next
        