# Modeling-For-Magnetic-Waves
Here is the repository for my postgraduation project. The main focus of this project is to build a mathematical model for magnetic waves and find its application in the engineering industry. Our signal is ditributed within a limited range in the time domain and does not have a explicit function form. 
This is what our reference signal looks like:<br>
![Ref_Signal](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
## Mathematical Model in the time domain
We consider our measured signal(the signal that has propagated through our samples) as a linear combination of reference signals(at different time points), based on the concept of "Time-of-Flight":<br>
![Time_of_Flight](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
And this is what our measured signal looks like:<br>
![Mea_Signal](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
Then the linear combination would be like:
$$E_{fit} (t)=k_1 E_r(t+\Delta t_1) + k_2 E_r(t+\Delta t_2)$$
with parameters $k_1, \Delta t_1, k_2, \Delta t_2$. <br>
The loss function we would like to minimize is the sum of squared residuals:
$$Fitting Error = \sum_{i=1}^n ( E_{mea} (t_i) - E_{fit}(t_i))$$
## Parameter Estimation for the time domain model
As mentioned before, this magnetic signal has no explicit function form, which means it would be hard to calculate the gradients of loss function $FittingError$ regarding the parameters $\Delta t_1, \Delta t_2$. In this case, we would like to leverage a heuristic algorithm called Genetic Algorithm(GA) to help us find the best estimated parameters. Here are some reasons why heuristic algorithm is applicble in our senerio:

### Some Adaptive Strategies for GA
In the 
## Mathematical Model in the frequency domain

