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
As mentioned before, this magnetic signal has no explicit function form, which means it would be hard to calculate the gradients of our loss function regarding the parameters $\Delta t_1, \Delta t_2$. In this case, we would like to leverage a heuristic algorithm called Genetic Algorithm(GA) to help us find the best estimated parameters. Here are some reasons why heuristic algorithm is applicable for our senerio:<br>
1) GA algorithm is not a gradient-based optimization method so it could be leveraged to optimize target functions which are not differentiable.<br>
2) The bounds for our parameter vector is known, which means the optimization searching would only happen within a specific solution space.<br>

For more discussions on the results and any further analysis regarding the measuremnt error, please check my paper:<br>
***H. Zhang, M. He, L. Shi. Terahertz Thickness Measurement Based on Stochastic Optimization Algorithm, Spectrosc. Spectral Anal. 40(2020) 3066-3070.***
### Some Adaptive Strategies for GA
To gurantee the diversity in the late iterations, some adaptive strategies might be applied when updating the population.
For mutation probabality:
$$p_{m,i}=p_{m,lower}+(p_{m,upper}-p_{m_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$
For cross probabality:
$$p_{c,i}=p_{c,lower}+(p_{c,upper}-p_{c_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$
## Mathematical Model in the frequency domain
We could also construct the model in frequency domain:
$$E_{fit}(\omega)=E_{r}(\omega)H(\omega)$$
Where $H(\omega)$ is the transfer function of our sample. Within this transfer function, we would like to estimate the reflection index ***n*** and the thickness ***d***. Here ***n*** is a function of frequency, and the Debye/Lorentz model could be leveraged to illustrate this function relationship. 
$$n=\sqrt{Debye/Lorentz model}$$
By IFFT(Inverse Fast Fourier Transformation):
$$E_{fit}(t)=IFFT(E_{fit}(\omega))=IFFT(E_{r}(\omega)H(\omega))$$
Then we could subsititue it into $Fitting Error$ to get the target optimization loss function.
