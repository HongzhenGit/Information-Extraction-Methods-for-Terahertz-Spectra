# Modeling-For-Magnetic-Waves
Here is the repository for my postgraduation project. The main focus of this project is to build a mathematical model for magnetic waves and find its application in the engineering industry. Our signal is ditributed within a limited range in the time domain and does not have a explicit function form. 
## Mathematical Model in the time domain
We consider our measured signal(the signal that has propagated through our samples) as a linear combination of reference signals(at different time points), based on the concept of "Time-of-Flight":<br>
![Time_of_Flight](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/time_of_flight.png)<br>
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
***[1] H. Zhang, M. He, L. Shi. Terahertz Thickness Measurement Based on Stochastic Optimization Algorithm, Spectrosc. Spectral Anal. 40(2020) 3066-3070.(in Chinese)***
### Some Adaptive Strategies for GA
To gurantee the diversity in the late iterations, some adaptive strategies might be applied when updating the population.
For mutation probabality:
$$p_{m,i}=p_{m,lower}+(p_{m,upper}-p_{m_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$
For cross probabality:
$$p_{c,i}=p_{c,lower}+(p_{c,upper}-p_{c_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$
The "lower" and "upper" mean the lower and upper bound of mutation/cross probability.<br>
Here is some fiiting results for the Time Domain Model:<br>
![Mea_Signal Time Domain](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/Sample%20Signals.png)<br>
## Mathematical Model in the frequency domain
We could also construct the model in frequency domain:
$$E_{fit}(\omega)=E_{r}(\omega)H(\omega)$$
Where $H(\omega)$ is the transfer function of our sample. Within this transfer function, we would like to estimate the reflection index ***n*** and the thickness ***d***. Here ***n*** is a function of frequency, and a Debye or Lorentz model could be leveraged to illustrate this function relationship. 
$$n=\sqrt{Debye/Lorentz Model}$$
By IFFT(Inverse Fast Fourier Transformation):
$$E_{fit}(t)=IFFT(E_{fit}(\omega))=IFFT(E_{r}(\omega)H(\omega))$$
Then we could subsititue it into $Fitting Error$ to get the target optimization loss function.<br>
Some references for Debye/Lorentz model:<br>
***[2] I. Kehuda, S. Khatun, K.J. Reza, M.M. Rahman, M.M. Fakir, Improved debye model for experimental approximation of human breast tissue properties at 6GHz ultra-wideband centre frequency, Int. J. Eng. Technol. 5 (2014) 4708–4717.***<br>
***[3] V.P. Drachev, U.K. Chettiar, A.V. Kildishev, The ag dielectric function in plasmonic metamaterials, Opt. Express 16 (2008) 1186–1195.***
## Parameter Estimation for the frequency domain model
Another heuristic algorithm named Differential Evolution(DE) is involved to estimate the parameters of our frequency model. We compared the performance of GA and DE, and it was observed that in most cases DE is the better one.<br>
Like GA, here are some adaptive strategies designed for DE to gurantee its population diversity during late iterations:<br>
1) DE/Rand/1. Random Selection
$$H_{i}^{k}=V_{p1}^{k}+F(V_{p2}^{k}-V_{p3}^{k})$$
2) DE/Best/1. Best Member and Random Selection
$$H_{i}^{k}=V_{best}^{k}+F(V_{p1}^{k}-V_{p2}^{k})$$
3) DE/Current to Best/1. Current Member, Best Member and Random Selection
$$H_{i}^{k}=V_{i}^{k}+F(V_{best}^{k}-V_{i}^{k})+F(V_{p1}^{k}-V_{p2}^{k})$$

Some references for Differential Algorithm:<br>
***[4] A.K. Qin, V.L. Huang, P.N. Suganthan, Differential evolution algorithm with strategy adaptation for global numerical optimization, IEEE Trans. Evol. Comput. 13 (2009) 398–417.***<br>
For more details and discussions regarding the research above, please check my paper:<br>
***[5] H. Zhang, L. Shi, M. He, Extension of terahertz time-domain spectroscopy: A micron-level thickness gauging technology, Opt. Commun. 506 (2022) 127597.***<br>
Here is some fiiting results for the Frequency Domain Model:<br>
![Mea_Signal Frequency Domain](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/Fitting_results_for_Frequency_Domain_Method.png)<br>
## Discussion about the convergence performance
When appling heuristic algorithms to find out a optimal solution, a large-enough number of iterations would be required to guarantee its convergence. In our case, if the fitting error < 145, it will be a converged case, which means, the estimated parameters are as our desiring compared to our benchmark. We counted the number of trials with a larger fitting error than 145 out of 200 trials by different number of iterations, and illustrate their relationship:<br>
![Convergence Performance](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/iteration_convergence_performance.png)<br>
There seems to be a Sigmoid-like curve between the number of ourtliers and the number of iterations. We could see that, in our case, at least 90 iterations are required to garantee the convergence of DE algorithm. But there are still about 25 outliers even though the number of iterations is relatively large. The next step of our research, is to enhance the stability and improve the convergence performance of our heuristic algorithm.
