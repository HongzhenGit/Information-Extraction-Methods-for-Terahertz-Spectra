# Information Extraction Methods for Teraheartz Spectra
Here is the repository for my postgraduation project. The main focus of this project is to extract information from the spectra of magnetic waves and find its application in the real engineering industry. Our signal is ditributed within a limited range in the time domain and does not have a explicit function form. 
## Extract information from the time domain spectrum
We consider our measured signal(the signal that has propagated through sample mediums) as a linear combination of reference signals(at different time points), based on the concept of "Time-of-Flight":<br>

![Time_of_Flight](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/time_of_flight.png)<br>
Then the time domain model would be like:

$$E_{fit} (t)=k_1 E_r(t+\Delta t_{1}) + k_{2} E_r(t+\Delta t_{2})$$
<br>
with parameters $k_{1}$, $\Delta t_{1}$, $k_{2}$, $\Delta t_{2}$. The loss function we would like to minimize is the sum of squared residuals:

$$Fitting Error = \sum_{i=1}^{n} ( E_{mea} (t_{i}) - E_{fit}(t_{i}))$$
### Parameter Estimation for the time domain model
As mentioned before, this magnetic signal has no explicit function form, which means it would be hard to calculate the gradients of our loss function regarding the parameters $\Delta t_{1}$, $\Delta t_{2}$.
In this case, we would like to leverage a heuristic algorithm called Genetic Algorithm(GA) to help us find the best estimated parameters. Here are some reasons why heuristic algorithm is applicable for our senerio:<br>
1) GA algorithm is not a gradient-based optimization method so it could be leveraged to optimize target functions which are not differentiable.<br>
2) The bounds for our parameter vector is known, which means the optimization searching would only happen within a specific solution space.<br>
Once the parameters are calibrated, the thickness and refractive index could be calculated according to Time of Flight (TOF) theory and Fresnel’s Equation.<br>
For more discussions on the results and any further analysis regarding the measuremnt error, please check my paper:<br>
***[1] H. Zhang, M. He, L. Shi. Terahertz Thickness Measurement Based on Stochastic Optimization Algorithm, Spectrosc. Spectral Anal. 40(2020) 3066-3070.(in Chinese)***
### Some Adaptive Strategies for GA
To keep the diversity in the late iterations, some adaptive strategies might be applied when updating the population.
For mutation probabality:

$$p_{m,i}=p_{m,lower}+(p_{m,upper}-p_{m_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$

For cross probabality:

$$p_{c,i}=p_{c,lower}+(p_{c,upper}-p_{c_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$

The "lower" and "upper" mean the lower and upper bound of mutation/cross probability. Here are some fiiting results for the Time Domain Model:<br>
![Mea_Signal_Time_Domain](https://github.com/HongzhenGit/Information-Extraction-Methods-for-Terahertz-Spectra/blob/main/Assets/Sample_Signals.png)<br>
## Extract information from the frequency domain spectrum
We could also construct the model in frequency domain:

$$E_{fit}(\omega)=E_{r}(\omega)H(\omega)$$

Where $H(\omega)$ is the transfer function of our sample. Within this transfer function, we would like to estimate the reflection index $n$ and the thickness $d$. Here $n$ is a function of frequency, and a Debye or Lorentz model could be leveraged to illustrate this function relationship. 

$$n=\sqrt{Debye/Lorentz Model}$$

By IFFT(Inverse Fast Fourier Transformation):

$$E_{fit}(t)=IFFT(E_{fit}(\omega))=IFFT(E_{r}(\omega)H(\omega))$$

Then we could subsititue it into $Fitting Error$ to get the target optimization loss function.<br>
Some references for Debye/Lorentz model:<br>
***[2] I. Kehuda, S. Khatun, K.J. Reza, M.M. Rahman, M.M. Fakir, Improved debye model for experimental approximation of human breast tissue properties at 6GHz ultra-wideband centre frequency, Int. J. Eng. Technol. 5 (2014) 4708–4717.***<br>
***[3] V.P. Drachev, U.K. Chettiar, A.V. Kildishev, The ag dielectric function in plasmonic metamaterials, Opt. Express 16 (2008) 1186–1195.***
### Parameter Estimation for the frequency domain model
Another heuristic algorithm named Differential Evolution(DE) is involved to estimate the parameters of our frequency model. Once the parameters are calibrated, the thickness and refractive index could be extracted simultaneously. We compared the performance of GA and DE, and it was observed that in most cases DE is the better one. <br>
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
## Sparse Deconvolution Approach for Pulse Position Extraction
Heuristic Optimization might not converge to a fixed point in every signal trail, and they are also time-consuming. To improve the stability and speed of our method, we reformulated our pulse detection problem into a sparse deconvolution problem:

$$\pmb{y}=\pmb{Ah}+\pmb{e}$$

Where $y$ is the vector of actually measured signal, $h$ is a parameter vector, $e$ is the noise term, and matrix $A$ is consisted of i-laged referece signal $r^{i}$, i.e $A=\[r^{0}, r^{1}, ..., r^{n}\]$. The parameter vector $h$ is sparse, becasue the amplitude of multi-reflected pulses attenuates quickly after serveral reflections. To get this sparse parameter vector, a cost function constrined by a L1 norm is applied:

$$\frac{1}{2} \parallel \pmb{Ah}-\pmb{y} \parallel _{2}^{2} + \lambda \parallel \pmb{h} \parallel _{1}$$

and minimized by a LASSO algorithm proposed by M. Tabassum (2018). 
For more information about sparse deconvolution, please refer:<br>
***[6] F. Bobmann, G. Plonka, T. Peter, O. Nemitz and T. Schmitte, Sparse Deconvolution Methods for Ultrasonic NDT, J. Nondestruct Eval. 31 (2012) 225–244.***<br>
***[7] J. Dong, J.B. Jackson, M.Melis, et al. Terahertz frequency-wavelet domain deconvolution for stratigraphic and subsurface investigation of art painting, Optics Express 24(2016) 26972.***<br>
***[8] M. N. Tabassuma and E. Ollila, Sequential adaptive elastic net approach for single-snapshot source localization, J. Acoust. Soc. Am. 143(2018) 3873-3882***<br>
Here are some priliminary results:<br>
![Sparse_Vector_and_Rebuild_Signal](https://github.com/HongzhenGit/Information-Extraction-Methods-for-Terahertz-Spectra/blob/main/Assets/vector_signal.png)<br>
The figure on the left is the sparse parameter vector for the sample on the right. We expected there would be only 2 non-zero peaks, however, there are 5 in the actual parameter vector, and 3 of them locate closely. In the following research, we need to improve the sensitivity of this sparse deconvolution method. 
