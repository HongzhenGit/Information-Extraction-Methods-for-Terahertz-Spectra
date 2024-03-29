# Model-based information extraction methods for terahertz spectra
Here is the repository for my master's reearch project. The main focus of this project is to extract information from the spectra of terahertz pulses, and our method was successfully validated in numerical experiments as well as non-destructive applications such as coating layer detection and material characterization. The spectrum of our observed terahertz pulses is ditributed within 0.2-2 THz in the frequency domain. The deviced used is terahertz time-domain spectroscopy.

## Information extraction model in the time domain
We consider our measured signal(the terahertz pulses that has propagated through sample mediums) as a linear combination of reference signals which locate at different time points, based on the concept of "Time-of-Flight":

![Time_of_Flight](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/time_of_flight.png)

Then the time domain model would be like:

$$E_{fit} (t)=k_1 E_r(t+\Delta t_{1}) + k_{2} E_r(t+\Delta t_{2}) \quad with \quad parameters \quad k_{1}, \quad \Delta t_{1}, \quad k_{2}, \quad \Delta t_{2}$$

The loss function we would like to minimize is the sum of squared residuals:

$$Fitting Error = \sum_{i=1}^{n} ( E_{mea} (t_{i}) - E_{fit}(t_{i}))$$

once the minimum value of this loss function is reached, the ToF of each terahertz poulse could be extracted.

### Parameter estimation for the time domain model
The terahertz pulse has no explicit function form, which means it would be hard to calculate the gradients of our loss function regarding its parameters. In this case, we would like to leverage a heuristic algorithm called Genetic Algorithm(GA) to help us find the best estimated parameters. Here are some reasons why heuristic algorithm is applicable for our scenario:
1) GA algorithm is not a gradient-based optimization method so it could be leveraged to optimize target functions which are not differentiable.<br>
2) The bounds for our parameter vector is known, which means the optimization searching would only happen within a specific solution space.<br>

Once the parameters are calibrated, the thickness and refractive index could be calculated according to Time of Flight (TOF) theory and Fresnel’s Equation. For more discussions on the results and any further analysis regarding the measuremnt error, please check my paper:<br>
***[1] H. Zhang, M. He, L. Shi. Terahertz Thickness Measurement Based on Stochastic Optimization Algorithm, Spectrosc. Spectral Anal. 40(2020) 3066-3070.(in Chinese)***

### Some adaptive strategies for GA
To keep the population diversity in the late iterations, some adaptive strategies might be applied when updating the population.<br>
For mutation probabality:

$$p_{m,i}=p_{m,lower}+(p_{m,upper}-p_{m_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$

For cross probabality:

$$p_{c,i}=p_{c,lower}+(p_{c,upper}-p_{c_lower})\frac{f_{i}-f_{min}}{f_{max}-f_{min}}$$

The "lower" and "upper" mean the lower and upper bound of mutation/cross probability. Here are some fiiting results for the Time Domain Model:

![Mea_Signal_Time_Domain](https://github.com/HongzhenGit/Information-Extraction-Methods-for-Terahertz-Spectra/blob/main/Assets/Sample_Signals.png)

## Information extraction model in the frequency domain
We could also construct the model in frequency domain:

$$E_{fit}(\omega)=E_{r}(\omega)H(\omega)$$ 

Where _H(w)_ is the transfer function of our sample. Within this transfer function, we would like to estimate the reflection index _n_ and the thickness _d_. Here _n_ is a function of frequency, and a Debye or Lorentz model could be leveraged to illustrate this function relationship. 

$$n=\sqrt{Debye/Lorentz Model}$$

By IFFT(Inverse Fast Fourier Transformation):

$$E_{fit}(t)=IFFT(E_{fit}(\omega))=IFFT(E_{r}(\omega)H(\omega)) \quad with \quad parameters \quad n, \quad d $$

Then we could subsititue it into _Fitting Error_ to get the target optimization loss function. Some references for Debye/Lorentz model:<br>
***[2] I. Kehuda, S. Khatun, K.J. Reza, M.M. Rahman, M.M. Fakir, Improved debye model for experimental approximation of human breast tissue properties at 6GHz ultra-wideband centre frequency, Int. J. Eng. Technol. 5 (2014) 4708–4717.***<br>
***[3] V.P. Drachev, U.K. Chettiar, A.V. Kildishev, The ag dielectric function in plasmonic metamaterials, Opt. Express 16 (2008) 1186–1195.***

### Parameter estimation for the frequency domain model<br>
Another heuristic algorithm named Differential Evolution(DE) is involved to estimate the parameters of our frequency model. Once the parameters are calibrated, the thickness and refractive index could be extracted simultaneously. We compared the performance of GA and DE, and it was observed that in most cases DE is the better one.<br>
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
***[5] H. Zhang, L. Shi, M. He, Extension of terahertz time-domain spectroscopy: A micron-level thickness gauging technology, Opt. Commun. 506 (2022) 127597.***

Here are some fiiting results for the Frequency Domain Model:

![Mea_Signal Frequency Domain](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/Fitting_results_for_Frequency_Domain_Method.png)

## Discussion about the convergence performance
When appling heuristic algorithms to find out a optimal solution, a large-enough number of iterations would be required to guarantee its convergence. In our case, if the fitting error _< 145_, it will be a converged case, which means, the estimated parameters are as our desiring compared to our benchmark. We counted the number of trials with a larger fitting error than 145 out of 200 trials by different number of iterations, and illustrate their relationship:

![Convergence Performance](https://github.com/HongzhenGit/Modeling-For-Magnetic-Waves/blob/main/Assets/iteration_convergence_performance.png)

There seems to be a Sigmoid-like curve between the number of ourtliers and the number of iterations. We could see that, in our case, at least 90 iterations are required to garantee the convergence of DE algorithm. But there are still about 25 outliers even though the number of iterations is relatively large. The next step of our research, is to enhance the stability and improve the convergence performance of our heuristic algorithm.

## Sparse deconvolution approach for pulse position extraction
Heuristic Optimization might not converge to a fixed point in every signal trail, and they are also time-consuming. To improve the stability and speed of our method, we reformulated our pulse detection problem into a sparse deconvolution problem:

$$\pmb{y}=\pmb{Ah}+\pmb{e}$$

Where _***y***_ is the vector of actually measured signal, _***h***_ is a parameter vector, _***e***_ is the noise term, and matrix _***A***_ is consisted of i-laged referece signal 

$$r^{i}, A=\[r^{0}, r^{1}, ..., r^{n}\]$$

The parameter vector _***h***_ is sparse, becasue the amplitude of multi-reflected pulses attenuates quickly after serveral reflections. To get this sparse parameter vector, a cost function constrined by a L1 norm is applied:

$$\frac{1}{2} \parallel \pmb{Ah}-\pmb{y} \parallel _{2}^{2} + \lambda \parallel \pmb{h} \parallel _{1}$$

and minimized by a LASSO algorithm proposed by M. Tabassum (2018). 

For more information about sparse deconvolution, please refer:<br>
***[6] F. Bobmann, G. Plonka, T. Peter, O. Nemitz and T. Schmitte, Sparse Deconvolution Methods for Ultrasonic NDT, J. Nondestruct Eval. 31 (2012) 225–244.***<br>
***[7] J. Dong, J.B. Jackson, M.Melis, et al. Terahertz frequency-wavelet domain deconvolution for stratigraphic and subsurface investigation of art painting, Optics Express 24(2016) 26972.***<br>
***[8] M. N. Tabassuma and E. Ollila, Sequential adaptive elastic net approach for single-snapshot source localization, J. Acoust. Soc. Am. 143(2018) 3873-3882***<br>
Here are some priliminary results:

![Sparse_Vector_and_Rebuild_Signal](https://github.com/HongzhenGit/Information-Extraction-Methods-for-Terahertz-Spectra/blob/main/Assets/vector_signal.png)

The figure on the left is the sparse parameter vector for the sample on the right. We expected there would be only 2 non-zero peaks, however, there are 5 in the actual parameter vector, and 3 of them locate closely. In the following research, we need to improve the sensitivity of this sparse deconvolution method. 
