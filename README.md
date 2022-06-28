# Modeling-For-Magnetic-Waves
Here is the repository for my postgraduation project. The main focus of this project is to build a mathematical model for magnetic waves and find its application in the engineering industry. Our signal is ditributed within a limited range in the time domain and does not have a explicit function form. 
This is what our reference signal looks like:<br>
![Ref_Signal](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
## Mathematical Model in the time domain
We consider our measured signal(the signal that has propagated through our samples) as a linear combination of the reference signal, based on the concept of "Time-of-Flight":<br>
![Time_of_Flight](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
And this is what our measured signal looks like:<br>
![Mea_Signal](/assets/img/philly-magic-garden.jpg "Magic Gardens")<br>
Then the linear combination would be like:
$$E_m(t)=k_1 E_r(t+\Delta t_1) + k_2 E_r(t+\Delta t_2)$$
The loss function we would like to minimize is the sum of squared residuals

