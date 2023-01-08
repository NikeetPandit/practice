# Spectral Analysis Experimentation
In this mini-project, I explore spectral analysis methods; how they excel and their limitations. I investigate both Fourier and least-squares spectral analysis methods, look at the autocovariance functions (and its derivation from power spectrum). I explore different representations, such as investigations in power, amplitude, and percentage variance; thereby allowing me to understand more effectively the most suitable representation given my signal I wish to represent compactly. 
I investigate aliasing, how I can detect aliased peaks (see function),  sub Nyquist effects, different windowing methods, lengths of the series and their effects on the separability close peaks. I highlight here also the benefit of least-squares methods over Fourier in this regard. 
I push Fourier over its strict limitations by introducing trends and gaps and demonstrate how least-squares based spectral analysis are superior in this regard. I use wavelets here and there for effective representations when experimenting with transient signals. 
All series are synthetic which I create (see function) in this mini-project which makes for very effective experimentation. 


<div class="container" style="display: inline-block;">  
  <figure>
  <div style="float: left; padding: 10px;">
    <img src='https://github.com/NikeetPandit/practice/blob/main/Spectral%20Analysis%20Work/functions/IM/read_me_IM.PNG' width="350" height="350" align="center"/>
    <figcaption align="center"><b>Demonstrating sub-Nyquist artefacts; Aliasing even when Nyquist condition is obeyed</b></figcaption>
  </div>

  <div style="float: right; padding: 10px;">
    <img src='https://github.com/NikeetPandit/practice/blob/main/Spectral%20Analysis%20Work/functions/IM/read_me_IM2.PNG' width="350" height="350" align="center"/>
    <figcaption align="center"><b>Creating a synthetic series which has a variable frequency which, due to aliasing, keeps folding over itself</b></figcaption>
  </div>
  </figure>
</div>


### Cite As
Nikeet Pandit (2023). Image Processing Work (https://github.com/NikeetPandit/practice)
* Use functions at own risk
