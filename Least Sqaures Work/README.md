# Least Squares Experimentation

In this project, I investigate least-squares and weighted-least squares both from a theoretical and application. I apply my algorithms to fit geometry to GRACE-FO satellite ephemeris. 

While applied for GRACE-FO, routines allow user to fit a geometry (line, plane, curve, ellipse) to any dataset using least-squares and weighted least-squares. The weighted least-squares was very experimental and needs more work. There are also routines for rigorous error analysis and modelling. 

Before this project, I had not investigated the estimator before. Not knowing any better, I thought the least-squares adjustment needed to be performed using Calculus! Hence, in the routines file, at the end of the code you can see 100s of lines of commented code where I attempt to use Python SymPy! It mostly worked but was not totally robust. In any case, it was a great exercise!

--------------------------------------------
![This is an image](https://github.com/NikeetPandit/practice/blob/main/Least%20Sqaures%20Work/IM/read_me_IM.PNG)
*Fit to GRACE-FO Position*

--------------------------------------------

### Cite As
Nikeet Pandit (2023). Image Processing Work (https://github.com/NikeetPandit/practice)
* Use functions at own risk
