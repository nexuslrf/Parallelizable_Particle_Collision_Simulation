https://www.soinside.com/question/FXNseYqxjkZhbrYKHisAgM

To compile code: 
```
gcc -fopenmp -o collision collision.c -lm
```

Formulas can be seen in [MS whiteboard](https://wbd.ms/share/v2/aHR0cHM6Ly93aGl0ZWJvYXJkLm1pY3Jvc29mdC5jb20vYXBpL3YxLjAvd2hpdGVib2FyZHMvcmVkZWVtLzhkMTA5YWQwZGRjMDQxZGZhZjA4MjhlNDY3ZTI3NzRhX0JCQTcxNzYyLTEyRTAtNDJFMS1CMzI0LTVCMTMxRjQyNEUzRA==)

(Tip: use [MathJax Plugin](https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima) with Chrome to render tex math expression)

**Q1**. Given 2 moving particles of size $r$, with $\vec{v_1},\vec{v_2}$, when collision happens?

![](./figs/colli.png)

$\iff$

$\vec{AE}=\lambda\vec{AC},\quad\vec{BF}=\lambda\vec{BD},\quad||\vec{EF}||=2r$

$A(x_1,y_1), C(x_1',y_1'), B(x_2,y_2), D(x_2',y_2')$  $\Rightarrow$

$\vec{EF} = -\vec{AE}+\vec{AB}+\vec{BF}$

$\vec{EF} = \left(\lambda(x_2'-x_1'+x_1-x_2)+x_2-x_1,\lambda(y_2'-y_1'+y_1-y_2)+y_2-y_1\right)$

given $||\vec{EF}||^2=4r^2$, set $\Delta_x=x_2-x_1$ ....; $D_x = \Delta_x'-\Delta_x$ ...

$\Rightarrow$

$\lambda^2(D_x^2+D_y^2)+2\lambda(\Delta_xD_x+\Delta_yD_y)+\Delta_x^2+\Delta_y^2-4r^2=0$

$\Rightarrow$

$\lambda = $ blabla...
