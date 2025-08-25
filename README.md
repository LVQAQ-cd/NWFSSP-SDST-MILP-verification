# NWFSSP-SDST-MILP-verification
NWFSSP-SDST MILP verification

$\text{Minimize}\quad C_{\max}=C_{n,m}$

$\sum_{j=1}^n X_{j,k}=1,\quad k\in\{1,2,\cdots,n\}$

$\sum_{k=1}^n X_{j,k}=1,\quad j\in\{1,2,\cdots,n\}$

$\color{blue}{C_{1,w}-\sum_{j=1}^n X_{j,1}\bigl(p_{j,w}+st_{0,j}^{\,w}\bigr)\ge 0},\quad w\in\{1,\cdots,m\}$

$\color{blue}{C_{k,w+1}-C_{k,w}-\sum_{j=1}^{n}X_{j,k}\,p_{j,w+1}=0},\quad k\in\{1,2,\cdots,n\},\ w\in\{1,2,\cdots,m-1\}$

$\color{blue}{C_{k+1,w}-C_{k,w}-\sum_{j=1}^n X_{j,k+1}\,p_{j,w}-\sum_{j=1}^n\sum_{i=1}^n X_{i,k}X_{j,k+1}\,st_{i,j}^{\,w}\ge 0},\quad k\in\{1,2,\cdots,n-1\},\ w\in\{1,2,\cdots,m\}$

$C_{k,w}\ge 0,\quad k\in\{1,2,\cdots,n\},\ w\in\{1,2,\cdots,m\}$

$X_{j,k}\in\{0,1\},\quad j,k\in\{1,2,\cdots,n\}$
