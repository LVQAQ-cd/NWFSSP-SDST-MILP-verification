# NWFSSP-SDST-MILP-verification
NWFSSP-SDST MILP verification


\begin{equation}
	\text{Minimize} \, C_{\max} = C_{n,m}
\end{equation}
\begin{equation}
	\sum_{j=1}^n X_{j,k}=1,\quad k\in\{1,2,\cdots,n\},
\end{equation}
\begin{equation}
\sum_{k=1}^n X_{j,k}=1,\quad j\in\{1,2,\cdots,n\},
\end{equation}
\begin{equation}\color{blue}
C_{1,w} - \sum_{j=1}^nX_{j,1}(p_{j,w}+st_{0,j}^w) \geq 0,  w\in\{1,\cdots,m\}
\end{equation}
\begin{equation}\color{blue}
	\begin{aligned}
		C_{k,w+1}& - C_{k,w} - \sum_{j=1}^{n}X_{j,k}p_{j,w+1} = 0, \\
		&	 k\in \{1,2,\cdots,n\}, w\in\{1,2,\cdots,m-1\}
	\end{aligned}
\end{equation}
\begin{equation}\color{blue}
	\begin{aligned}
		C_{k+1,w} - C_{k,w}& - \sum_{j=1}^nX_{j,k+1} p_{j,w}-\sum_{j=1}^n\sum_{i=1}^nX_{i,k}X_{j,k+1}st_{i,j}^{w}\geq 0 ,\\
		&	 k\in \{1,2,\cdots,n-1\},w\in\{1,2,\cdots,m\}
	\end{aligned}
\end{equation}
\begin{equation}
	C_{k, w}\geq 0,\quad k\in \{1,2,\cdots,n\},w\in\{1,2,\cdots,m\}
\end{equation}
\begin{equation}
	X_{j,k}\in \{0,1\},\quad  j, k\in \{1,2,\cdots,n\}
\end{equation}
